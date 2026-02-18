#!/usr/bin/env python3

import argparse
import csv
import gzip
import shutil
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence
import os
import requests
from Bio import SeqIO


ANNOTRIEVE_API = "https://genome.crg.es/annotrieve/api/v0"


@dataclass(frozen=True)
class AnnotationRow:
    annotation_id: str
    assembly_accession: str
    assembly_name: str
    organism_name: str
    taxid: str
    database: str
    provider: str
    source_url: str

    @property
    def species_dirname(self) -> str:
        return f"{self.taxid}_{self.organism_name.replace(' ', '_')}"


def read_annotations_report(tsv_path: Path) -> List[AnnotationRow]:
    rows: List[AnnotationRow] = []
    with tsv_path.open(newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for r in reader:
            rows.append(
                AnnotationRow(
                    annotation_id=r["annotation_id"],
                    assembly_accession=r["assembly_accession"],
                    assembly_name=r.get("assembly_name", ""),
                    organism_name=r["organism_name"],
                    taxid=r["taxid"],
                    database=r.get("database", ""),
                    provider=r.get("provider", ""),
                    source_url=r["source_url"],
                )
            )
    return rows


def annotrieve_genome_url(assembly_accession: str) -> str:
    r = requests.get(
        f"{ANNOTRIEVE_API}/assemblies/{assembly_accession}",
        headers={"Accept": "application/json"},
        timeout=60,
    )
    r.raise_for_status()
    return r.json()["download_url"]


def download_file(url: str, outfile: Path, retries: int = 5, chunk_size: int = 1024 * 1024) -> Path:
    tmp = outfile.with_suffix(outfile.suffix + ".part")
    pos = 0
    headers: Dict[str, str] = {}

    if tmp.exists():
        pos = tmp.stat().st_size
        headers["Range"] = f"bytes={pos}-"

    for attempt in range(1, retries + 1):
        try:
            with requests.get(url, stream=True, headers=headers, timeout=120) as resp:
                resp.raise_for_status()
                mode = "ab" if pos else "wb"
                with tmp.open(mode) as f:
                    for chunk in resp.iter_content(chunk_size=chunk_size):
                        if chunk:
                            f.write(chunk)
            tmp.replace(outfile)
            return outfile
        except Exception:
            if attempt == retries:
                raise

    raise RuntimeError(f"Failed downloading after {retries} attempts: {url}")


_GENE_KEYS = ("Name", "gene", "gene_name", "Alias")
_ID_KEY = "ID"


def _parse_attributes(attr_text: str) -> Dict[str, str]:
    attrs: Dict[str, str] = {}
    for part in attr_text.strip().split(";"):
        if not part:
            continue
        if "=" in part:
            k, v = part.split("=", 1)
        elif " " in part:
            k, v = part.split(" ", 1)
            v = v.strip('"')
        else:
            continue
        attrs[k] = v
    return attrs


def _iter_gff_lines(gff_path: Path) -> Iterable[str]:
    opener = gzip.open if gff_path.suffix.endswith("gz") else open
    with opener(gff_path, "rt") as f:  # type: ignore[arg-type]
        for line in f:
            yield line


def find_gene_feature_ids(gff_path: Path, gene: str) -> List[str]:
    aliases = [a.strip() for a in gene.split("|") if a.strip()]
    alias_set = {a.lower() for a in aliases}
    matched: List[str] = []

    for line in _iter_gff_lines(gff_path):
        if not line or line.startswith("#"):
            continue
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 9:
            continue
        ftype = fields[2]
        if ftype != "gene":
            continue
        attrs = _parse_attributes(fields[8])

        raw_vals = [attrs.get(k, "") for k in _GENE_KEYS]
        tokens: List[str] = []
        for rv in raw_vals:
            if not rv:
                continue
            for t in rv.split(","):
                t = t.strip()
                if t:
                    tokens.append(t)

        if any(t.lower() in alias_set for t in tokens):
            gid = attrs.get(_ID_KEY)
            if gid:
                matched.append(gid)

    return matched


def filter_gff_for_gene(gff_in: Path, gff_out: Path, gene_ids: Sequence[str]) -> None:
    gene_ids_set = set(gene_ids)
    keep_ids = set(gene_ids_set)

    for line in _iter_gff_lines(gff_in):
        if not line or line.startswith("#"):
            continue
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 9:
            continue
        attrs = _parse_attributes(fields[8])
        fid = attrs.get(_ID_KEY)
        parents = attrs.get("Parent", "")
        parent_ids = set(p for p in parents.split(",") if p)

        if fid and fid in keep_ids:
            keep = True
        elif parent_ids & keep_ids:
            keep = True
        elif fields[2] == "mRNA" and parent_ids & gene_ids_set:
            keep = True
        else:
            keep = False

        if keep and fid:
            keep_ids.add(fid)

    with gff_out.open("wt") as out:
        for line in _iter_gff_lines(gff_in):
            if not line:
                continue
            if line.startswith("#"):
                out.write(line)
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            attrs = _parse_attributes(fields[8])
            fid = attrs.get(_ID_KEY)
            parents = attrs.get("Parent", "")
            parent_ids = set(p for p in parents.split(",") if p)
            if (fid and fid in keep_ids) or (parent_ids & keep_ids):
                out.write(line)


def _detect_container_runtime(preferred: Optional[str] = None) -> Optional[str]:
    if preferred:
        if shutil.which(preferred):
            return preferred
        return None

    for candidate in ("docker", "apptainer", "singularity"):
        if shutil.which(candidate):
            return candidate
    return None


def run_gffread_extract(
    gff_path: Path,
    genome_fa_gz: Path,
    out_fa: Path,
    *,
    gffread_container: Optional[str] = None,
    container_runtime: Optional[str] = None,
    host_mount_dir: Optional[Path] = None,
) -> None:
    if gffread_container:
        runtime = _detect_container_runtime(container_runtime)
        if not runtime:
            raise RuntimeError(
                "No container runtime found. Install docker/apptainer/singularity or run with local gffread."
            )

        if host_mount_dir is None:
            raise RuntimeError("Internal error: host_mount_dir is required when using gffread container")

        host_mount_dir = host_mount_dir.resolve()

        def in_container(p: Path) -> str:
            p = p.resolve()
            rel = p.relative_to(host_mount_dir)
            return str(Path("/data") / rel)

        if runtime == "docker":
            cmd = [
                "docker",
                "run",
                "--rm",
                "-v",
                f"{host_mount_dir}:/data",
                "-w",
                "/data",
                gffread_container,
                "gffread",
                "-w",
                in_container(out_fa),
                "-g",
                in_container(genome_fa_gz),
                in_container(gff_path),
            ]
        else:
            # apptainer/singularity have near-identical CLI; both support docker:// URIs
            image_ref = gffread_container
            if not (image_ref.startswith("docker://") or image_ref.endswith(".sif")):
                image_ref = "docker://" + image_ref

            cmd = [
                runtime,
                "exec",
                "-B",
                f"{host_mount_dir}:/data",
                image_ref,
                "gffread",
                "-w",
                in_container(out_fa),
                "-g",
                in_container(genome_fa_gz),
                in_container(gff_path),
            ]

        subprocess.run(cmd, check=True)
        return

    cmd = [
        "gffread",
        "-w",
        str(out_fa),
        "-g",
        str(genome_fa_gz),
        str(gff_path),
    ]
    subprocess.run(cmd, check=True)


def select_longest_transcript(fa_path: Path) -> Optional[SeqIO.SeqRecord]:
    best = None
    best_len = -1
    for rec in SeqIO.parse(str(fa_path), "fasta"):
        l = len(rec.seq)
        if l > best_len:
            best = rec
            best_len = l
    return best


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--annotations-report", required=True)
    parser.add_argument("--gene", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--start", type=int, default=0, help="Start row index (0-based) into annotations report")
    parser.add_argument("--end", type=int, help="End row index (0-based, exclusive) into annotations report")
    parser.add_argument("--max-assemblies", type=int)
    parser.add_argument("--concurrency", type=int, default=2)
    parser.add_argument("--keep-genomes", action="store_true")
    parser.add_argument(
        "--gffread-container",
        help="Container image to run gffread, e.g. quay.io/biocontainers/gffread:0.12.7--h077b44d_6",
    )
    parser.add_argument(
        "--container-runtime",
        choices=["docker", "apptainer", "singularity"],
        help="Container runtime to use when --gffread-container is set",
    )
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    src_report = Path(args.annotations_report)
    (outdir / "_source_annotations_report.tsv").write_text(src_report.read_text())

    rows = read_annotations_report(src_report)
    start = max(0, args.start)
    end = args.end
    if end is not None and end < start:
        raise ValueError("--end must be >= --start")

    rows = rows[start:end]

    if args.max_assemblies is not None:
        rows = rows[: args.max_assemblies]

    annotations_dir = outdir / "annotations"
    genomes_dir = outdir / "genomes"
    work_dir = outdir / "work"
    gene_dir = outdir / "gene"
    annotations_dir.mkdir(parents=True, exist_ok=True)
    genomes_dir.mkdir(parents=True, exist_ok=True)
    work_dir.mkdir(parents=True, exist_ok=True)
    gene_dir.mkdir(parents=True, exist_ok=True)

    script = Path(__file__).resolve().parent / "async_gff_downloader.py"
    subprocess.run(
        [
            "python",
            str(script),
            "--tsv",
            str(src_report),
            "--start",
            str(start),
            "--end",
            str(start + len(rows)),
            "--outdir",
            str(annotations_dir),
            "--concurrency",
            str(args.concurrency),
        ],
        check=True,
    )
    
    out_fa = gene_dir / f"{row.species_dirname}.longest_transcripts.fa"
    selected_records: List[SeqIO.SeqRecord] = []

    for row in rows:
        species_dir = annotations_dir / row.species_dirname
        gffs = list(species_dir.glob("*.gff*"))
        if not gffs:
            continue
        gff_path = gffs[0]

        gene_ids = find_gene_feature_ids(gff_path, args.gene)
        print(f"Gene IDs for {row.species_dirname}: {gene_ids}")
        if not gene_ids:
            continue

        genome_path = genomes_dir / f"{row.assembly_accession}.genome.fa.gz"
        if not genome_path.exists():
            url = annotrieve_genome_url(row.assembly_accession)
            download_file(url, genome_path)
            
        genome_fa_uncompressed = genome_path.with_suffix('')  # Remove .gz
        if not genome_fa_uncompressed.exists():
            import gzip
            import shutil
            with gzip.open(genome_path, 'rb') as f_in:
                with open(genome_fa_uncompressed, 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)

        # Then use the uncompressed file
        with tempfile.TemporaryDirectory(dir=str(work_dir)) as tmpd:
            tmpd_path = Path(tmpd)
            filtered_gff = tmpd_path / "gene.gff3"
            filter_gff_for_gene(gff_path, filtered_gff, gene_ids)

            transcripts_fa = tmpd_path / "transcripts.fa"
            run_gffread_extract(
                filtered_gff,
                genome_fa_uncompressed,  # Use uncompressed file here
                transcripts_fa,
                gffread_container=args.gffread_container,
                container_runtime=args.container_runtime,
                host_mount_dir=outdir,
            )

            best = select_longest_transcript(transcripts_fa)
            if best is None:
                continue

            best.id = f"{row.assembly_accession}|{row.organism_name}|{best.id}"
            best.description = ""
            selected_records.append(best)

        if not args.keep_genomes and genome_path.exists():
            genome_path.unlink()

    with out_fa.open("w") as out_handle:
        SeqIO.write(selected_records, out_handle, "fasta")
        
    os.remove(genome_fa_uncompressed)
    os.remove(genome_path)
    os.remove(genome_fa_uncompressed + '.fai')
    os.remove(gff_path)


if __name__ == "__main__":
    main()
