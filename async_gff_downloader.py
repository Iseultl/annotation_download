
import asyncio
import aiohttp
import aiofiles
import csv
import os
import random
from pathlib import Path
import ssl

# =====================
# Config
# =====================
CONCURRENCY = 2          # safe for HPC
TIMEOUT = aiohttp.ClientTimeout(total=3600)
MAX_RETRIES = 5
BASE_DELAY = 5

# =====================
# Async download logic
# =====================
async def download_file(session, semaphore, row, base_dir: Path):
    async with semaphore:
        organism = row["organism_name"].replace(" ", "_")
        species_dir = base_dir / f"{row['taxid']}_{organism}"
        species_dir.mkdir(parents=True, exist_ok=True)

        url = row["source_url"]
        filename = species_dir / os.path.basename(url)

        if filename.exists():
            return

        for attempt in range(1, MAX_RETRIES + 1):
            try:
                async with session.get(url) as resp:
                    resp.raise_for_status()
                    async with aiofiles.open(filename, "wb") as f:
                        async for chunk in resp.content.iter_chunked(1024 * 1024):
                            await f.write(chunk)
                break  # success

            except (aiohttp.ClientError, asyncio.TimeoutError) as e:
                if attempt == MAX_RETRIES:
                    print(f"[FAILED] {url}: {e}")
                    return

                delay = BASE_DELAY * (2 ** (attempt - 1)) + random.uniform(0, 2)
                await asyncio.sleep(delay)

        # Write README once
        readme = species_dir / "README.txt"
        if not readme.exists():
            async with aiofiles.open(readme, "w") as f:
                await f.write(
                    f"""Organism name     : {row['organism_name']}
TaxID             : {row['taxid']}
Annotation ID     : {row['annotation_id']}
Assembly accession: {row['assembly_accession']}
Assembly name     : {row['assembly_name']}
Database          : {row['database']}
Provider          : {row['provider']}

Source URL:
{url}
"""
                )

# =====================
# Main async driver
# =====================
async def main(tsv_file, start, end, outdir: str, concurrency: int = CONCURRENCY):
    semaphore = asyncio.Semaphore(concurrency)
    base_dir = Path(outdir)

    with open(tsv_file, newline="") as f:
        reader = list(csv.DictReader(f, delimiter="\t"))
        rows = reader[start:end]

    ssl_context = ssl.create_default_context()
    ssl_context.check_hostname = False
    ssl_context.verify_mode = ssl.CERT_NONE
    
    connector = aiohttp.TCPConnector(ssl=ssl_context)
    
    async with aiohttp.ClientSession(timeout=TIMEOUT, connector=connector) as session:
        tasks = [
            download_file(session, semaphore, row, base_dir)
            for row in rows
        ]
        await asyncio.gather(*tasks)

# =====================
# Entry point
# =====================
if __name__ == "__main__":
    import sys

    if len(sys.argv) == 4:
        tsv = sys.argv[1]
        start = int(sys.argv[2])
        end = int(sys.argv[3])
        outdir = "."
    else:
        import argparse

        parser = argparse.ArgumentParser()
        parser.add_argument("--tsv", required=True)
        parser.add_argument("--start", type=int, default=0)
        parser.add_argument("--end", type=int, required=True)
        parser.add_argument("--outdir", required=True)
        parser.add_argument("--concurrency", type=int, default=CONCURRENCY)
        args = parser.parse_args()

        tsv = args.tsv
        start = args.start
        end = args.end
        outdir = args.outdir
        concurrency = args.concurrency

    asyncio.run(main(tsv, start, end, outdir, concurrency=locals().get("concurrency", CONCURRENCY)))
