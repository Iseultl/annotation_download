#!/usr/bin/env python3
"""
Fetch the best annotation and assembly for each taxid listed in taxid.txt.

Candidates are ranked primarily by:
    1. Highest BUSCO completeness
    2. Lowest BUSCO duplication

Further tie-breakers:
    3. Highest single-copy BUSCO percentage
    4. Lowest fragmented percentage
    5. Lowest missing percentage
"""

import requests
import sys
import csv
import logging


logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s"
)

logger = logging.getLogger(__name__)


ANNOTATIONS_API = (
    "https://genome.crg.es/annotrieve/api/v0/annotations"
)

ASSEMBLIES_API = (
    "https://genome.crg.es/annotrieve/api/v0/assemblies"
)

TAXID_FILE = "taxids.txt"
OUTPUT_FILE = "annotations.tsv"


def fetch_json(url):
    """Fetch JSON data from URL with error handling."""

    try:
        response = requests.get(url, timeout=30)
        response.raise_for_status()
        return response.json()

    except requests.exceptions.RequestException as e:
        logger.error(f"Failed to fetch {url}: {e}")
        return None


def load_taxids(filename):
    """Load one taxid per line from a text file."""

    taxids = set()

    try:
        with open(filename, "r") as infile:

            for line in infile:

                taxid = line.strip()

                if taxid:
                    taxids.add(taxid)

    except FileNotFoundError:

        logger.error(
            f"Could not find taxid file: {filename}"
        )

        sys.exit(1)

    logger.info(
        f"Loaded {len(taxids)} taxids"
    )

    return taxids


def busco_ranking(annotation):
    """
    Return a ranking tuple for an annotation.

    Higher tuple values are considered better.
    """

    busco = annotation.get("busco", {})

    complete = float(
        busco.get("complete", 0)
    )

    duplicated = float(
        busco.get("duplicated", 100)
    )

    single_copy = float(
        busco.get("single_copy", 0)
    )

    fragmented = float(
        busco.get("fragmented", 100)
    )

    missing = float(
        busco.get("missing", 100)
    )

    return (
        complete,
        -duplicated,
        single_copy,
        -fragmented,
        -missing
    )


def fetch_best_annotations(selected_taxids):
    """
    Fetch annotation records and retain the best annotation
    for each selected taxid.
    """

    logger.info(
        "Fetching annotations and selecting the best "
        "annotation for each taxid..."
    )

    best_annotations = {}

    offset = 0
    limit = 1000

    while True:

        url = (
            f"{ANNOTATIONS_API}"
            f"?offset={offset}&limit={limit}"
        )

        logger.info(
            f"Fetching annotation page "
            f"offset={offset}, limit={limit}"
        )

        data = fetch_json(url)

        if not data:
            logger.error(
                "Failed to fetch annotations data"
            )
            sys.exit(1)

        results = data.get("results", [])

        if not results:
            break

        for annotation in results:

            taxid = str(
                annotation.get("taxid", "")
            )

            # Skip taxids we are not interested in
            if taxid not in selected_taxids:
                continue

            # Compare against the current best annotation
            if taxid not in best_annotations:

                best_annotations[taxid] = annotation

            else:

                current_best = (
                    best_annotations[taxid]
                )

                if (
                    busco_ranking(annotation)
                    > busco_ranking(current_best)
                ):

                    best_annotations[taxid] = annotation

        total = data.get("total", 0)

        logger.info(
            f"Processed "
            f"{min(offset + len(results), total)}/{total} "
            f"annotations. "
            f"Taxids found so far: "
            f"{len(best_annotations)}/"
            f"{len(selected_taxids)}"
        )

        offset += limit

        if offset >= total:
            break

    logger.info(
        f"Found annotations for "
        f"{len(best_annotations)}/"
        f"{len(selected_taxids)} taxids"
    )

    return best_annotations


def fetch_required_assemblies(required_accessions):
    """
    Fetch assembly records and retain URLs for only
    the required assembly accessions.
    """

    logger.info(
        f"Searching for "
        f"{len(required_accessions)} required assemblies..."
    )

    assemblies_dict = {}

    offset = 0
    limit = 1000

    while True:

        url = (
            f"{ASSEMBLIES_API}"
            f"?offset={offset}&limit={limit}"
        )

        logger.info(
            f"Fetching assembly page "
            f"offset={offset}, limit={limit}"
        )

        data = fetch_json(url)

        if not data:

            logger.error(
                "Failed to fetch assemblies data"
            )

            break

        results = data.get("results", [])

        if not results:
            break

        for assembly in results:

            accession = assembly.get(
                "assembly_accession",
                ""
            )

            if accession in required_accessions:

                download_url = assembly.get(
                    "download_url",
                    ""
                )

                assemblies_dict[
                    accession
                ] = download_url

        total = data.get("total", 0)

        logger.info(
            f"Found {len(assemblies_dict)}/"
            f"{len(required_accessions)} "
            f"required assemblies"
        )

        # Stop early once all assemblies are found
        if (
            len(assemblies_dict)
            == len(required_accessions)
        ):
            logger.info(
                "All required assemblies found."
            )
            break

        offset += limit

        if offset >= total:
            break

    return assemblies_dict


def main():

    # ----------------------------------------
    # 1. Load target taxids
    # ----------------------------------------

    selected_taxids = load_taxids(
        TAXID_FILE
    )

    if not selected_taxids:

        logger.error(
            "No taxids found in taxid.txt"
        )

        sys.exit(1)

    # ----------------------------------------
    # 2. Find the best annotation per taxid
    # ----------------------------------------

    best_annotations = (
        fetch_best_annotations(
            selected_taxids
        )
    )

    if not best_annotations:

        logger.error(
            "No matching annotations found"
        )

        sys.exit(1)

    # ----------------------------------------
    # 3. Extract required assemblies
    # ----------------------------------------

    required_accessions = {

        annotation.get(
            "assembly_accession",
            ""
        )

        for annotation in best_annotations.values()

        if annotation.get(
            "assembly_accession",
            ""
        )
    }

    logger.info(
        f"Need {len(required_accessions)} "
        f"unique assemblies"
    )

    # ----------------------------------------
    # 4. Fetch assembly URLs
    # ----------------------------------------

    assemblies_dict = (
        fetch_required_assemblies(
            required_accessions
        )
    )

    # ----------------------------------------
    # 5. Write output
    # ----------------------------------------

    logger.info(
        f"Writing output to {OUTPUT_FILE}"
    )

    missing_assemblies = set()

    with open(
        OUTPUT_FILE,
        "w",
        newline=""
    ) as tsvfile:

        writer = csv.writer(
            tsvfile,
            delimiter="\t",
            lineterminator="\n"
        )

        # Sort taxids numerically where possible
        for taxid in sorted(
            best_annotations,
            key=int
        ):

            annotation = (
                best_annotations[taxid]
            )

            organism_name = annotation.get(
                "organism_name",
                ""
            )

            annotation_id = annotation.get(
                "annotation_id",
                ""
            )

            annotation_url = (
                annotation
                .get("source_file_info", {})
                .get("url_path", "")
            )

            assembly_accession = annotation.get(
                "assembly_accession",
                ""
            )

            assembly_url = assemblies_dict.get(
                assembly_accession,
                ""
            )

            if (
                assembly_accession
                and not assembly_url
            ):
                missing_assemblies.add(
                    assembly_accession
                )

            busco = annotation.get(
                "busco",
                {}
            )

            writer.writerow([
                taxid,
                organism_name,
                annotation_id,
                annotation_url,
                assembly_accession,
                assembly_url,
                busco.get("complete", ""),
                busco.get("duplicated", ""),
                busco.get("single_copy", ""),
                busco.get("fragmented", ""),
                busco.get("missing", "")
            ])

    # ----------------------------------------
    # 6. Report summary
    # ----------------------------------------

    logger.info(
        f"Successfully wrote "
        f"{len(best_annotations)} records "
        f"to {OUTPUT_FILE}"
    )

    missing_taxids = (
        selected_taxids
        - set(best_annotations.keys())
    )

    if missing_taxids:

        logger.warning(
            f"No annotation found for "
            f"{len(missing_taxids)} taxids"
        )

        logger.warning(
            "Missing taxids:"
        )

        for taxid in sorted(
            missing_taxids,
            key=int
        ):
            logger.warning(
                f"  {taxid}"
            )

    if missing_assemblies:

        logger.warning(
            f"Missing assembly URLs for "
            f"{len(missing_assemblies)} assemblies"
        )

        for accession in sorted(
            missing_assemblies
        ):
            logger.warning(
                f"  {accession}"
            )


if __name__ == "__main__":
    main()