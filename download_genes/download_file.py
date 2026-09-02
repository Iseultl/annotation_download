#!/usr/bin/env python3

import argparse
import csv
import logging
import sys
import time
import urllib.error
import urllib.request
from datetime import datetime
from pathlib import Path


LOGGER = logging.getLogger("download_file")
RETRY_HEADER = ["taxid", "artifact", "url", "attempt", "timestamp", "error"]


def configure_logging():
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")


def download_file(url, dest_path, max_attempts=3, retry_delay=15):
    """Download *url* to *dest_path*. Returns (ok: bool, error_message: str)."""
    last_error = "NA"
    for attempt in range(1, max_attempts + 1):
        LOGGER.info("Downloading %s -> %s (attempt %s/%s)", url, dest_path, attempt, max_attempts)
        try:
            urllib.request.urlretrieve(url, dest_path)
            LOGGER.info("Download complete: %s", dest_path)
            return True, "NA"
        except urllib.error.HTTPError as e:
            last_error = f"HTTP {e.code} downloading {url}: {e.reason}"
        except urllib.error.URLError as e:
            last_error = f"URL error downloading {url}: {e.reason}"
        except Exception as e:
            last_error = f"Unexpected error downloading {url}: {e}"

        LOGGER.warning("%s", last_error)
        if attempt < max_attempts:
            sleep_seconds = retry_delay * attempt
            LOGGER.info("Sleeping %ss before retry", sleep_seconds)
            time.sleep(sleep_seconds)

    return False, last_error


def append_retry_row(log_tsv, taxid, artifact, url, attempt, error):
    log_tsv = Path(log_tsv)
    log_tsv.parent.mkdir(parents=True, exist_ok=True)
    exists = log_tsv.exists()
    with open(log_tsv, "a", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        if not exists:
            writer.writerow(RETRY_HEADER)
        writer.writerow([taxid, artifact, url, attempt, datetime.now().isoformat(timespec="seconds"), error])


def parse_args(argv):
    parser = argparse.ArgumentParser(description="Download and stage annotation and fasta files")
    parser.add_argument("--taxid", required=True)
    parser.add_argument("--annotation-url", required=True)
    parser.add_argument("--fasta-url", required=True)
    parser.add_argument("--retry-log", default="download_retry.tsv")
    parser.add_argument("--download-delay", type=int, default=8)
    parser.add_argument("--retry-delay", type=int, default=15)
    parser.add_argument("--max-attempts", type=int, default=3)
    return parser.parse_args(argv)


def main(argv=None):
    configure_logging()
    args = parse_args(argv or sys.argv[1:])

    annotation_file = Path("annotation.gff.gz")
    fasta_file = Path("annotation.fasta.gz")

    LOGGER.info("Starting download for taxid %s", args.taxid)

    ok, err = download_file(args.annotation_url, annotation_file, args.max_attempts, args.retry_delay)
    if not ok:
        append_retry_row(args.retry_log, args.taxid, "annotation", args.annotation_url, args.max_attempts, err)
        return 1

    if args.download_delay > 0:
        LOGGER.info("Sleeping %ss before FASTA download", args.download_delay)
        time.sleep(args.download_delay)

    ok, err = download_file(args.fasta_url, fasta_file, args.max_attempts, args.retry_delay)
    if not ok:
        append_retry_row(args.retry_log, args.taxid, "fasta", args.fasta_url, args.max_attempts, err)
        return 1

    return 0


if __name__ == "__main__":
    raise SystemExit(main())