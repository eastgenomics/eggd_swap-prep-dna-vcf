#!/usr/bin/env python3
"""Search DNAnexus for a TSO500 DNA VCF partner for each specimen in a list.

For each specimen ID (accession) in the input file, searches every
DNAnexus project whose name starts with "002_" for a file matching
    /output/<run>/eggd_tso500/gather/Results/<sample>/<sample>_MergedSmallVariants.genome.vcf
where <sample> contains the specimen ID.

Queries run in parallel (one thread per specimen) via `dx find data
--all-projects`, then results are filtered locally to projects named
"002_*" and to the expected TSO500 gather/Results folder path -- this is
one broad indexed metadata search per specimen rather than one search per
(specimen, project) pair, which would not scale to the number of 002_
projects in the org.

Usage:
    python3 find_dna_partners.py unmatched_rna_specimen_ids.txt > dna_partners.tsv

Requires `dx` (dxpy) to be installed and already logged in.
"""

from __future__ import annotations

import argparse
import json
import re
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

EXPECTED_FOLDER = re.compile(r"^/output/[^/]+/eggd_tso500/gather/Results/[^/]+$")


def read_specimen_ids(path: Path) -> list[str]:
    ids = [line.strip() for line in path.read_text().splitlines()]
    return [i for i in ids if i and not i.startswith("#")]


def dx_find_data_all_projects(name_pattern: str) -> list[dict]:
    """Run `dx find data --all-projects --name <pattern> --json` and parse it.

    Raises CalledProcessError on a real dx failure; a pattern with zero
    matches still exits 0 with an empty JSON array, so that is not an error.
    """

    result = subprocess.run(
        ["dx", "find", "data", "--all-projects", "--name", name_pattern, "--json"],
        capture_output=True,
        text=True,
        check=True,
    )
    return json.loads(result.stdout)


def search_specimen(specimen_id: str) -> dict:
    """Search for one specimen's TSO500 DNA genome VCF across 002_ projects."""

    name_pattern = f"*{specimen_id}*_MergedSmallVariants.genome.vcf"
    try:
        raw_matches = dx_find_data_all_projects(name_pattern)
    except subprocess.CalledProcessError as error:
        return {
            "specimen_id": specimen_id,
            "status": "DX_QUERY_FAILED",
            "error": error.stderr.strip(),
            "matches": [],
        }

    matches = []
    for item in raw_matches:
        describe = item.get("describe") or {}
        project_id = item.get("project")
        folder = describe.get("folder", "")
        filename = describe.get("name", "")

        matches.append(
            {
                "project_id": project_id,
                "folder": folder,
                "filename": filename,
                "file_id": item.get("id"),
                "folder_matches_expected_pattern": bool(EXPECTED_FOLDER.match(folder)),
            }
        )

    return {"specimen_id": specimen_id, "status": "OK", "matches": matches}


def resolve_project_names(project_ids: set[str]) -> dict[str, str]:
    """Batch-resolve project IDs to names via one `dx describe` per project.

    Kept simple (sequential) since the number of *distinct* projects touched
    is normally much smaller than the number of specimens once real DNA
    projects start recurring across specimens.
    """

    names: dict[str, str] = {}
    for project_id in sorted(project_ids):
        try:
            result = subprocess.run(
                ["dx", "describe", project_id, "--json"],
                capture_output=True,
                text=True,
                check=True,
            )
            names[project_id] = json.loads(result.stdout).get("name", "")
        except subprocess.CalledProcessError:
            names[project_id] = ""
    return names


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("specimen_list", type=Path, help="Plain text file, one specimen ID per line")
    parser.add_argument("--max-workers", type=int, default=8)
    parser.add_argument(
        "--project-prefix",
        default="002_",
        help="Only keep matches from projects whose name starts with this prefix",
    )
    args = parser.parse_args(argv)

    specimen_ids = read_specimen_ids(args.specimen_list)
    print(f"Searching {len(specimen_ids)} specimens across all accessible projects...", file=sys.stderr)

    results: dict[str, dict] = {}
    with ThreadPoolExecutor(max_workers=args.max_workers) as pool:
        futures = {pool.submit(search_specimen, sid): sid for sid in specimen_ids}
        for future in as_completed(futures):
            sid = futures[future]
            results[sid] = future.result()
            print(f"  done: {sid}", file=sys.stderr)

    # Resolve every distinct project ID seen, once, then filter by name prefix.
    all_project_ids = {
        m["project_id"]
        for r in results.values()
        for m in r["matches"]
        if r["status"] == "OK"
    }
    project_names = resolve_project_names(all_project_ids)

    print(
        "\t".join(
            [
                "specimen_id",
                "status",
                "n_matches_total",
                "n_matches_002_and_expected_folder",
                "project_id",
                "project_name",
                "folder",
                "filename",
                "file_id",
            ]
        )
    )
    for sid in specimen_ids:
        r = results[sid]
        if r["status"] != "OK":
            print("\t".join([sid, r["status"], "0", "0", "", "", "", "", ""]))
            continue

        filtered = [
            m
            for m in r["matches"]
            if project_names.get(m["project_id"], "").startswith(args.project_prefix)
            and m["folder_matches_expected_pattern"]
        ]

        if not filtered:
            print(
                "\t".join(
                    [sid, "NO_MATCH_IN_002_PROJECTS", str(len(r["matches"])), "0", "", "", "", "", ""]
                )
            )
            continue

        for m in filtered:
            print(
                "\t".join(
                    [
                        sid,
                        "FOUND",
                        str(len(r["matches"])),
                        str(len(filtered)),
                        m["project_id"],
                        project_names.get(m["project_id"], ""),
                        m["folder"],
                        m["filename"],
                        m["file_id"],
                    ]
                )
            )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
