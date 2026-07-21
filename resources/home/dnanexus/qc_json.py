#!/usr/bin/env python3
"""Assemble the prepare-sites funnel/QC report as one small, flat JSON document.

Takes the build-detection result plus a list of ``step=count`` pairs and
writes them out verbatim. Contains no selection/liftover logic of its own —
every count is computed by the caller (usually via ``bcftools view -H | wc -l``)
so this script's only job is producing one consistently-shaped JSON file.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample-id", default=None)
    parser.add_argument("--build-json", required=True, type=Path)
    parser.add_argument(
        "--count",
        action="append",
        default=[],
        metavar="STEP=N",
        help="Repeatable. Site count remaining after a named pipeline step.",
    )
    args = parser.parse_args(argv)

    build_info = json.loads(args.build_json.read_text())

    counts: dict[str, int] = {}
    for item in args.count:
        if "=" not in item:
            print(f"Malformed --count value (expected STEP=N): {item!r}", file=sys.stderr)
            return 2
        step, _, value = item.partition("=")
        if not value.isdigit():
            print(
                f"Malformed --count value (N must be a non-negative integer): {item!r}",
                file=sys.stderr,
            )
            return 2
        counts[step] = int(value)

    report = {
        "schema_version": "prepare-qc.v1",
        "sample_id": args.sample_id,
        "build_detection": build_info,
        "site_counts_by_step": counts,
        "final_site_count": counts.get("final"),
    }
    print(json.dumps(report, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
