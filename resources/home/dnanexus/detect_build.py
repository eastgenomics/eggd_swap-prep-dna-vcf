#!/usr/bin/env python3
"""Detect the genome build of a DNA VCF from its own header, without trusting
a filename or an assumption baked into the pipeline.

Primary signal: the ``##contig=<ID=chr1,length=...>`` (or ``##contig=<ID=1,...>``)
line for the first autosome. Contig *names* do not distinguish hg19/GRCh37 from
GRCh38 (both commonly use ``chr1``), but the declared length does:

    hg19 / GRCh37 chr1 length: 249250621
    GRCh38        chr1 length: 248956422

Secondary, non-authoritative signal: the free-text ``##reference=`` header
line, if present (Illumina TSO500/Pisces VCFs populate this with a filesystem
path such as ``/opt/illumina/resources/genomes/hg19_hardPAR``). This is only
ever used as a cross-check. If it disagrees with the contig-length signal, or
the contig-length signal does not match either known value, this script fails
closed rather than guessing.

The supplied GRCh38 reference's own ``.fai`` is checked the same way, so a
misconfigured "GRCh38" reference input is caught rather than silently
propagated.
"""

from __future__ import annotations

import argparse
import gzip
import json
import re
import sys
from pathlib import Path

HG19_GRCH37_CHR1_LENGTH = 249250621
GRCH38_CHR1_LENGTH = 248956422

_KNOWN_CHR1_LENGTHS = {
    HG19_GRCH37_CHR1_LENGTH: "HG19_OR_GRCH37",
    GRCH38_CHR1_LENGTH: "GRCH38",
}

_CONTIG_LINE = re.compile(r"^##contig=<ID=(?P<id>[^,>]+),length=(?P<length>\d+)")

_REFERENCE_LINE = re.compile(r"^##reference=(?P<value>.+)$")

_BUILD_STRING_HINTS = {
    "HG19_OR_GRCH37": ("hg19", "grch37", "b37", "human_g1k_v37"),
    "GRCH38": ("hg38", "grch38", "b38"),
}


def _open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    return path.open("rt")


def _first_autosome_chr1_length(path: Path) -> int | None:
    """Return the declared length of the first chr1/1 ``##contig`` line."""

    with _open_text(path) as handle:
        for line in handle:
            if not line.startswith("#"):
                break
            match = _CONTIG_LINE.match(line.rstrip("\n"))
            if match and match.group("id") in ("chr1", "1"):
                return int(match.group("length"))
    return None


def _reference_header_value(path: Path) -> str | None:
    with _open_text(path) as handle:
        for line in handle:
            if not line.startswith("#"):
                break
            match = _REFERENCE_LINE.match(line.rstrip("\n"))
            if match:
                return match.group("value")
    return None


def _reference_string_signal(value: str | None) -> str:
    if value is None:
        return "ABSENT"
    lowered = value.lower()
    matches = {
        build
        for build, hints in _BUILD_STRING_HINTS.items()
        if any(hint in lowered for hint in hints)
    }
    if len(matches) == 1:
        return matches.pop()
    if len(matches) > 1:
        return "AMBIGUOUS"
    return "ABSENT"


def _fai_chr1_length(fai_path: Path) -> int | None:
    with fai_path.open("rt") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if fields and fields[0] in ("chr1", "1"):
                return int(fields[1])
    return None


def detect(dna_vcf: Path, reference_fasta_fai: Path) -> dict:
    """Return the build-detection result, or raise ValueError if ambiguous."""

    observed_length = _first_autosome_chr1_length(dna_vcf)
    if observed_length is None:
        raise ValueError(
            "AMBIGUOUS_REFERENCE_BUILD: dna_vcf has no ##contig=<ID=chr1|1,length=...> "
            "header line; genome build cannot be determined from a contig name alone."
        )

    length_signal = _KNOWN_CHR1_LENGTHS.get(observed_length)
    if length_signal is None:
        raise ValueError(
            f"AMBIGUOUS_REFERENCE_BUILD: dna_vcf chr1 length {observed_length} matches "
            f"neither known build (hg19/GRCh37={HG19_GRCH37_CHR1_LENGTH}, "
            f"GRCh38={GRCH38_CHR1_LENGTH})."
        )

    reference_value = _reference_header_value(dna_vcf)
    reference_signal = _reference_string_signal(reference_value)
    if reference_signal not in ("ABSENT", length_signal):
        raise ValueError(
            f"AMBIGUOUS_REFERENCE_BUILD: dna_vcf ##contig length indicates {length_signal} "
            f"but ##reference={reference_value!r} indicates {reference_signal}."
        )

    fai_length = _fai_chr1_length(reference_fasta_fai)
    if fai_length != GRCH38_CHR1_LENGTH:
        raise ValueError(
            f"AMBIGUOUS_REFERENCE_BUILD: supplied reference_fasta_fai chr1 length "
            f"{fai_length} does not match GRCh38 ({GRCH38_CHR1_LENGTH}); refusing to "
            "guess whether liftover is needed against a reference that is not "
            "confirmed GRCh38."
        )

    return {
        "detected_build": length_signal,
        "chr1_length_observed": observed_length,
        "reference_header_value": reference_value,
        "reference_header_signal": reference_signal,
        "reference_fasta_fai_chr1_length": fai_length,
        "liftover_required": length_signal == "HG19_OR_GRCH37",
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("dna_vcf", type=Path)
    parser.add_argument("reference_fasta_fai", type=Path)
    args = parser.parse_args(argv)

    try:
        result = detect(args.dna_vcf, args.reference_fasta_fai)
    except ValueError as error:
        print(str(error), file=sys.stderr)
        return 3

    print(json.dumps(result))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
