"""Unit tests for resources/home/dnanexus/detect_build.py."""

from __future__ import annotations

import pytest

import detect_build as d


def _write_vcf(tmp_path, name, header_lines):
    path = tmp_path / name
    path.write_text("\n".join(header_lines) + "\n")
    return path


def _write_fai(tmp_path, name, chr1_length):
    path = tmp_path / name
    path.write_text(f"chr1\t{chr1_length}\t6\t60\t61\n")
    return path


def _hg19_header():
    return [
        "##fileformat=VCFv4.1",
        "##reference=/opt/illumina/resources/genomes/hg19_hardPAR",
        "##contig=<ID=chr1,length=249250621>",
        "##contig=<ID=chr2,length=243199373>",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
    ]


def _grch38_header(reorder_attrs=False):
    if reorder_attrs:
        contig_line = "##contig=<ID=chr1,assembly=GRCh38,length=248956422>"
    else:
        contig_line = "##contig=<ID=chr1,length=248956422,assembly=GRCh38>"
    return [
        "##fileformat=VCFv4.2",
        "##reference=GRCh38",
        contig_line,
        "##contig=<ID=chr2,length=242193529,assembly=GRCh38>",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
    ]


def test_hg19_header_detects_hg19_and_requires_liftover(tmp_path):
    vcf = _write_vcf(tmp_path, "hg19.vcf", _hg19_header())
    fai = _write_fai(tmp_path, "ref.fa.fai", d.GRCH38_CHR1_LENGTH)

    result = d.detect(vcf, fai)

    assert result["detected_build"] == "HG19_OR_GRCH37"
    assert result["liftover_required"] is True


def test_grch38_header_bare_contig_detects_grch38_no_liftover(tmp_path):
    vcf = _write_vcf(tmp_path, "grch38.vcf", _grch38_header())
    fai = _write_fai(tmp_path, "ref.fa.fai", d.GRCH38_CHR1_LENGTH)

    result = d.detect(vcf, fai)

    assert result["detected_build"] == "GRCH38"
    assert result["liftover_required"] is False


def test_contradictory_header_raises_ambiguous(tmp_path):
    header = [
        "##fileformat=VCFv4.1",
        "##reference=GRCh38",  # contradicts the hg19-length contig line below
        "##contig=<ID=chr1,length=249250621>",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
    ]
    vcf = _write_vcf(tmp_path, "contradiction.vcf", header)
    fai = _write_fai(tmp_path, "ref.fa.fai", d.GRCH38_CHR1_LENGTH)

    with pytest.raises(ValueError, match="AMBIGUOUS_REFERENCE_BUILD"):
        d.detect(vcf, fai)


def test_no_contig_line_raises_value_error(tmp_path):
    header = [
        "##fileformat=VCFv4.1",
        "##reference=hg19",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
    ]
    vcf = _write_vcf(tmp_path, "no_contig.vcf", header)
    fai = _write_fai(tmp_path, "ref.fa.fai", d.GRCH38_CHR1_LENGTH)

    with pytest.raises(ValueError, match="AMBIGUOUS_REFERENCE_BUILD"):
        d.detect(vcf, fai)


def test_reordered_contig_attributes_still_detected(tmp_path):
    vcf = _write_vcf(tmp_path, "reordered.vcf", _grch38_header(reorder_attrs=True))
    fai = _write_fai(tmp_path, "ref.fa.fai", d.GRCH38_CHR1_LENGTH)

    result = d.detect(vcf, fai)

    assert result["detected_build"] == "GRCH38"
    assert result["liftover_required"] is False


def test_mismatched_reference_fai_raises_value_error(tmp_path):
    vcf = _write_vcf(tmp_path, "grch38.vcf", _grch38_header())
    fai = _write_fai(tmp_path, "bad_ref.fa.fai", d.HG19_GRCH37_CHR1_LENGTH)

    with pytest.raises(ValueError, match="AMBIGUOUS_REFERENCE_BUILD"):
        d.detect(vcf, fai)
