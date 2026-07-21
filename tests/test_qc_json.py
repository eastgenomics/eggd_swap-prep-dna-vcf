"""Unit tests for resources/home/dnanexus/qc_json.py."""

from __future__ import annotations

import json

import pytest

import qc_json as q


def _write_build_json(tmp_path, **overrides):
    build_info = {
        "detected_build": "GRCH38",
        "chr1_length_observed": 248956422,
        "reference_header_value": "GRCh38",
        "reference_header_signal": "GRCH38",
        "reference_fasta_fai_chr1_length": 248956422,
        "liftover_required": False,
    }
    build_info.update(overrides)
    path = tmp_path / "build.json"
    path.write_text(json.dumps(build_info))
    return path


def test_cli_produces_valid_json_with_expected_shape(tmp_path, capsys):
    build_json = _write_build_json(tmp_path)

    exit_code = q.main(
        [
            "--sample-id",
            "sample-001",
            "--build-json",
            str(build_json),
            "--count",
            "raw=100",
            "--count",
            "final=42",
        ]
    )

    assert exit_code == 0
    output = json.loads(capsys.readouterr().out)

    assert output["schema_version"] == "prepare-qc.v1"
    assert output["sample_id"] == "sample-001"
    assert output["build_detection"]["detected_build"] == "GRCH38"
    assert output["site_counts_by_step"] == {"raw": 100, "final": 42}
    assert output["final_site_count"] == 42


def test_count_without_equals_sign_exits_cleanly(tmp_path, capsys):
    build_json = _write_build_json(tmp_path)

    exit_code = q.main(
        [
            "--build-json",
            str(build_json),
            "--count",
            "raw_no_equals",
        ]
    )

    assert exit_code != 0
    assert "Malformed --count value" in capsys.readouterr().err


def test_count_with_non_integer_value_exits_cleanly(tmp_path, capsys):
    build_json = _write_build_json(tmp_path)

    exit_code = q.main(
        [
            "--build-json",
            str(build_json),
            "--count",
            "raw=not-a-number",
        ]
    )

    assert exit_code != 0
    assert "Malformed --count value" in capsys.readouterr().err
