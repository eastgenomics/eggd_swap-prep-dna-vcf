"""Shared pytest fixtures/path setup for eggd_swap-prep-dna-vcf tests."""

from __future__ import annotations

import sys
from pathlib import Path

RESOURCES_DIR = Path(__file__).resolve().parent.parent / "resources" / "home" / "dnanexus"
sys.path.insert(0, str(RESOURCES_DIR))
