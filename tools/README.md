# DNA/RNA partner-finding tools

These two scripts are **not** part of the DNAnexus app — the app itself
takes only plain file inputs and has no orchestration or pairing logic
(see the top-level README). They are standalone operational tools, run
locally against the `dx` CLI, for the one recurring task the app design
deliberately pushes outside DNAnexus: finding which DNA VCF and RNA BAM in
the organisation actually belong to the same specimen, when the two were
not run together on the same sample sheet.

Both scripts were written and used to support real-world testing of this
app and its sibling `eggd_swap-compare` — see the *Real DNAnexus job
execution* testing sections on both apps' Confluence pages for the results
they produced.

## Why these exist

TSO500 DNA and PanCan RNA sequencing for the same patient specimen are not
always run together, or even in the same sequencing run. When an RNA
sample has no DNA (or vice versa) in the run it actually came from, the
only way to find a candidate partner is to search the whole organisation
by specimen accession. Doing that by hand, one `dx find data` call per
specimen per project, does not scale once there are dozens of specimens and
dozens of `002_*` projects. Both scripts automate exactly that search.

## `find_dna_partners.py` — RNA-led search

Given a list of RNA specimen accessions with no DNA sample in their own
run, searches every DNAnexus project named `002_*` for a TSO500 DNA genome
VCF belonging to the same accession.

- Expected file: `<sample>_MergedSmallVariants.genome.vcf`
- Expected folder shape: `/output/<run>/eggd_tso500/gather/Results/<sample>/`

```bash
python3 tools/find_dna_partners.py unmatched_rna_specimen_ids.txt > dna_partners_query_results.tsv
```

## `find_rna_partners.py` — DNA-led search

The reverse: given a list of DNA specimen accessions with no RNA sample in
their own run, searches every `002_*` project for a PanCan RNA BAM
belonging to the same accession.

- Expected file: `<sample>.mark_duplicates.star.Processed.out.bam`
- Expected folder shape: `/output/<run>/eggd_staraligner-<version>/`

```bash
python3 tools/find_rna_partners.py unmatched_dna_specimen_ids.txt > rna_partners_query_results.tsv
```

## Common behaviour

- **Input**: a plain text file, one specimen accession per line (blank
  lines and lines starting with `#` are skipped).
- **Search strategy**: one broad `dx find data --all-projects --name
  '*<specimen>*<suffix>'` call per specimen, run in parallel across
  specimens (`--max-workers`, default 8) via a thread pool — not one
  search per (specimen, project) pair, which does not scale to the number
  of `002_*` projects in the org.
- **Filtering**: results are first fetched broadly (`--all-projects`),
  then filtered locally to (a) projects whose name starts with `002_`
  (`--project-prefix`, default `002_`) and (b) files sitting in the
  expected pipeline output folder shape for that data type. A raw match
  outside either filter is not reported as a candidate.
- **Project names are resolved once per distinct project ID** seen across
  all specimens (not once per specimen), via a single `dx describe` call
  per project, to avoid redundant API calls when many specimens share the
  same source run.
- **Output**: tab-separated, one row per candidate match (a specimen can
  have zero, one, or several candidates — e.g. the same accession
  appearing in two different real sequencing runs). Status values:
  `FOUND`, `NO_MATCH_IN_002_PROJECTS`, or `DX_QUERY_FAILED` (a genuine `dx`
  CLI error, distinct from a normal zero-match search).

## Known caveat: `002_*` is not a reliable "valid production run" filter

Both scripts restrict candidates to projects named `002_*`, on the
assumption that prefix means a valid, non-failed production run. Real-world
use of these scripts surfaced a case where a failed sequencing run had
never had its `002_` prefix corrected/removed, so a candidate returned by
`find_rna_partners.py` pointed at data from a run that should not have been
treated as valid. **Treat every candidate these scripts return as
provisional** — confirm the source run's actual QC/pass status before
relying on a match, rather than trusting the `002_` prefix alone. This is
an organisational data-hygiene gap, not a defect in either script.

## Requirements

- `dx` (dxpy) installed and already logged in (`dx login`) with read access
  to the projects being searched.
- Python 3.9+ (standard library only — no extra dependencies).
