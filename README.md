# eggd_swap-prep-dna-vcf

DNAnexus app: selects germline heterozygous SNVs from a TSO500 DNA genome
VCF, lifts them over to GRCh38 only if the input VCF actually needs it, and
restricts the result to the RNA capture region — producing a reusable site
list for the RNA-swap identity-QC `compare` app (`ASEReadCounter` against a
matched PanCan RNA BAM).

This app deliberately has **no manifest, no candidate ID, no validation-mode
flag, and no coordinator dependency**. Every input is a plain file; every
output is a plain file. Orchestration (deciding which DNA/RNA pairs to run,
running a cross-BAM sweep on a suspected swap) happens outside this app.

See the design rationale and the reboot this app is part of:
<https://cuhbioinformatics.atlassian.net/wiki/spaces/URA/pages/4753260639>

## Pipeline

0. **Detect genome build** from `dna_vcf`'s own header (`##contig=<ID=chr1|1,length=...>`
   contig length — the only reliable signal, since contig *names* don't
   distinguish hg19/GRCh37 from GRCh38). Cross-checked against the VCF's
   `##reference=` line (if present) and against `reference_fasta_fai`'s own
   chr1 length. Fails closed (non-zero exit, no output) on any disagreement
   or unrecognised length — never guesses.
1. Select PASS, biallelic SNV, heterozygous, `FORMAT/DP>=30` candidates.
2. Restrict to the diploid-het VAF band `0.35–0.65`.
3. Normalise contig naming to `chr`-prefixed (idempotent if already
   `chr`-prefixed; a plain rename, run unconditionally regardless of build).
4. If (and only if) Step 0 detected hg19/GRCh37: liftover to GRCh38 via GATK
   `LiftoverVcf` (`RECOVER_SWAPPED_REF_ALT true`, `-XX:MaxRAMPercentage=70.0`).
   If already GRCh38: skip this step; `rejected_liftover_vcf` is still
   produced but header-only.
5. Exclude MNV-type records (`strlen(REF)>1 && strlen(REF)==strlen(ALT)`) —
   these can crash `ASEReadCounter` downstream with an overlapping
   `VariantContext` error.
6. Restrict to the RNA capture region (`bcftools view -R capture_bed`).
7. `bgzip` + `tabix` the final site VCF.

Every step count is recorded in `prepare_qc_json` (a funnel report) and
also echoed as a single summary line to the job log, along with the Step 0
build-detection evidence — so the funnel is visible without downloading the
JSON output.

## Inputs

| Input | Class | Notes |
|---|---|---|
| `dna_vcf` | file | Raw TSO500 genome VCF (PASS + non-PASS calls) |
| `reference_fasta` / `_fai` / `reference_dict` | file | GRCh38 reference matching the RNA aligner (e.g. CTAT) |
| `chain_file` | file | hg19→GRCh38 chain; used only if needed |
| `capture_bed` | file | GRCh38 RNA capture manifest |
| `gatk_jar` | file | GATK jar, supplied directly (no asset-selection indirection) |
| `sample_id` | string, optional | Label only — never used for pairing logic |

`bcftools`/`samtools`/`bgzip`/`tabix` come from the team's formal, **Approved**
`htslib_suite_asset` v1.22 DNAnexus AssetBundle (`record-J1YBvy049yKpP7kk1j4ggqxZ`,
approved [DI-2083](https://cuhbioinformatics.atlassian.net/browse/DI-2083)),
declared under `assetDepends` rather than pulled from apt — not a project-specific
choice. **Caveat:** the asset itself was built in an Ubuntu 20.04 context, and
the RNA-swap-specific qualification record for it
([draft](https://cuhbioinformatics.atlassian.net/wiki/spaces/DV/pages/4751786032))
states the required Ubuntu 24.04 smoke test has not yet been performed — this
app runs on Ubuntu 24.04, so that smoke test should be completed (or the
existing one from the asset's own approval page treated as sufficient
precedent) before treating this dependency as release-qualified.

## Outputs

| Output | Class |
|---|---|
| `site_vcf` / `site_vcf_tbi` | file |
| `prepare_qc_json` | file |

`LiftoverVcf`'s `--REJECT` output (and its per-contig rejection counts) is
logged to the job log and summarised as a count in `prepare_qc_json`, not
persisted as a declared file output — nothing downstream consumes the
rejected records themselves, and in every real run to date the count has
been zero. Keeping it as a log line rather than a file avoids the
header-only-placeholder branch that would otherwise be needed when
liftover is skipped, and matches how other apps in this project surface
low-value diagnostic counts.

## Known limitation, tracked for the `compare` app

`ASEReadCounter`'s own output row count is **not** a reliable proxy for the
number of candidate sites offered to it — GATK's `LocusWalker` traversal
silently omits a locus with zero surviving reads (raw zero coverage, or all
reads removed by engine-level filters) rather than emitting a
`totalCount=0` row. Only sites where reads survive engine filtering but fail
ASEReadCounter's own `--min-mapping-quality`/`--min-base-quality` thresholds
get an explicit `totalCount=0` row. The `compare` app's coverage-breadth
calculation must use this app's site list / `prepare_qc_json` count as the
true denominator, not the `ASEReadCounter` TSV row count.

## Built (unpublished) on DNAnexus

| | |
|---|---|
| App ID | `app-J9PBZZ04g9qZ3Y4zvXY2Kxp2` |
| Name / version | `eggd_swap-prep-dna-vcf` / `1.0.0` |
| Jira story | [DI-3661](https://cuhbioinformatics.atlassian.net/browse/DI-3661) |
| Billed to | `org-emee_1` |
| Region | `aws:eu-central-1` |
| Built from project | `project-J9PBZG04b0g16X701bjjF6VF` (`004_260720_swap_check_apps`) |
| Developers | `org-emee_1` only |
| Authorized users | `org-emee_1` only |
| Published | No — built via `dx build --app --bill-to org-emee_1 .`, not `dx publish` |

Built following the
[DRAFT Bioinformatics Development Manual v4](https://cuhbioinformatics.atlassian.net/wiki/spaces/DV/pages/4734550017/DRAFT+Bioinformatics+Development+Manual+v4)
app requirements: app (not applet), `eggd_` prefix, `org-emee_1`-only
developers/users, Ubuntu 24.04, `aws:eu-central-1`, timeout policy set,
asset (`htslib_suite_asset`) over apt/manual compilation, `set -eo
pipefail`.

An earlier `0.1.0` build (`app-J9PBQK846y016X701bjjF5Bb`, built in
`003_260718_rna-swap-qc`) predates the DI-3661 story and the version bump
to `1.0.0`; it is superseded by this build and left as-is (unpublished,
harmless) rather than deleted.
