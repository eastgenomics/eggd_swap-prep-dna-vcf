#!/bin/bash
# eggd_swap-prep-dna-vcf — select germline het SNVs from a TSO500 DNA VCF,
# lift over to GRCh38 only when the input actually needs it, restrict to
# the RNA capture region, and report a per-step site-count funnel.
#
# Deliberately simple: no asset manifest, no candidate_id, no validation_mode.
# Every input is a plain file this app is directly given; every failure is a
# real tool/data problem, not a platform-contract negotiation.
set -eo pipefail

main() {
    local work_dir
    work_dir="$(mktemp -d)"
    trap 'rm -rf "$work_dir"' EXIT

    dx-download-all-inputs --parallel

    local dna_vcf_path chain_path capture_bed_path gatk_jar_path
    dna_vcf_path="$HOME/in/dna_vcf/$(ls "$HOME/in/dna_vcf")"
    chain_path="$HOME/in/chain_file/$(ls "$HOME/in/chain_file")"
    capture_bed_path="$HOME/in/capture_bed/$(ls "$HOME/in/capture_bed")"
    gatk_jar_path="$HOME/in/gatk_jar/$(ls "$HOME/in/gatk_jar")"

    # GATK requires the .fai and .dict to sit next to the .fasta with matching
    # basenames; dx-download-all-inputs keeps each declared input in its own
    # subfolder, so build one consistent reference bundle by symlink here
    # rather than relying on any coincidental naming from the three inputs.
    mkdir -p "$work_dir/ref"
    ln -s "$HOME/in/reference_fasta/$(ls "$HOME/in/reference_fasta")" "$work_dir/ref/reference.fasta"
    ln -s "$HOME/in/reference_fasta_fai/$(ls "$HOME/in/reference_fasta_fai")" "$work_dir/ref/reference.fasta.fai"
    ln -s "$HOME/in/reference_dict/$(ls "$HOME/in/reference_dict")" "$work_dir/ref/reference.dict"
    local reference_fasta_path="$work_dir/ref/reference.fasta"
    local reference_fai_path="$work_dir/ref/reference.fasta.fai"

    # --- Step 0: determine genome build from the VCF's own header ----------
    # Never assumed from a filename. Fails closed (non-zero exit, no output
    # produced) if the header's contig-length signal is missing, unrecognised,
    # or contradicted by its own ##reference= line or by the supplied
    # reference_fasta_fai.
    local build_json="$work_dir/build.json"
    python3 /home/dnanexus/detect_build.py "$dna_vcf_path" "$reference_fai_path" >"$build_json"
    local liftover_required
    liftover_required="$(python3 -c "import json,sys; print(json.load(open(sys.argv[1]))['liftover_required'])" "$build_json")"
    echo "Detected genome build: $(python3 -c "import json,sys; print(json.load(open(sys.argv[1]))['detected_build'])" "$build_json"); liftover_required=${liftover_required}"

    local raw_count selected_count vafband_count renamed_count lifted_count rejected_count post_mnv_count capture_count final_count

    raw_count="$(bcftools view -H "$dna_vcf_path" | wc -l)"

    # --- Step 1+2: PASS, biallelic SNV, het GT, DP>=30, VAF band 0.35-0.65 --
    bcftools view -f PASS -v snps -m2 -M2 "$dna_vcf_path" \
        | bcftools view -i 'GT="het"' -Ou \
        | bcftools view -i 'FORMAT/DP>=30 && FORMAT/VF>=0.35 && FORMAT/VF<=0.65' -Oz \
            -o "$work_dir/selected.vcf.gz"
    bcftools index -t "$work_dir/selected.vcf.gz"
    selected_count="$(bcftools view -H "$work_dir/selected.vcf.gz" | wc -l)"
    vafband_count="$selected_count"

    # --- Step 3: normalise contig naming (idempotent either way) ------------
    # A chr-prefixed source is required by both the liftover chain (hg19 side)
    # and the CTAT/GRCh38 reference (already-GRCh38 side), so this step runs
    # unconditionally rather than being gated on the detected build.
    {
        for n in $(seq 1 22); do printf '%s\tchr%s\n' "$n" "$n"; done
        printf 'X\tchrX\nY\tchrY\nMT\tchrM\nM\tchrM\n'
    } >"$work_dir/rename-chrs.tsv"
    bcftools annotate --rename-chrs "$work_dir/rename-chrs.tsv" "$work_dir/selected.vcf.gz" -Oz \
        -o "$work_dir/renamed.vcf.gz"
    bcftools index -t "$work_dir/renamed.vcf.gz"
    renamed_count="$(bcftools view -H "$work_dir/renamed.vcf.gz" | wc -l)"

    local lifted_vcf="$work_dir/lifted.vcf.gz"
    local rejected_vcf="$work_dir/rejected.vcf.gz"

    if [[ "$liftover_required" == "True" ]]; then
        # --- Step 4: liftover hg19 -> GRCh38 --------------------------------
        # GATK requires --REJECT to point at a real path; its content is not
        # persisted as a declared output (see rejected_count logging below) —
        # only the count and GATK's own per-contig log lines, already visible
        # in the job log, are kept.
        java -XX:MaxRAMPercentage=70.0 -jar "$gatk_jar_path" LiftoverVcf \
            -I "$work_dir/renamed.vcf.gz" \
            -O "$lifted_vcf" \
            --CHAIN "$chain_path" \
            --REJECT "$rejected_vcf" \
            -R "$reference_fasta_path" \
            --RECOVER_SWAPPED_REF_ALT true \
            --WRITE_ORIGINAL_POSITION true
        rejected_count="$(bcftools view -H "$rejected_vcf" 2>/dev/null | wc -l || echo 0)"
    else
        # dna_vcf is already GRCh38 (confirmed in Step 0 against the supplied
        # reference_fasta_fai) — skip GATK LiftoverVcf entirely. There is
        # nothing to reject, so no rejected_vcf is produced at all.
        cp "$work_dir/renamed.vcf.gz" "$lifted_vcf"
        rejected_count=0
    fi
    bcftools index -t "$lifted_vcf" 2>/dev/null || true
    lifted_count="$(bcftools view -H "$lifted_vcf" | wc -l)"
    echo "LiftoverVcf rejected ${rejected_count} record(s)$( [[ "$liftover_required" == "True" ]] || echo ' (liftover skipped — dna_vcf already GRCh38)' )."

    # --- Step 5: exclude MNV-type records (overlapping-record crash fix) ----
    bcftools view -e '(strlen(REF)>1 && strlen(REF)=strlen(ALT))' "$lifted_vcf" -Oz \
        -o "$work_dir/no-mnv.vcf.gz"
    bcftools index -t "$work_dir/no-mnv.vcf.gz"
    post_mnv_count="$(bcftools view -H "$work_dir/no-mnv.vcf.gz" | wc -l)"

    # --- Step 6: restrict to the RNA capture region --------------------------
    bcftools view -R "$capture_bed_path" "$work_dir/no-mnv.vcf.gz" -Oz \
        -o "$work_dir/site.vcf.gz"
    bcftools index -t "$work_dir/site.vcf.gz"
    capture_count="$(bcftools view -H "$work_dir/site.vcf.gz" | wc -l)"
    final_count="$capture_count"

    tabix -p vcf "$work_dir/site.vcf.gz"

    echo "Site funnel: raw=${raw_count} selected(PASS,biallelic,het,DP>=30)=${selected_count} vaf_band=${vafband_count} renamed=${renamed_count} lifted_or_grch38_input=${lifted_count} rejected_by_liftover=${rejected_count} post_mnv_exclusion=${post_mnv_count} capture_region=${capture_count} final=${final_count}"

    # --- QC / funnel report ---------------------------------------------------
    local qc_json="$work_dir/prepare_qc.json"
    python3 /home/dnanexus/qc_json.py \
        --sample-id "${sample_id:-}" \
        --build-json "$build_json" \
        --count "raw=${raw_count}" \
        --count "selected_pass_biallelic_het_dp30=${selected_count}" \
        --count "vaf_band=${vafband_count}" \
        --count "renamed=${renamed_count}" \
        --count "lifted_or_grch38_input=${lifted_count}" \
        --count "rejected_by_liftover=${rejected_count}" \
        --count "post_mnv_exclusion=${post_mnv_count}" \
        --count "capture_region=${capture_count}" \
        --count "final=${final_count}" \
        >"$qc_json"

    local site_vcf_id site_vcf_tbi_id qc_json_id
    site_vcf_id="$(dx upload "$work_dir/site.vcf.gz" --brief)"
    site_vcf_tbi_id="$(dx upload "$work_dir/site.vcf.gz.tbi" --brief)"
    qc_json_id="$(dx upload "$qc_json" --brief)"

    dx-jobutil-add-output site_vcf "$site_vcf_id" --class=file
    dx-jobutil-add-output site_vcf_tbi "$site_vcf_tbi_id" --class=file
    dx-jobutil-add-output prepare_qc_json "$qc_json_id" --class=file
}

main
