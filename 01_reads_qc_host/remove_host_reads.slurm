#!/usr/bin/env bash
#SBATCH --job-name=host_remove
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=7-00:00:00
#SBATCH --output=host_removal_logs/%x_%A_%a.out
#SBATCH --error=host_removal_logs/%x_%A_%a.err

# 每个 Slurm array task 处理一个物种：
# 1) 每次只复制该物种的一个样本对到 scratch；
# 2) BWA-MEM 比对，SAMtools 提取双端均未比对的 reads；
# 3) 删除该样本在 scratch 中的输入副本和 BAM；
# 4) 整个物种完成并验证后，删除原 clean reads；
# 5) 将结果从 scratch 放回当前目录的 02_host_removed_reads；
# 6) 最终验证成功后删除该物种 scratch 目录。

if ! command -v module >/dev/null 2>&1; then
    [[ -f /etc/profile.d/modules.sh ]] && source /etc/profile.d/modules.sh
    [[ -f /usr/share/Modules/init/bash ]] && source /usr/share/Modules/init/bash
fi
command -v module >/dev/null 2>&1 || { echo "ERROR: module unavailable" >&2; exit 1; }
module load BWA
module load SAMtools

set -Eeuo pipefail
shopt -s nullglob

SCRIPT_PATH="$(readlink -f "${BASH_SOURCE[0]}")"
WORK_ROOT="${WORK_ROOT:-$PWD}"

if [[ -z "${READS_ROOT:-}" ]]; then
    if [[ -d "$WORK_ROOT/01_clean_reads" ]]; then
        READS_ROOT="$WORK_ROOT/01_clean_reads"
    elif [[ -d "$WORK_ROOT/01_clean_data/01_clean_reads" ]]; then
        READS_ROOT="$WORK_ROOT/01_clean_data/01_clean_reads"
    else
        READS_ROOT="$WORK_ROOT/01_clean_reads"
    fi
fi

GENOME_DIR="${GENOME_DIR:-/path/to/metadata/genome_masked/compress}"
SCRATCH_ROOT="${SCRATCH_ROOT:-/path/to/scratch/host_removal_work}"
REF_CACHE="${REF_CACHE:-/path/to/scratch/host_removal_reference_cache}"
RUN_NAME="${RUN_NAME:-$(basename "$WORK_ROOT")}"
RUN_SCRATCH="${RUN_SCRATCH:-$SCRATCH_ROOT/$RUN_NAME}"

# 最终结果始终回到当前运行目录。
OUT_ROOT="${OUT_ROOT:-$WORK_ROOT/02_host_removed_reads}"
REPORT_DIR="${REPORT_DIR:-$OUT_ROOT/reports}"
LOG_DIR="${LOG_DIR:-$WORK_ROOT/host_removal_logs}"

THREADS="${THREADS:-24}"
SAMTOOLS_THREADS="${SAMTOOLS_THREADS:-8}"
INDEX_MAX_CONCURRENT="${INDEX_MAX_CONCURRENT:-2}"
MAX_SPECIES_CONCURRENT="${MAX_SPECIES_CONCURRENT:-1}"
KEEP_BAM="${KEEP_BAM:-0}"
OVERWRITE="${OVERWRITE:-0}"
DELETE_CLEAN_READS="${DELETE_CLEAN_READS:-1}"
VERIFY_GZIP="${VERIFY_GZIP:-1}"

GENOME_INVENTORY="$REPORT_DIR/genome_inventory.tsv"
SPECIES_MATCH="$REPORT_DIR/species_genome_matching.tsv"
FUZZY_MATCH="$REPORT_DIR/fuzzy_genome_matches.tsv"
UNMATCHED="$REPORT_DIR/unmatched_species.tsv"
PAIRING_ISSUES="$REPORT_DIR/pairing_issues.tsv"
REFERENCE_MANIFEST="$REPORT_DIR/reference_manifest.tsv"
SPECIES_MANIFEST="$REPORT_DIR/species_manifest.tsv"
SAMPLE_MANIFEST="$REPORT_DIR/sample_manifest.tsv"
SAMPLE_STATUS_DIR="$REPORT_DIR/sample_status"
SPECIES_STATUS_DIR="$REPORT_DIR/species_status"
SAMPLE_SUMMARY="$REPORT_DIR/host_removal_sample_summary.tsv"
SPECIES_SUMMARY="$REPORT_DIR/host_removal_species_summary.tsv"
MATCH_SUMMARY="$REPORT_DIR/matching_summary.txt"

mkdir -p "$REPORT_DIR" "$SAMPLE_STATUS_DIR" "$SPECIES_STATUS_DIR"

log(){ printf '[%s] %s\n' "$(date '+%F %T')" "$*" >&2; }
die(){ log "ERROR: $*"; exit 1; }
need(){ command -v "$1" >/dev/null 2>&1 || die "Command not found: $1"; }

normalise_binomial(){
    local genus species rest
    IFS='_' read -r genus species rest <<< "$1"
    [[ -n "${genus:-}" && -n "${species:-}" ]] || return 1
    printf '%s_%s\n' "${genus^}" "${species,,}"
}

genome_taxon(){
    local b p
    b="$(basename "$1")"; p="${b%%_GCF_*}"
    [[ "$p" != "$b" ]] || p="${b%%_GCA_*}"
    [[ "$p" != "$b" ]] || return 1
    printf '%s\n' "$p"
}

manual_target(){
    case "${1,,}" in
        cephalophacus_bancanus) printf 'Cephalopachus_bancanus\n' ;;
        xanthonycticebus_pygmaeus|nycticebus_pygmaeus) printf 'Nycticebus_bengalensis\n' ;;
        rhinopithecus_roxellanae) printf 'Rhinopithecus_roxellana\n' ;;
        macaca_fuscata) printf 'Macaca_cyclopis\n' ;;
        trachypithecus_delacouri) printf 'Trachypithecus_poliocephalus\n' ;;
        cebus_capucinus) printf 'Cebus_imitator\n' ;;
        *) return 1 ;;
    esac
}

manifest_row(){ awk -F '\t' -v n="$((${2}+1))" 'NR==n{print;exit}' "$1"; }

sample_complete(){
    local s="$1" d="$2"
    [[ -s "$d/${s}_1.nonhost.fastq.gz" &&
       -s "$d/${s}_2.nonhost.fastq.gz" &&
       -s "$d/${s}.host_mapping.flagstat.txt" &&
       -s "$d/${s}.host_removed.done" ]]
}

validate_sample(){
    local s="$1" d="$2"
    sample_complete "$s" "$d" || return 1
    if [[ "$VERIFY_GZIP" == 1 ]]; then
        gzip -t "$d/${s}_1.nonhost.fastq.gz"
        gzip -t "$d/${s}_2.nonhost.fastq.gz"
    fi
}

copy_atomic(){
    local src="$1" dst="$2" gz="${3:-0}" tmp
    [[ -s "$src" ]] || die "Missing source: $src"
    mkdir -p "$(dirname "$dst")"
    tmp="${dst}.copying.${SLURM_JOB_ID:-$$}"
    rm -f "$tmp"
    rsync -a --partial "$src" "$tmp"
    [[ -s "$tmp" ]] || die "Copy failed: $src"
    [[ "$(stat -c '%s' "$src")" -eq "$(stat -c '%s' "$tmp")" ]] || die "Size mismatch: $src"
    [[ "$gz" == 1 ]] && gzip -t "$tmp"
    mv -f "$tmp" "$dst"
}

build_manifests(){
    need find; need awk; need sort
    [[ -d "$READS_ROOT" ]] || die "Reads directory not found: $READS_ROOT"
    [[ -d "$GENOME_DIR" ]] || die "Genome directory not found: $GENOME_DIR"
    mkdir -p "$REPORT_DIR" "$SAMPLE_STATUS_DIR" "$SPECIES_STATUS_DIR"

    printf 'genome_gz\tgenome_taxon\tcanonical_binomial\tgenus\tpriority\n' > "$GENOME_INVENTORY"
    local g tax canon genus pri base
    while IFS= read -r g; do
        tax="$(genome_taxon "$g")" || continue
        canon="$(normalise_binomial "$tax")" || continue
        genus="${canon%%_*}"; base="$(basename "$g")"
        [[ "$base" == *_GCF_* ]] && pri=0 || pri=1
        printf '%s\t%s\t%s\t%s\t%s\n' "$g" "$tax" "$canon" "$genus" "$pri" >> "$GENOME_INVENTORY"
    done < <(find "$GENOME_DIR" -maxdepth 1 -type f -name '*.fa.masked.gz' -print | LC_ALL=C sort)
    [[ $(wc -l < "$GENOME_INVENTORY") -gt 1 ]] || die "No genomes found"

    printf 'species_directory\tcanonical_species\tmatch_type\trequested_target\tmatched_genome_taxon\tgenome_gz\tcandidate_count\tnote\n' > "$SPECIES_MATCH"
    local sp target manual type matched genome n note exact genline
    while IFS= read -r sp; do
        canon="$(normalise_binomial "$sp")" || {
            printf '%s\tNA\tNOT_FOUND\tNA\tNA\tNA\t0\tCannot parse directory\n' "$sp" >> "$SPECIES_MATCH"; continue; }
        target="$canon"; type=""; note=""; manual="$(manual_target "$canon" 2>/dev/null || true)"
        if [[ -n "$manual" ]]; then target="$manual"; type="MANUAL_ALIAS"; note="Manual: $canon -> $target"; fi

        n="$(awk -F '\t' -v t="${target,,}" 'NR>1&&tolower($3)==t{n++}END{print n+0}' "$GENOME_INVENTORY")"
        exact="$(awk -F '\t' -v t="${target,,}" 'NR>1&&tolower($3)==t{print $5"\t"$1"\t"$2}' "$GENOME_INVENTORY" | LC_ALL=C sort -t $'\t' -k1,1n -k2,2 | head -n1)"
        if [[ -n "$exact" ]]; then
            IFS=$'\t' read -r pri genome matched <<< "$exact"
            if [[ -z "$type" ]]; then type="EXACT_BINOMIAL"; note="First two Latin-name fields matched"; fi
        else
            genus="${target%%_*}"
            n="$(awk -F '\t' -v x="${genus,,}" 'NR>1&&tolower($4)==x{n++}END{print n+0}' "$GENOME_INVENTORY")"
            genline="$(awk -F '\t' -v x="${genus,,}" 'NR>1&&tolower($4)==x{print $1"\t"$2}' "$GENOME_INVENTORY" | LC_ALL=C sort | head -n1)"
            if [[ -n "$genline" ]]; then
                IFS=$'\t' read -r genome matched <<< "$genline"
                if [[ -n "$manual" ]]; then type="MANUAL_ALIAS_GENUS_FALLBACK"; else type="GENUS_FALLBACK"; fi
                note="No exact genome; alphabetically first genome in genus $genus used"
            else
                type="NOT_FOUND"; matched="NA"; genome="NA"; n=0; note="No exact or genus genome"
            fi
        fi
        printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$sp" "$canon" "$type" "$target" "$matched" "$genome" "$n" "$note" >> "$SPECIES_MATCH"
    done < <(find "$READS_ROOT" -mindepth 1 -maxdepth 1 -type d -printf '%f\n' | LC_ALL=C sort)

    awk -F '\t' 'BEGIN{OFS="\t"}NR==1||($3!="EXACT_BINOMIAL"&&$3!="NOT_FOUND")' "$SPECIES_MATCH" > "$FUZZY_MATCH"
    awk -F '\t' 'BEGIN{OFS="\t"}NR==1||$3=="NOT_FOUND"' "$SPECIES_MATCH" > "$UNMATCHED"

    printf 'species_directory\tproject\tsample_id\tr1\texpected_r2\tissue\n' > "$PAIRING_ISSUES"
    printf 'sample_id\tspecies_directory\tcanonical_species\tproject\tr1\tr2\tmatch_type\tmatched_genome_taxon\tgenome_gz\treference_fasta\tfinal_output_directory\n' > "$SAMPLE_MANIFEST"
    local project sample r1 r2 rel out ref
    while IFS=$'\t' read -r sp canon type target matched genome n note; do
        [[ "$sp" != species_directory && "$type" != NOT_FOUND ]] || continue
        ref="$REF_CACHE/$(basename "${genome%.gz}")"
        while IFS= read -r -d '' r1; do
            r2="${r1%_1.clean.fastq.gz}_2.clean.fastq.gz"
            rel="${r1#${READS_ROOT}/${sp}/}"; project="$(dirname "$rel")"; [[ "$project" != . ]] || project=no_project
            sample="$(basename "${r1%_1.clean.fastq.gz}")"
            if [[ ! -s "$r2" ]]; then
                printf '%s\t%s\t%s\t%s\t%s\tmissing_R2\n' "$sp" "$project" "$sample" "$r1" "$r2" >> "$PAIRING_ISSUES"; continue
            fi
            out="$OUT_ROOT/reads/$sp/$project"
            printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$sample" "$sp" "$canon" "$project" "$r1" "$r2" "$type" "$matched" "$genome" "$ref" "$out" >> "$SAMPLE_MANIFEST"
        done < <(find "$READS_ROOT/$sp" -type f -name '*_1.clean.fastq.gz' -print0 | LC_ALL=C sort -z)
    done < "$SPECIES_MATCH"

    printf 'reference_id\tgenome_gz\treference_fasta\tmatched_genome_taxon\n' > "$REFERENCE_MANIFEST"
    local rid=0
    while IFS=$'\t' read -r genome matched; do
        ((rid+=1)); ref="$REF_CACHE/$(basename "${genome%.gz}")"
        printf '%s\t%s\t%s\t%s\n' "$rid" "$genome" "$ref" "$matched" >> "$REFERENCE_MANIFEST"
    done < <(awk -F '\t' 'NR>1&&$3!="NOT_FOUND"{print $6"\t"$5}' "$SPECIES_MATCH" | LC_ALL=C sort -u)

    printf 'species_task_id\tspecies_directory\tcanonical_species\tmatch_type\tmatched_genome_taxon\tgenome_gz\treference_fasta\texpected_samples\tpairing_issue_count\tsource_directory\tspecies_scratch_directory\n' > "$SPECIES_MANIFEST"
    local sid=0 expected issues
    while IFS=$'\t' read -r sp canon type target matched genome n note; do
        [[ "$sp" != species_directory && "$type" != NOT_FOUND ]] || continue
        expected="$(awk -F '\t' -v x="$sp" 'NR>1&&$2==x{n++}END{print n+0}' "$SAMPLE_MANIFEST")"
        issues="$(awk -F '\t' -v x="$sp" 'NR>1&&$1==x{n++}END{print n+0}' "$PAIRING_ISSUES")"
        ((sid+=1)); ref="$REF_CACHE/$(basename "${genome%.gz}")"
        printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$sid" "$sp" "$canon" "$type" "$matched" "$genome" "$ref" "$expected" "$issues" "$READS_ROOT/$sp" "$RUN_SCRATCH/$sp" >> "$SPECIES_MANIFEST"
    done < "$SPECIES_MATCH"

    {
        printf 'READS_ROOT=%s\nGENOME_DIR=%s\nRUN_SCRATCH=%s\nREF_CACHE=%s\nOUT_ROOT=%s\n\n' "$READS_ROOT" "$GENOME_DIR" "$RUN_SCRATCH" "$REF_CACHE" "$OUT_ROOT"
        awk -F '\t' 'NR>1{n[$3]++}END{for(k in n)print k"\t"n[k]}' "$SPECIES_MATCH" | sort
        printf '\nreferences=%s\nspecies_tasks=%s\npaired_samples=%s\npairing_issues=%s\nMAX_SPECIES_CONCURRENT=%s\n' "$rid" "$sid" "$(( $(wc -l < "$SAMPLE_MANIFEST")-1 ))" "$(( $(wc -l < "$PAIRING_ISSUES")-1 ))" "$MAX_SPECIES_CONCURRENT"
    } > "$MATCH_SUMMARY"
    cat "$MATCH_SUMMARY" >&2
}

index_task(){
    need bwa; need samtools; need gzip
    local task="${SLURM_ARRAY_TASK_ID:-${2:-}}" line rid gz ref matched tmp
    [[ -n "$task" ]] || die "Missing reference task ID"
    line="$(manifest_row "$REFERENCE_MANIFEST" "$task")"; [[ -n "$line" ]] || die "Invalid reference task $task"
    IFS=$'\t' read -r rid gz ref matched <<< "$line"
    mkdir -p "$REF_CACHE"; gzip -t "$gz"
    if [[ ! -s "$ref" ]]; then
        tmp="${ref}.tmp.${SLURM_JOB_ID:-$$}"; rm -f "$tmp"; gzip -dc "$gz" > "$tmp"; [[ -s "$tmp" ]] || die "Reference decompression failed"; mv -f "$tmp" "$ref"
    fi
    if [[ ! -s "$ref.amb" || ! -s "$ref.ann" || ! -s "$ref.bwt" || ! -s "$ref.pac" || ! -s "$ref.sa" ]]; then
        rm -f "$ref.amb" "$ref.ann" "$ref.bwt" "$ref.pac" "$ref.sa"; bwa index "$ref"
    fi
    [[ -s "$ref.fai" ]] || samtools faidx "$ref"
    log "Reference ready: $matched"
}

species_task(){
    need bwa; need samtools; need rsync; need gzip; need stat
    local task="${SLURM_ARRAY_TASK_ID:-${2:-}}" line sid sp canon type matched genome ref expected issues source scratch
    [[ -n "$task" ]] || die "Missing species task ID"
    line="$(manifest_row "$SPECIES_MANIFEST" "$task")"; [[ -n "$line" ]] || die "Invalid species task $task"
    IFS=$'\t' read -r sid sp canon type matched genome ref expected issues source scratch <<< "$line"
    [[ -s "$ref" && -s "$ref.bwt" && -s "$ref.sa" ]] || die "Reference/index missing: $ref"

    local sin="$scratch/inputs" sout="$scratch/results" stmp="$scratch/tmp"
    mkdir -p "$sin" "$sout" "$stmp" "$SAMPLE_STATUS_DIR" "$SPECIES_STATUS_DIR"
    local sp_status="$SPECIES_STATUS_DIR/${sid}_${sp}.tsv" start="$(date +%s)" files_deleted=0 bytes_deleted=0
    status(){ printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$sid" "$sp" "$canon" "$type" "$matched" "$expected" "$files_deleted" "$bytes_deleted" "$1" "$(( $(date +%s)-start ))" "$scratch" "$2" > "$sp_status"; }

    ((issues==0)) || { status RETAINED_PAIRING_ISSUE "$issues pairing issue(s)"; return 1; }
    ((expected>0)) || { status RETAINED_NO_SAMPLES "No paired samples"; return 1; }

    log "Species $sp; scratch=$scratch; final=$OUT_ROOT/reads/$sp"
    local sample project r1 r2 final sproj sinput sr1 sr2 bam o1 o2 tmp1 tmp2 flag done statfile
    local primary pairs unmapped nonhost removed pct t0 elapsed

    while IFS=$'\t' read -r sample _sp _canon project r1 r2 _type _matched _genome _ref final; do
        [[ "$_sp" == "$sp" ]] || continue
        sproj="$sout/$project"; sinput="$sin/$project/$sample"; mkdir -p "$sproj" "$sinput" "$final"
        statfile="$SAMPLE_STATUS_DIR/${sp}_${project//\//_}_${sample}.tsv"

        if [[ "$OVERWRITE" != 1 ]] && validate_sample "$sample" "$final"; then
            printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\tNA\tNA\tNA\tNA\tREUSED_FINAL\t0\n' "$sample" "$sp" "$canon" "$project" "$type" "$matched" "$ref" > "$statfile"; continue
        fi
        if [[ "$OVERWRITE" != 1 ]] && validate_sample "$sample" "$sproj"; then
            printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\tNA\tNA\tNA\tNA\tREUSED_SCRATCH\t0\n' "$sample" "$sp" "$canon" "$project" "$type" "$matched" "$ref" > "$statfile"; continue
        fi
        [[ -s "$r1" && -s "$r2" ]] || { status RETAINED_MISSING_SOURCE "Missing $project/$sample"; die "Missing source reads: $sp/$project/$sample"; }

        sr1="$sinput/$(basename "$r1")"; sr2="$sinput/$(basename "$r2")"; t0="$(date +%s)"
        copy_atomic "$r1" "$sr1" 1; copy_atomic "$r2" "$sr2" 1

        bam="$stmp/${project//\//_}_${sample}.all.bam"
        tmp1="$stmp/${project//\//_}_${sample}_1.nonhost.fastq.gz"
        tmp2="$stmp/${project//\//_}_${sample}_2.nonhost.fastq.gz"
        o1="$sproj/${sample}_1.nonhost.fastq.gz"; o2="$sproj/${sample}_2.nonhost.fastq.gz"
        flag="$sproj/${sample}.host_mapping.flagstat.txt"; done="$sproj/${sample}.host_removed.done"
        rm -f "$bam" "$tmp1" "$tmp2" "$o1" "$o2" "$flag" "$done"

        bwa mem -t "$THREADS" "$ref" "$sr1" "$sr2" | samtools view -@ "$SAMTOOLS_THREADS" -b -o "$bam" -
        samtools quickcheck -v "$bam"; samtools flagstat -@ "$SAMTOOLS_THREADS" "$bam" > "$flag"
        primary="$(samtools view -@ "$SAMTOOLS_THREADS" -c -f 1 -F 0x900 "$bam")"; pairs=$((primary/2))
        unmapped="$(samtools view -@ "$SAMTOOLS_THREADS" -c -f 12 -F 0x900 "$bam")"; nonhost=$((unmapped/2)); removed=$((pairs-nonhost))
        pct="$(awk -v n="$nonhost" -v d="$pairs" 'BEGIN{if(d>0)printf "%.4f",100*n/d;else print "NA"}')"

        samtools collate -@ "$SAMTOOLS_THREADS" -u -O -T "$stmp/${project//\//_}_${sample}.collate" "$bam" |
            samtools fastq -@ "$SAMTOOLS_THREADS" -f 12 -F 0x900 -n -1 "$tmp1" -2 "$tmp2" -0 /dev/null -s /dev/null -
        [[ -s "$tmp1" && -s "$tmp2" ]] || die "Empty non-host output: $sp/$sample"
        gzip -t "$tmp1"; gzip -t "$tmp2"; mv -f "$tmp1" "$o1"; mv -f "$tmp2" "$o2"

        if [[ "$KEEP_BAM" == 1 ]]; then mkdir -p "$sout/bam/$project"; mv -f "$bam" "$sout/bam/$project/${sample}.all.bam"; else rm -f "$bam"; fi
        printf 'sample=%s\nnonhost_pairs=%s\ncompleted=%s\n' "$sample" "$nonhost" "$(date -Is)" > "$done"
        validate_sample "$sample" "$sproj" || die "Scratch validation failed: $sp/$sample"
        rm -f "$sr1" "$sr2" "$stmp/${project//\//_}_${sample}.collate"* 2>/dev/null || true
        find "$sinput" -depth -type d -empty -delete 2>/dev/null || true
        elapsed=$(( $(date +%s)-t0 ))
        printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\tSUCCESS_SCRATCH\t%s\n' "$sample" "$sp" "$canon" "$project" "$type" "$matched" "$ref" "$pairs" "$nonhost" "$removed" "$pct" "$elapsed" > "$statfile"
    done < <(tail -n +2 "$SAMPLE_MANIFEST")

    local completed=0 scratch_ready=0 final_ready=0
    while IFS=$'\t' read -r sample _sp _canon project r1 r2 _type _matched _genome _ref final; do
        [[ "$_sp" == "$sp" ]] || continue; sproj="$sout/$project"
        if validate_sample "$sample" "$final"; then ((final_ready+=1)); elif validate_sample "$sample" "$sproj"; then ((scratch_ready+=1)); else status RETAINED_INCOMPLETE_OUTPUT "Incomplete $project/$sample"; return 1; fi
        ((completed+=1))
    done < <(tail -n +2 "$SAMPLE_MANIFEST")
    [[ "$completed" -eq "$expected" ]] || { status RETAINED_COMPLETION_MISMATCH "$completed/$expected complete"; return 1; }

    if [[ "$DELETE_CLEAN_READS" == 1 ]]; then
        local -a del=(); local f size
        if [[ -d "$source" ]]; then
            while IFS= read -r -d '' f; do del+=("$f"); size="$(stat -c '%s' "$f")"; bytes_deleted=$((bytes_deleted+size)); done < <(find "$source" -type f \( -name '*_1.clean.fastq.gz' -o -name '*_2.clean.fastq.gz' \) -print0)
        fi
        files_deleted="${#del[@]}"
        if [[ "$files_deleted" -ne 0 && "$files_deleted" -ne $((expected*2)) ]]; then files_deleted=0; bytes_deleted=0; status RETAINED_DELETE_LIST_MISMATCH "Found ${#del[@]}, expected $((expected*2))"; return 1; fi
        if ((files_deleted>0)); then rm -f -- "${del[@]}"; find "$source" -depth -type d -empty -delete 2>/dev/null || true; fi
    fi

    while IFS=$'\t' read -r sample _sp _canon project r1 r2 _type _matched _genome _ref final; do
        [[ "$_sp" == "$sp" ]] || continue
        validate_sample "$sample" "$final" && continue
        sproj="$sout/$project"; validate_sample "$sample" "$sproj" || { status TRANSFER_FAILED_SCRATCH_RETAINED "Missing scratch $project/$sample"; return 1; }
        copy_atomic "$sproj/${sample}_1.nonhost.fastq.gz" "$final/${sample}_1.nonhost.fastq.gz" 1
        copy_atomic "$sproj/${sample}_2.nonhost.fastq.gz" "$final/${sample}_2.nonhost.fastq.gz" 1
        copy_atomic "$sproj/${sample}.host_mapping.flagstat.txt" "$final/${sample}.host_mapping.flagstat.txt" 0
        copy_atomic "$sproj/${sample}.host_removed.done" "$final/${sample}.host_removed.done" 0
        validate_sample "$sample" "$final" || { status TRANSFER_FAILED_SCRATCH_RETAINED "Final validation failed $project/$sample"; return 1; }
    done < <(tail -n +2 "$SAMPLE_MANIFEST")

    if [[ "$KEEP_BAM" == 1 && -d "$sout/bam" ]]; then mkdir -p "$OUT_ROOT/bam/$sp"; rsync -a "$sout/bam/" "$OUT_ROOT/bam/$sp/"; fi

    completed=0
    while IFS=$'\t' read -r sample _sp _canon project r1 r2 _type _matched _genome _ref final; do
        [[ "$_sp" == "$sp" ]] || continue; validate_sample "$sample" "$final" || { status FINAL_VALIDATION_FAILED_SCRATCH_RETAINED "$project/$sample"; return 1; }; ((completed+=1))
    done < <(tail -n +2 "$SAMPLE_MANIFEST")
    [[ "$completed" -eq "$expected" ]] || { status FINAL_COMPLETION_MISMATCH_SCRATCH_RETAINED "$completed/$expected"; return 1; }

    rm -rf "$scratch"
    if [[ "$DELETE_CLEAN_READS" == 1 ]]; then status SUCCESS_DELETED_SOURCE "Final: $OUT_ROOT/reads/$sp"; else status SUCCESS_RETAINED_SOURCE "Final: $OUT_ROOT/reads/$sp"; fi
    log "Completed $sp -> $OUT_ROOT/reads/$sp"
}

collect_reports(){
    printf 'sample_id\tspecies_directory\tcanonical_species\tproject\tmatch_type\tmatched_genome_taxon\treference_fasta\tinput_pairs\tnonhost_pairs\tremoved_pairs\tretained_percent\tstatus\telapsed_seconds\n' > "$SAMPLE_SUMMARY"
    local a=("$SAMPLE_STATUS_DIR"/*.tsv); ((${#a[@]}==0)) || LC_ALL=C sort "${a[@]}" >> "$SAMPLE_SUMMARY"
    printf 'species_task_id\tspecies_directory\tcanonical_species\tmatch_type\tmatched_genome_taxon\texpected_samples\tfiles_deleted\tbytes_deleted\tstatus\telapsed_seconds\tscratch_directory\tnote\n' > "$SPECIES_SUMMARY"
    local b=("$SPECIES_STATUS_DIR"/*.tsv); ((${#b[@]}==0)) || LC_ALL=C sort -t $'\t' -k1,1n "${b[@]}" >> "$SPECIES_SUMMARY"
    { printf 'Final output: %s\nScratch: %s\n\n' "$OUT_ROOT" "$RUN_SCRATCH"; printf 'Match types:\n'; awk -F '\t' 'NR>1{n[$3]++}END{for(k in n)print k"\t"n[k]}' "$SPECIES_MATCH" | sort; printf '\nSpecies statuses:\n'; awk -F '\t' 'NR>1{n[$9]++}END{for(k in n)print k"\t"n[k]}' "$SPECIES_SUMMARY" | sort; } > "$MATCH_SUMMARY"
    cat "$MATCH_SUMMARY"
}

export_vars(){
    printf 'ALL,WORK_ROOT=%s,READS_ROOT=%s,GENOME_DIR=%s,SCRATCH_ROOT=%s,REF_CACHE=%s,RUN_NAME=%s,RUN_SCRATCH=%s,OUT_ROOT=%s,REPORT_DIR=%s,LOG_DIR=%s,THREADS=%s,SAMTOOLS_THREADS=%s,INDEX_MAX_CONCURRENT=%s,MAX_SPECIES_CONCURRENT=%s,KEEP_BAM=%s,OVERWRITE=%s,DELETE_CLEAN_READS=%s,VERIFY_GZIP=%s' "$WORK_ROOT" "$READS_ROOT" "$GENOME_DIR" "$SCRATCH_ROOT" "$REF_CACHE" "$RUN_NAME" "$RUN_SCRATCH" "$OUT_ROOT" "$REPORT_DIR" "$LOG_DIR" "$THREADS" "$SAMTOOLS_THREADS" "$INDEX_MAX_CONCURRENT" "$MAX_SPECIES_CONCURRENT" "$KEEP_BAM" "$OVERWRITE" "$DELETE_CLEAN_READS" "$VERIFY_GZIP"
}

submit_pipeline(){
    need sbatch; cd "$WORK_ROOT"; mkdir -p "$LOG_DIR" "$RUN_SCRATCH" "$REF_CACHE"; build_manifests
    local nr ns ex ij sj cj ev; nr=$(( $(wc -l < "$REFERENCE_MANIFEST")-1 )); ns=$(( $(wc -l < "$SPECIES_MANIFEST")-1 )); ((nr>0&&ns>0)) || die "No runnable tasks"; ev="$(export_vars)"
    ex="$(sbatch --parsable --export="$ev" --job-name=host_index --cpus-per-task=1 --mem=48G --array="1-${nr}%${INDEX_MAX_CONCURRENT}" "$SCRIPT_PATH" index)"; ij="${ex%%;*}"
    ex="$(sbatch --parsable --export="$ev" --job-name=host_species --dependency="afterok:$ij" --array="1-${ns}%${MAX_SPECIES_CONCURRENT}" "$SCRIPT_PATH" species)"; sj="${ex%%;*}"
    ex="$(sbatch --parsable --export="$ev" --job-name=host_collect --dependency="afterany:$sj" --cpus-per-task=1 --mem=2G --time=00:15:00 "$SCRIPT_PATH" collect)"; cj="${ex%%;*}"
    printf '3 Slurm job IDs submitted:\n  index array:   %s (%s tasks; max %s)\n  species array: %s (%s tasks; max %s)\n  collect:       %s\n' "$ij" "$nr" "$INDEX_MAX_CONCURRENT" "$sj" "$ns" "$MAX_SPECIES_CONCURRENT" "$cj"
}

resume_pipeline(){
    need sbatch; [[ -s "$SPECIES_MANIFEST" && -s "$SAMPLE_MANIFEST" ]] || die "Existing manifests required"
    local ns ex sj cj ev; ns=$(( $(wc -l < "$SPECIES_MANIFEST")-1 )); ev="$(export_vars)"; mkdir -p "$LOG_DIR"
    ex="$(sbatch --parsable --export="$ev" --job-name=host_resume --array="1-${ns}%${MAX_SPECIES_CONCURRENT}" "$SCRIPT_PATH" species)"; sj="${ex%%;*}"
    ex="$(sbatch --parsable --export="$ev" --job-name=host_collect --dependency="afterany:$sj" --cpus-per-task=1 --mem=2G --time=00:15:00 "$SCRIPT_PATH" collect)"; cj="${ex%%;*}"
    printf '2 resume job IDs submitted:\n  species array: %s\n  collect: %s\n' "$sj" "$cj"
}

help(){ cat <<EOF
Usage:
  bash $(basename "$SCRIPT_PATH") manifest
  bash $(basename "$SCRIPT_PATH") submit
  bash $(basename "$SCRIPT_PATH") resume
  bash $(basename "$SCRIPT_PATH") collect

Final output:
  $OUT_ROOT/reads/<species>/<project>/

Scratch:
  $RUN_SCRATCH/<species>/

Recommended serial run:
  MAX_SPECIES_CONCURRENT=1 \\
  DELETE_CLEAN_READS=1 \\
  VERIFY_GZIP=1 \\
  bash $(basename "$SCRIPT_PATH") submit

Two species at the same time:
  MAX_SPECIES_CONCURRENT=2 \\
  DELETE_CLEAN_READS=1 \\
  VERIFY_GZIP=1 \\
  bash $(basename "$SCRIPT_PATH") submit

If transfer was interrupted after source deletion, do not rebuild manifests:
  MAX_SPECIES_CONCURRENT=1 bash $(basename "$SCRIPT_PATH") resume
EOF
}

case "${1:-help}" in
    manifest|prepare) build_manifests ;;
    index) index_task "$@" ;;
    species|run) species_task "$@" ;;
    collect) collect_reports ;;
    submit) submit_pipeline ;;
    resume) resume_pipeline ;;
    help|-h|--help) help ;;
    *) help; die "Unknown mode: ${1:-}" ;;
esac
