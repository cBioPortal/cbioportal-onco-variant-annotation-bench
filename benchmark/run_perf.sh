#!/bin/bash
# Timed re-runs of the local batch annotators, for the benchmark's perf table.
#
# Each tool is measured at three input sizes per repetition:
#
#   launch  10 variants in one or two places    - process start, nothing else
#   index   one variant per 1 Mb window the     - everything the tool must load
#           full input touches (381 here)         to cover the panel's footprint
#   full    the whole track input                - index cost + per-variant work
#
# A ten-variant input is NOT enough to separate these. A tool that loads its
# transcript index lazily by region only pays for the regions it is asked about,
# so ten variants understate its index cost by an order of magnitude and that
# cost then silently lands in its per-variant figure. The `index` point carries
# the same genomic footprint as the full input, so (full - index) is per-variant
# work and nothing else.
#
# Every tool is pinned to the same two cores at the same degree of parallelism,
# and every run is a cold process against a warm OS page cache: each tool gets a
# discarded warm-up run, then its repetitions back to back, so one tool's working
# set does not evict the next tool's.
#
# The two vibe-vep builds get their own HOME, because they share one
# ~/.vibe-vep transcript cache and each invalidates the other's, which otherwise
# charges whichever ran second a 12 s rebuild.
#
# vep-rs and fastVEP are each run twice -- every transcript, and --pick -- so
# the two engines are comparable in both modes rather than at different ones.
#
# $1 track key  $2 assembly  $3 track tools dir  $4 vep-rs cache name
# $5 fastVEP gff3 basename  $6 output TSV
set -u
D=${BENCH_DATA_DIR:?set BENCH_DATA_DIR to the directory holding refs/, bin/, vep_rs_cache/ and the per-track inputs}
TRACK=$1; ASM=$2; TD=$3; VCACHE=$4; GFF=$5; OUT=$6
REPS=${REPS:-3}
SCOPE=(systemd-run --user --scope --quiet -p MemoryMax=4G --)
PIN=(taskset -c 0,1)

TINY_VCF=$D/perf/tiny_$TRACK.vcf
TINY_MAF=$D/perf/tiny_$TRACK.maf
FP_VCF=$D/perf/footprint_$TRACK.vcf
FP_MAF=$D/perf/footprint_$TRACK.maf
mkdir -p "$D/perf/out"
{ grep '^#' "$TD/vep_rs/input.vcf"; grep -v '^#' "$TD/vep_rs/input.vcf" | head -10; } > "$TINY_VCF"
{ head -1 "$TD/vibe_vep/input.maf"; sed -n '2,11p' "$TD/vibe_vep/input.maf"; } > "$TINY_MAF"
python3 -I "$(dirname "$(readlink -f "$0")")/make_footprint.py" "$TD/vep_rs/input.vcf" "$FP_VCF" "$TD/vibe_vep/input.maf" "$FP_MAF"
N=$(grep -vc '^#' "$TD/vep_rs/input.vcf")
[ -s "$OUT" ] || printf 'track\ttool\tphase\trep\twall_s\tuser_s\tsys_s\tmaxrss_kb\tn_variants\n' > "$OUT"

run() { # $1 tool  $2 phase  $3 n_variants  $4 rep  ...command
  local tool=$1 phase=$2 nv=$3 rep=$4; shift 4
  local line
  line=$("${SCOPE[@]}" /usr/bin/time -f "PERF %e %U %S %M" "${PIN[@]}" "$@" 2>&1 >/dev/null \
         | grep '^PERF ' | tail -1)
  if [ -z "$line" ]; then echo "  !! FAILED $tool/$phase rep$rep" >&2; return 1; fi
  # shellcheck disable=SC2086
  set -- $line
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$TRACK" "$tool" "$phase" "$rep" "$2" "$3" "$4" "$5" "$nv" >> "$OUT"
  echo "  $tool/$phase($nv) rep$rep  wall ${2}s  cpu $(awk -v u="$3" -v s="$4" 'BEGIN{printf "%.2f", u+s}')s  rss $(( $5 / 1024 ))MB"
}

veprs=("$D/../vep-rs/target/release/vep" --json_cache "$D/vep_rs_cache/$VCACHE"
       --assembly "$ASM" --offline --fasta "$D/refs/$ASM.fa"
       --vcf --symbol --hgvs --numbers --canonical --flag_pick_allele
       --fork 2 --force_overwrite -o "$D/perf/out/veprs.vcf" -i)
veprs_pick=("${veprs[@]:0:${#veprs[@]}-1}" --pick -i)
fastvep_all=("$D/src/fastVEP/target/release/fastvep" annotate -o "$D/perf/out/fastvep.json"
         --gff3 "$D/refs/$GFF" --fasta "$D/refs/$ASM.fa" --output-format json
         --symbol --hgvs --canonical --no-progress -i)
fastvep=("${fastvep_all[@]:0:${#fastvep_all[@]}-1}" --pick -i)
vibe() { # $1 build  $2 extra-annotations(0|1)  -> echoes the command
  local cfg=(--config "$D/perf/noann.yaml"); [ "$2" = 1 ] && cfg=()
  printf '%s\0' env "HOME=$D/perf/home-$1" GOMAXPROCS=2 "$D/bin/vibe-vep-$1" \
    "${cfg[@]}" annotate maf --replace --assembly "$ASM" -o "$D/perf/out/vibe.maf"
}

# $1 tool  $2 tiny  $3 footprint  $4 full  ...command (input path appended)
block() {
  local tool=$1 tiny=$2 fp=$3 full=$4; shift 4
  local nt nf
  nt=$(count "$tiny"); nf=$(count "$fp")
  echo "=== $TRACK / $tool ==="
  "${SCOPE[@]}" "${PIN[@]}" "$@" "$full" >/dev/null 2>&1   # warm-up, discarded
  for rep in $(seq 1 "$REPS"); do
    run "$tool" launch "$nt" "$rep" "$@" "$tiny"
    run "$tool" index  "$nf" "$rep" "$@" "$fp"
    run "$tool" full   "$N"  "$rep" "$@" "$full"
  done
}

count() { case "$1" in *.vcf) grep -vc '^#' "$1";; *) echo $(( $(wc -l < "$1") - 1 ));; esac; }

V=$TD/vep_rs/input.vcf
M=$TD/vibe_vep/input.maf
block vep_rs        "$TINY_VCF" "$FP_VCF" "$V" "${veprs[@]}"
block vep_rs_pick   "$TINY_VCF" "$FP_VCF" "$V" "${veprs_pick[@]}"
block fastvep       "$TINY_VCF" "$FP_VCF" "$V" "${fastvep[@]}"
block fastvep_all   "$TINY_VCF" "$FP_VCF" "$V" "${fastvep_all[@]}"
mapfile -d '' -t vmain < <(vibe main 0)
mapfile -d '' -t vpr17 < <(vibe pr17 0)
mapfile -d '' -t vann  < <(vibe main 1)
block vibe_vep           "$TINY_MAF" "$FP_MAF" "$M" "${vmain[@]}"
block vibe_vep_pr17      "$TINY_MAF" "$FP_MAF" "$M" "${vpr17[@]}"
block vibe_vep_annotated "$TINY_MAF" "$FP_MAF" "$M" "${vann[@]}"
echo "PERFDONE $TRACK"
