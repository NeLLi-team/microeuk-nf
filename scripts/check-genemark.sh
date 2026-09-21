#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then
    printf 'usage: %s PROOF_DIR\n' "${0##*/}" >&2
    exit 2
fi

project_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
manifest="${project_root}/envs/gene/longrna/pixi.toml"
prefix="${project_root}/envs/gene/longrna/.pixi/envs/default"
proof_dir="$1"
mkdir -p "${proof_dir}"
proof_dir="$(cd "${proof_dir}" && pwd -P)"
if [[ -n "$(find "${proof_dir}" -mindepth 1 -maxdepth 1 -print -quit)" ]]; then
    printf 'proof directory must be empty: %s\n' "${proof_dir}" >&2
    exit 2
fi
mkdir "${proof_dir}/empty-home"

[[ -x "${prefix}/bin/braker3-longrna-preflight" ]]
[[ -x "${prefix}/share/braker3-longrna/ETP/bin/gmes/gmhmme3" ]]

PIXI_NO_INSTALL=true PIXI_FROZEN=true pixi run --as-is \
    --manifest-path "${manifest}" bash -s -- \
    "${proof_dir}" "${prefix}" <<'SMOKE'
set -euo pipefail

proof_dir="$1"
prefix="$2"
[[ "$(realpath "$(command -v braker3-longrna-preflight)")" \
    == "$(realpath "${prefix}/bin/braker3-longrna-preflight")" ]]
test_root="${prefix}/share/braker3-longrna/ETP/bin/gmes/GeneMark-E-tests/GeneMark.hmm"
gmhmme3="${prefix}/share/braker3-longrna/ETP/bin/gmes/gmhmme3"
comparator="${prefix}/share/braker3-longrna/ETP/bin/gmes/compare_intervals_exact.pl"
expected="${test_root}/output/genemark.gff3"
predicted="${proof_dir}/genemark.gff3"
repeated="${proof_dir}/genemark.repeat.gff3"

for required in \
    "${gmhmme3}" \
    "${comparator}" \
    "${test_root}/input/athaliana.mod" \
    "${test_root}/input/sequence.fasta" \
    "${expected}"; do
    [[ -s "${required}" ]]
done

HOME="${proof_dir}/empty-home" "${gmhmme3}" \
    -o "${predicted}" \
    -m "${test_root}/input/athaliana.mod" \
    -f gff3 \
    "${test_root}/input/sequence.fasta" \
    > "${proof_dir}/prediction.stdout" \
    2> "${proof_dir}/prediction.stderr"
[[ -s "${predicted}" ]]

HOME="${proof_dir}/empty-home" "${gmhmme3}" \
    -o "${repeated}" \
    -m "${test_root}/input/athaliana.mod" \
    -f gff3 \
    "${test_root}/input/sequence.fasta" \
    > "${proof_dir}/prediction.repeat.stdout" \
    2> "${proof_dir}/prediction.repeat.stderr"
[[ -s "${repeated}" ]]

sequence_length=$(awk '!/^>/ {gsub(/[[:space:]]/, ""); total += length} END {print total}' \
    "${test_root}/input/sequence.fasta")
awk -v sequence_length="${sequence_length}" '
    BEGIN {FS = "\t"}
    /^#/ || /^[[:space:]]*$/ {next}
    NF != 9 {exit 10}
    $1 != "seq" || $4 !~ /^[0-9]+$/ || $5 !~ /^[0-9]+$/ {exit 11}
    $4 < 1 || $5 < $4 || $5 > sequence_length {exit 12}
    $6 != "." || $7 !~ /^[+-]$/ {exit 13}
    tolower($3) !~ /^(gene|mrna|exon|cds|intron)$/ {exit 14}
    tolower($3) ~ /^(cds|intron)$/ && $8 !~ /^[012]$/ {exit 15}
    tolower($3) !~ /^(cds|intron)$/ && $8 != "." {exit 16}
    $9 !~ /(^|;)ID=[^;]+/ {exit 17}
    tolower($3) == "mrna" && $9 !~ /(^|;)Parent=gene[0-9]+/ {exit 18}
    tolower($3) ~ /^(exon|cds|intron)$/ && $9 !~ /(^|;)Parent=mRNA[0-9]+/ {exit 19}
    {count[tolower($3)]++}
    END {
        if (count["gene"] != 313 || count["mrna"] != 313 ||
            count["exon"] != 1527 || count["cds"] != 1527 ||
            count["intron"] != 1214) exit 20
        print "feature\tcount"
        print "gene\t" count["gene"]
        print "mRNA\t" count["mrna"]
        print "exon\t" count["exon"]
        print "CDS\t" count["cds"]
        print "intron\t" count["intron"]
    }
' "${predicted}" > "${proof_dir}/structure.tsv"

grep -v '^#' "${predicted}" > "${proof_dir}/predicted.feature-rows.tsv"
grep -v '^#' "${repeated}" > "${proof_dir}/repeated.feature-rows.tsv"
cmp "${proof_dir}/predicted.feature-rows.tsv" \
    "${proof_dir}/repeated.feature-rows.tsv"

canonical_cds() {
    awk 'BEGIN {FS=OFS="\t"} !/^#/ && $3 == "CDS" {print $1, $3, $4, $5, $7, $8}' "$1" \
        | LC_ALL=C sort
}

canonical_cds "${expected}" > "${proof_dir}/expected.cds.tsv"
canonical_cds "${predicted}" > "${proof_dir}/predicted.cds.tsv"
[[ -s "${proof_dir}/predicted.cds.tsv" ]]
comm -23 "${proof_dir}/expected.cds.tsv" "${proof_dir}/predicted.cds.tsv" \
    > "${proof_dir}/expected-only.cds.tsv"
comm -13 "${proof_dir}/expected.cds.tsv" "${proof_dir}/predicted.cds.tsv" \
    > "${proof_dir}/predicted-only.cds.tsv"

mapfile -t expected_only < "${proof_dir}/expected-only.cds.tsv"
mapfile -t predicted_only < "${proof_dir}/predicted-only.cds.tsv"
[[ "${#expected_only[@]}" -eq 2 ]]
[[ "${#predicted_only[@]}" -eq 2 ]]
[[ "${expected_only[0]}" == $'seq\tCDS\t264\t864\t-\t0' ]]
[[ "${expected_only[1]}" == $'seq\tCDS\t78\t148\t-\t2' ]]
[[ "${predicted_only[0]}" == $'seq\tCDS\t256\t864\t-\t0' ]]
[[ "${predicted_only[1]}" == $'seq\tCDS\t998925\t999045\t-\t0' ]]

diff_status=0
diff -u "${proof_dir}/expected.cds.tsv" "${proof_dir}/predicted.cds.tsv" \
    > "${proof_dir}/expected-vs-predicted.cds.diff" || diff_status=$?
[[ "${diff_status}" -eq 1 ]]

perl "${comparator}" --f1 "${predicted}" --f2 "${expected}" --v \
    > "${proof_dir}/compare.txt"

predicted_version=$(sed -n \
    's/^# Eukaryotic GeneMark\.hmm version //p' "${predicted}")
expected_version=$(sed -n \
    's/^# Eukaryotic GeneMark\.hmm version //p' "${expected}")
[[ "${predicted_version}" == "3.68" ]]
[[ "${expected_version}" == "3.62_lic" ]]

mapfile -t comparison_rows < <(
    awk -F '\t' '$1 ~ /^[0-9]+$/ {print $1, $2, $3, $4}' \
        "${proof_dir}/compare.txt"
)
[[ "${#comparison_rows[@]}" -eq 2 ]]
[[ "${comparison_rows[0]}" == "1527 1525 2 99.87" ]]
[[ "${comparison_rows[1]}" == "1527 1525 2 99.87" ]]

sha256sum \
    "${predicted}" \
    "${repeated}" \
    "${proof_dir}/predicted.feature-rows.tsv" \
    "${proof_dir}/expected.cds.tsv" \
    "${proof_dir}/predicted.cds.tsv" \
    "${proof_dir}/expected-only.cds.tsv" \
    "${proof_dir}/predicted-only.cds.tsv" \
    "${proof_dir}/expected-vs-predicted.cds.diff" \
    "${proof_dir}/compare.txt" \
    "${proof_dir}/structure.tsv" \
    > "${proof_dir}/sha256sums.txt"
printf 'GeneMark %s prediction smoke passed: 1527 CDS; 1525 match the %s upstream reference; repeated feature rows are exact.\n' \
    "${predicted_version}" "${expected_version}"
SMOKE
