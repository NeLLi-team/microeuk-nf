#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob

[[ $# -eq 3 && $1 == 1 ]] || exit 1
native_output=$2
input_bins=$3
[[ -f "$native_output/checkm2.log" ]] || exit 1
[[ ! -e "$native_output/quality_report.tsv" && ! -L "$native_output/quality_report.tsv" ]] || exit 1

error_count=0
error_pattern='(^|[[:space:]])ERROR([[:space:]:]|$)'
while IFS= read -r line || [[ -n "$line" ]]; do
    if [[ $line =~ $error_pattern ]]; then
        ((error_count += 1))
        [[ "$line" == *'] ERROR: No DIAMOND annotation was generated. Exiting' ]] || exit 1
    fi
done < "$native_output/checkm2.log"
[[ $error_count -eq 1 ]] || exit 1

diamond_files=("$native_output"/diamond_output/DIAMOND_RESULTS*.tsv)
[[ ${#diamond_files[@]} -gt 0 ]] || exit 1
for path in "${diamond_files[@]}"; do
    [[ -f "$path" && ! -L "$path" && ! -s "$path" ]] || exit 1
done

bin_files=("$input_bins"/*.fa)
protein_files=("$native_output"/protein_files/*.faa)
[[ ${#bin_files[@]} -gt 0 && ${#bin_files[@]} -eq ${#protein_files[@]} ]] || exit 1
for path in "${bin_files[@]}"; do
    [[ -f "$path" && -s "$path" ]] || exit 1
    bin_name=${path##*/}
    protein="$native_output/protein_files/${bin_name%.fa}.faa"
    [[ -f "$protein" && -s "$protein" ]] || exit 1
done
