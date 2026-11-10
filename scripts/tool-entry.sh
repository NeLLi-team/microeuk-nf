#!/usr/bin/env bash
# Print one `tools` map entry of a merged stage.json from a tool process report.
#
# The report is one tab-separated line: version, status, reason, output.
# A reason of `null` means no reason. A value with a double quote, backslash or
# control character is refused, because the entry is not JSON-escaped. The entry
# names the saved tool command as <NAME>/command.sh, relative to the stage
# directory.
set -euo pipefail

if [[ $# -ne 4 ]]; then
    echo 'Usage: tool-entry.sh NAME REPORT_TSV DATABASE_NAME DATABASE_VERSION' >&2
    exit 2
fi
name=$1
report=$2
database_name=$3
database_version=$4

IFS=$'\t' read -r version status reason output < "$report"
[[ -n $version && -n $status && -n $reason && -n $output ]] || {
    printf 'Incomplete tool report: %s\n' "$report" >&2
    exit 1
}
[[ $status == completed || $status == skipped ]] || {
    printf 'Unexpected tool status %s in %s\n' "$status" "$report" >&2
    exit 1
}
# The values go into JSON unescaped, so refuse any that would break it.
for field in "$version" "$status" "$reason" "$output" "$database_name" "$database_version"; do
    [[ $field != *[\"\\[:cntrl:]]* ]] || {
        printf 'Tool entry value must not contain a double quote, backslash or control character: %q\n' "$field" >&2
        exit 1
    }
done
if [[ $reason == null ]]; then
    reason_json=null
else
    reason_json=\"$reason\"
fi
printf '"%s": {"version": "%s", "status": "%s", "reason": %s, "database_name": "%s", "database_version": "%s", "output": "%s", "command_file": "%s/command.sh"}' \
    "$name" "$version" "$status" "$reason_json" "$database_name" "$database_version" \
    "$output" "$name"
