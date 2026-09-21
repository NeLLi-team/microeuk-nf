#!/usr/bin/env bash
set -euo pipefail

if [[ -z "${INTERPROSCAN_DATA_DIR:-}" ]]; then
  echo 'INTERPROSCAN_DATA_DIR must name the release-matched external data directory.' >&2
  exit 2
fi
if [[ ! -d "${INTERPROSCAN_DATA_DIR}" ]]; then
  echo "InterProScan data directory does not exist: ${INTERPROSCAN_DATA_DIR}" >&2
  exit 2
fi

prefix="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)"
share_dir="${prefix}/share/interproscan-5.76-107.0"
data_dir="$(cd "${INTERPROSCAN_DATA_DIR}" && pwd -P)"
user_dir="${PWD}"
for executable in java perl python3; do
  executable_path="$(command -v "${executable}")"
  if [[ "${executable_path}" != "${prefix}/bin/"* ]]; then
    echo "${executable} must resolve inside the Pixi prefix: ${executable_path}" >&2
    exit 2
  fi
done

cpu="${INTERPROSCAN_CPU:-1}"
arguments=("$@")
for ((index = 0; index < ${#arguments[@]}; index++)); do
  if [[ "${arguments[index]}" == '-cpu' || "${arguments[index]}" == '--cpu' ]]; then
    if ((index + 1 >= ${#arguments[@]})); then
      echo 'InterProScan CPU option requires a value.' >&2
      exit 2
    fi
    cpu="${arguments[index + 1]}"
  fi
done
if [[ ! "${cpu}" =~ ^[1-9][0-9]*$ ]]; then
  echo "InterProScan CPU count must be a positive integer: ${cpu}" >&2
  exit 2
fi
java_initial_heap="${INTERPROSCAN_JAVA_INITIAL_HEAP:-2028M}"
java_max_heap="${INTERPROSCAN_JAVA_MAX_HEAP:-14G}"
if [[ ! "${java_initial_heap}" =~ ^[1-9][0-9]*[MmGg]$ \
  || ! "${java_max_heap}" =~ ^[1-9][0-9]*[MmGg]$ ]]; then
  echo 'InterProScan JVM heap values must use an integer M or G suffix.' >&2
  exit 2
fi

temp_root="${INTERPROSCAN_TEMP_DIR:-${TMPDIR:-${PWD}}}"
if [[ ! -d "${temp_root}" || ! -w "${temp_root}" ]]; then
  echo "InterProScan temporary root is not a writable directory: ${temp_root}" >&2
  exit 2
fi
runtime_dir="$(mktemp -d "${temp_root%/}/interproscan.XXXXXXXX")"
cleanup() {
  find "${runtime_dir}" -depth -delete
}
trap cleanup EXIT

# InterProScan opens its entry LevelDB with a write lock. Give each invocation
# a private copy while keeping the package seed immutable and concurrency-safe.
entry_kv_dir="${runtime_dir}/entry-kv"
cp -a "${share_dir}/work/kvs/idb" "${entry_kv_dir}"

properties="${runtime_dir}/interproscan.properties"
while IFS= read -r line || [[ -n "${line}" ]]; do
  case "${line}" in
    data.directory=*) printf 'data.directory=%s\n' "${data_dir}" ;;
    bin.directory=*) printf 'bin.directory=%s\n' "${share_dir}/bin" ;;
    temporary.file.directory=*)
      printf 'temporary.file.directory=%s\n' "${runtime_dir}/temp/[UNIQUE]"
      ;;
    number.of.embedded.workers=*)
      printf 'number.of.embedded.workers=%s\n' "${cpu}"
      ;;
    maxnumber.of.embedded.workers=*)
      printf 'maxnumber.of.embedded.workers=%s\n' "${cpu}"
      ;;
    worker.number.of.embedded.workers=*)
      printf 'worker.number.of.embedded.workers=%s\n' "${cpu}"
      ;;
    worker.maxnumber.of.embedded.workers=*)
      printf 'worker.maxnumber.of.embedded.workers=%s\n' "${cpu}"
      ;;
    worker.command=*|worker.high.memory.command=*)
      key="${line%%=*}"
      printf '%s=java -XX:ActiveProcessorCount=%s -XX:ParallelGCThreads=%s -Xms%s -Xmx%s -jar %s/interproscan-5.jar\n' \
        "${key}" "${cpu}" "${cpu}" "${java_initial_heap}" "${java_max_heap}" \
        "${share_dir}"
      ;;
    jms.broker.temp.directory=*)
      printf 'jms.broker.temp.directory=%s\n' "${runtime_dir}/activemq-data"
      ;;
    *) printf '%s\n' "${line}" ;;
  esac
done < "${share_dir}/interproscan.properties" > "${properties}"
printf '%s\n' \
  "freemarker.path=${share_dir}/work/freemarker" \
  "i5.h2.database.original.location=${share_dir}/work/template/interpro.zip" \
  "kvstore.entrydb.path=${entry_kv_dir}" >> "${properties}"

java \
  "-XX:ActiveProcessorCount=${cpu}" \
  "-XX:ParallelGCThreads=${cpu}" \
  "-Xms${java_initial_heap}" \
  "-Xmx${java_max_heap}" \
  "-Dsystem.interproscan.properties=${properties}" \
  -jar "${share_dir}/interproscan-5.jar" \
  "$@" -u "${user_dir}"
