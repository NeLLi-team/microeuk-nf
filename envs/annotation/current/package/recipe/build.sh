#!/usr/bin/env bash
set -euo pipefail

source_root="${INTERPROSCAN_SOURCE_ROOT:-/clusterfs/jgi/scratch/science/mgs/nelli/databases/interproscan}"
share_dir="${PREFIX}/share/interproscan-5.76-107.0"
bin_dir="${PREFIX}/bin"

cd "${source_root}"
sha256sum --check "${RECIPE_DIR}/source.sha256"
sha256sum --check "${RECIPE_DIR}/source-static.sha256"
grep -Fx 'version: 5.76-107.0' .nellidb_version/version.txt >/dev/null
grep -Fx 'md5: 8c9a8b153e527f8cfc7bf24ee1652d78' \
  .nellidb_version/version.txt >/dev/null
grep -F 'interproscan-management-5.76-107.0.jar' \
  < <(unzip -p interproscan-5.jar META-INF/MANIFEST.MF) >/dev/null

install -d "${share_dir}/.nellidb_version" "${bin_dir}"
cp -a bin lib "${share_dir}/"
install -d "${share_dir}/work/template" "${share_dir}/work/kvs"
cp -a work/freemarker "${share_dir}/work/"
install -m 0644 work/template/interpro.zip \
  "${share_dir}/work/template/interpro.zip"
cp -a work/kvs/idb "${share_dir}/work/kvs/"
# LOCK is runtime state. Keep the remaining LevelDB files together as the
# pinned distribution seed.
rm "${share_dir}/work/kvs/idb/entryDB/LOCK"
install -m 0644 interproscan-5.jar interproscan.properties "${share_dir}/"
install -m 0755 interproscan.sh "${share_dir}/interproscan.sh"
install -m 0644 .nellidb_version/version.txt \
  "${share_dir}/.nellidb_version/version.txt"
install -m 0644 "${RECIPE_DIR}/source.sha256" "${share_dir}/source.sha256"
install -m 0644 "${RECIPE_DIR}/source-static.sha256" \
  "${share_dir}/source-static.sha256"
install -m 0755 "${RECIPE_DIR}/interproscan.sh" "${bin_dir}/interproscan.sh"
install -m 0755 "${RECIPE_DIR}/interproscan-current-preflight" \
  "${bin_dir}/interproscan-current-preflight"

test ! -e "${share_dir}/data"
test ! -e "${share_dir}/temp"
test ! -e "${share_dir}/work/idb"
test -d "${share_dir}/work/freemarker/WEB-INF/freemarker"
test -s "${share_dir}/work/template/interpro.zip"
test -s "${share_dir}/work/kvs/idb/entries.json"
test ! -e "${share_dir}/work/kvs/idb/entryDB/LOCK"

cd "${share_dir}"
{
  find .nellidb_version bin lib work \( -type f -o -type l \) -print0
  printf '%s\0' interproscan-5.jar interproscan.properties interproscan.sh \
    source.sha256 source-static.sha256
} | LC_ALL=C sort -z | xargs -0 sha256sum > copied-files.sha256
