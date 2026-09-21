#!/usr/bin/env bash
set -euo pipefail

share_dir="${PREFIX}/share/braker3-longrna"
bin_dir="${PREFIX}/bin"

grep -F 'my $version = "3.0.8";' BRAKER/scripts/braker.pl >/dev/null
grep -F 'system("stringtie -L -p $cores -o $out_gff $in_bam");' ETP/bin/gmetp.pl >/dev/null
grep -F "VERSION = '2.6.0'" ETP/bin/gmes/ProtHint/bin/prothint.py >/dev/null

install -d "${share_dir}/BRAKER" "${share_dir}/ETP" "${share_dir}/licenses" "${bin_dir}"
cp -R BRAKER/scripts "${share_dir}/BRAKER/"
cp BRAKER/LICENSE.TXT "${share_dir}/licenses/BRAKER-LICENSE.TXT"
cp -R ETP/bin ETP/tools "${share_dir}/ETP/"
cp ETP/INSTALL ETP/INSTALL.details ETP/README.md ETP/check_install.pl "${share_dir}/ETP/"
cp ETP/License-Creative-Commons-Attribution-NonCommercial-ShareAlike-4.0-International.txt \
  "${share_dir}/licenses/GeneMark-ETP-LICENSE.txt"
cp ProtHint-license/LICENSE "${share_dir}/licenses/ProtHint-LICENSE.txt"

sed \
  -e "s|@BRAKER_COMMIT@|35280d0c2f56dde026ab346c0255740c11ab80a1|g" \
  -e "s|@ETP_COMMIT@|64a69cbc11e15d51779da021d794cbe9c6d695a7|g" \
  "${RECIPE_DIR}/source.json.in" > "${share_dir}/source.json"

for wrapper in braker.pl gmetp.pl prothint.py braker3-longrna-preflight; do
  install -m 0755 "${RECIPE_DIR}/${wrapper}" "${bin_dir}/${wrapper}"
done

