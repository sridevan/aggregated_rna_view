#!/usr/bin/env bash
# Local Arpeggio proof of concept (pdbe-arpeggio 1.4.4, openbabel-wheel 3.1.1.23, gemmi 0.7.5; macOS arm64, Python 3.11).
# Inputs: PDBe updated mmCIF files already in the project cache (cache/www_ebi_ac_uk/pdbe_entry-files_download_<pdb>_updated.cif-*.cif).
set -eu; ROOT="$(cd "$(dirname "$0")/../../.." && pwd)"; W="${1:-/tmp/arpeggio_poc}"; mkdir -p "$W"; cd "$W"
python3 -m venv venv && ./venv/bin/pip install -q pdbe-arpeggio openbabel-wheel gemmi
cp "$ROOT"/cache/www_ebi_ac_uk/pdbe_entry-files_download_1y26_updated.cif-*.cif 1y26.cif
cp "$ROOT"/cache/www_ebi_ac_uk/pdbe_entry-files_download_3rkf_updated.cif-*.cif 3rkf.cif
time ./venv/bin/pdbe-arpeggio 1y26.cif -s /X/90/ -o out_1y26                       # 0.46 s wall, 67 MB
time ./venv/bin/pdbe-arpeggio 3rkf.cif -s /A/91/ /B/91/ /C/91/ /D/91/ -o out_3rkf  # 1.4 s wall, 107 MB
# outputs: out_<pdb>/<pdb>.json, copied to ../cache/arpeggio_local_<pdb>.json for compare_contacts.py
