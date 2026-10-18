# Isolated fpocketR environment (Apple M3 Pro, macOS 26.6, arm64 host) — created 2026-10-16

Everything is project-local; nothing was installed system-wide, no sudo, no Homebrew/pyenv/shell changes.

    # x86-64 micromamba binary (runs under the already-present Rosetta 2)
    mkdir -p tools/micromamba && cd tools/micromamba
    curl -Ls https://micro.mamba.pm/api/micromamba/osx-64/latest | tar -xj bin/micromamba     # Mach-O x86_64, micromamba 2.9.0
    export MAMBA_ROOT_PREFIX=$PWD/root CONDA_SUBDIR=osx-64
    ./bin/micromamba create -y -p ../fpocketr-env -c conda-forge --platform osx-64 fpocket=4.0.3 pymol-open-source=2.5 python=3.11 pip
    # pinned Python dependencies as conda-forge osx-64 binaries (pip could not build prody 2.3.0 from source under Rosetta)
    ./bin/micromamba install -y -p ../fpocketr-env -c conda-forge --platform osx-64 prody=2.3.0 rdkit=2023.9.6 scipy=1.12.0 trimesh=4.2.0 biopython=1.83 pandas seaborn "setuptools<80"
    ./bin/micromamba run -p ../fpocketr-env pip install --no-deps fpocketr==1.3.6

Verified: `file tools/fpocketr-env/bin/fpocket` → Mach-O 64-bit executable x86_64; Python 3.11.17 reports `platform.machine() == 'x86_64'`;
`python -m fpocketR --help` runs; PyMOL 2.5.0 imports; `pip check` reports no broken requirements.
Pinned versions (see `package_list.txt`, `fpocketr-env.yml`): fpocket 4.0.3 (conda-forge h8a118fb_2), fpocketR 1.3.6 (PyPI), pymol-open-source 2.5.0, prody 2.3.0, rdkit 2023.09.6,
scipy 1.12.0, trimesh 4.2.0, biopython 1.83, numpy 1.26.4, pandas 2.1.4, seaborn 0.13.2, setuptools 79.0.1 (ProDy 2.3.0 imports pkg_resources).
Deviation from the README "conda + pip" route: fpocketR's PyPI dependencies were installed from conda-forge at the exact pinned versions instead of being built by pip; versions are identical to the Bioconda recipe pins.

Run: `tools/micromamba/bin/micromamba run -p tools/fpocketr-env python -W ignore -m fpocketR -pdb <file.pdb> -c <chain> ...`
