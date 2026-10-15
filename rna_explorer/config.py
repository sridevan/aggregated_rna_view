"""Paths and external tool discovery."""
from __future__ import annotations
import os, shutil
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_CACHE_DIR = PROJECT_ROOT / "cache"
DEFAULT_OUT_DIR = PROJECT_ROOT / "data"
BUNDLED_CMALIGN = PROJECT_ROOT / "tools" / "infernal-1.1.5" / "bin" / "cmalign"

RFAM_BASE = "https://rfam.org"
RFAM_FTP = "https://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT"
PDBE_API = "https://www.ebi.ac.uk/pdbe/api"
PDBE_SEARCH = "https://www.ebi.ac.uk/pdbe/search/pdb/select"
PDBE_FILES = "https://www.ebi.ac.uk/pdbe/entry-files/download"


def find_cmalign() -> Path:
    """Locate cmalign: env RNA_EXPLORER_CMALIGN, bundled build, then PATH."""
    env = os.environ.get("RNA_EXPLORER_CMALIGN")
    if env and Path(env).is_file():
        return Path(env)
    if BUNDLED_CMALIGN.is_file():
        return BUNDLED_CMALIGN
    found = shutil.which("cmalign")
    if found:
        return Path(found)
    raise FileNotFoundError("cmalign not found: set RNA_EXPLORER_CMALIGN or install Infernal 1.1.x")
