"""Download and cache expert-labeled vessel masks from public datasets.

Currently supports the healthy subset of the HRF (High-Resolution Fundus)
image database: 15 fundus photographs of healthy eyes with binary vessel
segmentations hand-labeled by experts.

    Budai A, Bock R, Maier A, Hornegger J, Michelson G. Robust Vessel
    Segmentation in Fundus Images. Int J Biomed Imaging, 2013.
    https://www5.cs.fau.de/research/data/fundus-images/

The dataset is free for research use. Masks are cached locally (default
``~/.cache/vivarium_eye_vessels``, override with the ``VEV_DATA_DIR``
environment variable) and are not committed to this repository.
"""

import os
import urllib.request
import zipfile
from pathlib import Path

import numpy as np
from PIL import Image

HRF_SEGMENTATION_URL = (
    "https://www5.cs.fau.de/fileadmin/research/datasets/fundus-images/"
    "healthy_manualsegm.zip"
)


def get_cache_dir() -> Path:
    root = os.environ.get("VEV_DATA_DIR", "~/.cache/vivarium_eye_vessels")
    return Path(root).expanduser()


def fetch_hrf_masks() -> list[Path]:
    """Download (if needed) and return paths to the HRF healthy vessel masks."""
    cache = get_cache_dir() / "hrf_healthy_manualsegm"
    masks = sorted(cache.glob("*.tif"))
    if masks:
        return masks

    cache.mkdir(parents=True, exist_ok=True)
    archive = cache / "healthy_manualsegm.zip"
    print(f"Downloading HRF healthy vessel masks to {cache} ...")
    urllib.request.urlretrieve(HRF_SEGMENTATION_URL, archive)
    with zipfile.ZipFile(archive) as zf:
        zf.extractall(cache)
    archive.unlink()

    masks = sorted(cache.glob("*.tif"))
    if not masks:
        raise RuntimeError(f"No .tif masks found after extracting HRF archive to {cache}")
    return masks


def load_mask(path: Path) -> np.ndarray:
    """Load a vessel mask as a binary-ish uint8 array."""
    with Image.open(path) as image:
        return np.array(image.convert("L"))
