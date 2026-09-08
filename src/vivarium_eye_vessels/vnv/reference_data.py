"""Download and cache expert-labeled vessel masks from public datasets.

Currently supports the healthy subset of the HRF (High-Resolution Fundus)
image database: 15 fundus photographs of healthy eyes with binary vessel
segmentations hand-labeled by experts, and the artery/vein labels of the
same 15 eyes.

    Budai A, Bock R, Maier A, Hornegger J, Michelson G. Robust Vessel
    Segmentation in Fundus Images. Int J Biomed Imaging, 2013.
    https://www5.cs.fau.de/research/data/fundus-images/

    Hemelings R, Elen B, Stalmans I, Van Keer K, De Boever P, Blaschko MB.
    Artery-vein segmentation in fundus images using a fully convolutional
    network. Comput Med Imaging Graph, 2019.
    https://github.com/rubenhx/av-segmentation

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


ROSE_URL = "https://imed.nimte.ac.cn/dataofrose.html"


def fetch_rose_images(projection: str = "SVC") -> list[Path]:
    """Paths of the ROSE-1 en-face angiograms (train and test) for one projection.

    ROSE (Retinal OCT-Angiography vessel SEgmentation) is released on
    request, so nothing is downloaded: extract ROSE.zip under the cache
    directory so that ``<cache>/rose/ROSE/ROSE-1/<projection>/...`` exists.
    ``projection`` is SVC (superficial vascular complex), DVC or SVC_DVC.
    """
    return rose_files(projection, "img")


def rose_files(projection: str, kind: str) -> list[Path]:
    """Sorted ROSE-1 files of one ``kind`` (``img``, ``gt``, ``thick_gt``, ``thin_gt``)."""
    root = get_cache_dir() / "rose" / "ROSE" / "ROSE-1" / projection
    paths = sorted(
        path
        for split in ("train", "test")
        for path in (root / split / kind).glob("*")
        if path.suffix.lower() in (".tif", ".tiff", ".png")
    )
    if not paths:
        raise FileNotFoundError(
            f"No ROSE-1 {projection} {kind} files under {root}. Request the dataset at "
            f"{ROSE_URL} and extract ROSE.zip into {get_cache_dir() / 'rose'}."
        )
    return paths


def fetch_rose_labels(projection: str = "SVC", kind: str = "gt") -> list[Path]:
    """Paths of the ROSE-1 expert labels, paired by name with :func:`fetch_rose_images`.

    ``gt`` is the full pixel-level label (large vessels filled, capillaries
    as centerlines); ``thick_gt`` and ``thin_gt`` split it by class.
    """
    return rose_files(projection, kind)


def load_mask(path: Path) -> np.ndarray:
    """Load a vessel mask as a binary-ish uint8 array."""
    with Image.open(path) as image:
        return np.array(image.convert("L"))


HRF_AV_URL = "https://raw.githubusercontent.com/rubenhx/av-segmentation/master/HRF_AV_GT/"


def fetch_hrf_av_labels() -> list[Path]:
    """Download (if needed) and return the artery/vein labels of the HRF healthy eyes.

    Hemelings et al. labeled every vessel pixel of the 45 HRF images as
    artery or vein (pixels where the two cross in a third class) and
    published the labels as supplementary material to their 2019 paper (see
    the module docstring). The 15 healthy eyes' labels partition the HRF
    healthy masks pixel for pixel -- same images, same vessel pixels -- so
    the artery/vein targets are read on the very eyes the other fundus
    targets come from. Cached beside the masks; label ``NN_h_AVmanual.png``
    pairs with mask ``NN_h.tif``.
    """
    cache = get_cache_dir() / "hrf_av_gt"
    labels = [cache / f"{number:02d}_h_AVmanual.png" for number in range(1, 16)]
    if all(path.exists() for path in labels):
        return labels
    cache.mkdir(parents=True, exist_ok=True)
    print(f"Downloading HRF artery/vein labels to {cache} ...")
    for path in labels:
        if not path.exists():
            urllib.request.urlretrieve(HRF_AV_URL + path.name, path)
    return labels


def load_av_label(path: Path) -> tuple[np.ndarray, np.ndarray]:
    """An artery/vein label as two boolean masks, (artery, vein).

    The labels color arteries red, veins blue and the pixels where the two
    cross green; a crossing pixel belongs to both vessels.
    """
    with Image.open(path) as image:
        rgb = np.asarray(image.convert("RGB")) > 127
    red, green, blue = rgb[..., 0], rgb[..., 1], rgb[..., 2]
    return red | green, blue | green
