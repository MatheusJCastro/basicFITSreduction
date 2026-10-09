"""Tools for reducing raw astronomical FITS images."""

from __future__ import annotations

import re
from collections import defaultdict
from pathlib import Path
from typing import Callable, Iterable

import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from astropy.visualization import PowerDistStretch, ZScaleInterval
from scipy.ndimage import shift
from scipy.spatial import cKDTree
from skimage.feature import blob_log
from skimage.registration import phase_cross_correlation

__version__ = "0.3.2"


class ReduceRaw:
    """Reduce, calibrate, align, and combine astronomical FITS images."""

    def __init__(self) -> None:
        """Initialize the reduction configuration and result containers."""
        self.dir_name: str | Path = ""
        self.obj_name = ""
        self.combine_filters: list[str] = ["R", "G", "B"]
        self.combine = True
        self.reduce_bias = True
        self.reduce_flat = True
        self.show_blobs = False

        self.data_calibrated: dict[str, np.ndarray] = {}
        self.data_uncalibrated: dict[str, np.ndarray] = {}
        self.data_aligned: dict[str, np.ndarray] = {}

    def run(self) -> None:
        """Run the configured reduction pipeline and save its outputs."""
        self._validate_configuration()
        files = self._discover_files()
        calibration = self._load_calibration_data(files)
        data_by_filter = self._group_data_files(files["data"])

        self._reduce_images(calibration, data_by_filter)
        self._save_calibrated_images(data_by_filter)

        if self.combine:
            self.data_aligned = self._align(
                {
                    name: self.data_calibrated[name]
                    for name in self.combine_filters
                }
            )
            self._save_combined_images()

    def _validate_configuration(self) -> None:
        """Validate the minimum configuration required by :meth:`run`."""
        if not self.dir_name:
            raise ValueError("dir_name must identify a data directory")
        if not self.obj_name:
            raise ValueError("obj_name must identify the target object")
        if self.combine and not self.combine_filters:
            raise ValueError("combine_filters must contain at least one filter")

        data_directory = Path(self.dir_name)
        if not data_directory.is_dir():
            raise FileNotFoundError(f"Data directory not found: {data_directory}")

    def _discover_files(self) -> dict[str, list[Path]]:
        """Find calibration and science FITS files.

        Returns:
            A mapping with ``bias``, ``flat``, ``astrometry``, and ``data`` files.
        """
        data_directory = Path(self.dir_name)
        calibration_directory = data_directory.parent
        all_files = [
            file
            for file in calibration_directory.iterdir()
            if file.is_file() and file.suffix.lower() == ".fits"
        ]
        all_files.extend(
            file
            for file in data_directory.iterdir()
            if file.is_file() and file.suffix.lower() == ".fits"
        )

        files = {
            "bias": [],
            "flat": [],
            "astrometry": [],
            "data": [],
        }
        object_name = self.obj_name.lower()

        for file in all_files:
            filename = file.name.lower()
            if "bias" in filename:
                files["bias"].append(file)
            elif "flat" in filename:
                files["flat"].append(file)
            elif "astrometry" in filename or "field" in filename:
                files["astrometry"].append(file)
            elif object_name in file.stem.lower():
                files["data"].append(file)

        return files

    @staticmethod
    def _filter_name(path: Path, object_name: str) -> str:
        """Extract a normalized filter name from a FITS filename."""
        name = path.stem.upper()
        name = re.sub(object_name, "", name, flags=re.IGNORECASE)
        name = re.sub(r"\d+", "", name)
        name = re.sub(r"[^A-Z]+", " ", name).strip()
        return name

    def _group_data_files(self, files: Iterable[Path]) -> dict[str, list[Path]]:
        """Group science FITS files by filter name."""
        grouped = defaultdict(list)
        for file in files:
            filter_name = self._filter_name(file, self.obj_name)
            if filter_name:
                grouped[filter_name.upper()].append(file)
        return dict(grouped)

    def _group_flat_files(self, files: Iterable[Path]) -> dict[str, list[Path]]:
        """Group flat-field FITS files by filter name."""
        grouped = defaultdict(list)
        for file in files:
            filter_name = self._filter_name(file, "flat")
            if filter_name:
                grouped[filter_name.upper()].append(file)
        return dict(grouped)

    @staticmethod
    def _read_fits(files: Iterable[Path]) -> np.ndarray:
        """Read the primary image from each FITS file as float64 data."""
        images = []
        for file in files:
            with fits.open(file) as hdul:
                images.append(hdul[0].data.astype(np.float64))
        return np.asarray(images)

    def _load_calibration_data(self, files: dict[str, list[Path]]) -> dict[str, object]:
        """Load bias and flat data needed by the reduction pipeline."""
        flat_by_filter = self._group_flat_files(files["flat"])
        bias_images = self._read_fits(files["bias"]) if files["bias"] else None
        flat_images = {
            name: self._read_fits(paths)
            for name, paths in flat_by_filter.items()
        }
        return {
            "bias": bias_images,
            "flats": flat_images,
            "flat_files": flat_by_filter,
        }

    def _reduce_images(self, calibration: dict[str, object], data_by_filter: dict[str, list[Path]]) -> None:
        """Create calibrated and uncalibrated median images by filter."""
        bias_images = calibration["bias"]
        flat_images = calibration["flats"]
        master_bias = None
        if self.reduce_bias:
            if bias_images is None or not len(bias_images):
                raise ValueError("Bias reduction requested, but no bias files exist")
            master_bias = np.median(bias_images, axis=0)

        for filter_name, files in data_by_filter.items():
            data_images = self._read_fits(files)
            self.data_uncalibrated[filter_name] = np.median(data_images, axis=0)
            reduced_images = data_images

            if master_bias is not None:
                reduced_images = reduced_images - master_bias

            if self.reduce_flat:
                if filter_name not in flat_images:
                    raise ValueError(
                        f"Flat reduction requested, but no flat exists for "
                        f"filter {filter_name}"
                    )
                flat_corrected = np.median(
                    flat_images[filter_name]
                    - (master_bias if master_bias is not None else 0),
                    axis=0,
                )
                master_flat = flat_corrected / np.median(flat_corrected)
                reduced_images = reduced_images / master_flat

            self.data_calibrated[filter_name] = np.median(reduced_images, axis=0)

    def _save_fits(self, input_file: Path, image: np.ndarray, output_file: Path) -> None:
        """Save an image while preserving its FITS header."""
        with fits.open(input_file) as hdul:
            hdul[0].data = image.astype(np.float64)
            hdul[0].header["CALIB"] = (
                True,
                "Image has been calibrated",
            )
            hdul[0].header["BIASCOR"] = (
                self.reduce_bias,
                "Master bias subtraction applied",
            )
            hdul[0].header["FLATCOR"] = (
                self.reduce_flat,
                "Master flat correction applied",
            )
            hdul.writeto(output_file, overwrite=True)

    def _save_calibrated_images(self, data_by_filter: dict[str, list[Path]]) -> None:
        """Save one calibrated FITS file for every available filter."""
        output_directory = Path(self.dir_name) / "calibrated"
        output_directory.mkdir(exist_ok=True)

        for filter_name, image in self.data_calibrated.items():
            print(f"Saving {self.obj_name} {filter_name}")
            self._save_fits(
                data_by_filter[filter_name][0],
                image,
                output_directory
                / f"{self.obj_name.upper()}_{filter_name}.fits")


    @staticmethod
    def _normalize(image: np.ndarray, low: float = 1, high: float = 99) -> np.ndarray:
        """Normalize an image to the range from zero to one."""
        vmin, vmax = np.percentile(image, [low, high])
        if vmax == vmin:
            return np.zeros_like(image, dtype=np.float64)
        return np.clip((image - vmin) / (vmax - vmin), 0, 1)

    @staticmethod
    def _zscale(image: np.ndarray) -> np.ndarray:
        """Normalize an image using an astronomical z-scale stretch."""
        interval = ZScaleInterval()
        finite = np.isfinite(image)
        if not finite.any():
            return np.zeros_like(image, dtype=np.float64)

        background = np.nanpercentile(image[finite], 20)
        scale = np.nanpercentile(image[finite], 99) - background
        if not np.isfinite(scale) or scale <= 0:
            return np.zeros_like(image, dtype=np.float64)

        normalized = np.nan_to_num((image - background) / scale, nan=0.0)
        normalized = np.maximum(normalized, 0.0)
        vmin, vmax = interval.get_limits(normalized)
        if vmax == vmin:
            return np.zeros_like(normalized, dtype=np.float64)
        normalized = np.clip(normalized, vmin, vmax)
        normalized = (normalized - vmin) / (vmax - vmin)
        return PowerDistStretch()(normalized)

    def _create_rgb(self, data: dict[str, np.ndarray], method: Callable[[np.ndarray], np.ndarray]) -> np.ndarray:
        """Create an RGB array from the configured filter images."""
        channels = [method(data[name.upper()]) for name in self.combine_filters]
        return np.stack(channels, axis=-1)

    def _find_stars(self, image: np.ndarray) -> np.ndarray:
        """Detect star-like blobs and return their row-column coordinates."""
        image = np.asarray(image, dtype=np.float64)
        finite = np.isfinite(image)
        if not finite.any():
            return np.empty((0, 2), dtype=np.float64)

        background = np.nanpercentile(image[finite], 20)
        scale = np.nanpercentile(image[finite], 99) - background
        if not np.isfinite(scale) or scale <= 0:
            return np.empty((0, 2), dtype=np.float64)

        normalized = np.nan_to_num((image - background) / scale, nan=0.0)
        detections = blob_log(
            np.maximum(normalized, 0.0),
            min_sigma=1,
            max_sigma=4,
            num_sigma=10,
            threshold=3,
            overlap=0.5,
        )

        if self.show_blobs:
            figure, axis = plt.subplots()
            axis.imshow(normalized, cmap="gray", origin="upper")
            for row, column, sigma in detections:
                axis.add_patch(
                    plt.Circle(
                        (column, row),
                        radius=np.sqrt(2) * sigma,
                        fill=False,
                        color="red",
                    )
                )
            axis.set_title("Detected stars")
            plt.show()
            plt.close(figure)

        return (
            detections[:, :2]
            if len(detections)
            else np.empty((0, 2), dtype=np.float64)
        )

    def _find_shift(self, reference: np.ndarray, image: np.ndarray) -> tuple[float, float]:
        """Estimate the row and column shift needed to align an image."""
        reference_stars = self._find_stars(reference)
        image_stars = self._find_stars(image)
        if len(reference_stars) < 2 or len(image_stars) < 2:
            raise ValueError(
                "Unable to detect enough stars to align the images"
            )

        coarse_shift, _, _ = phase_cross_correlation(
            reference,
            image,
            normalization="phase",
            upsample_factor=10,
        )
        shifted_stars = image_stars + coarse_shift
        distances, indices = cKDTree(shifted_stars).query(
            reference_stars,
            distance_upper_bound=10,
        )
        valid = np.isfinite(distances) & (indices < len(image_stars))
        if valid.sum() < 2:
            return tuple(coarse_shift)

        unique_matches = {}
        for distance, reference_star, image_index in zip(
            distances[valid],
            reference_stars[valid],
            indices[valid],
        ):
            if (
                image_index not in unique_matches
                or distance < unique_matches[image_index][0]
            ):
                unique_matches[image_index] = (distance, reference_star)

        if len(unique_matches) < 2:
            return tuple(coarse_shift)

        matched_reference = np.array(
            [match[1] for match in unique_matches.values()]
        )
        matched_image = image_stars[list(unique_matches.keys())]
        shifts = matched_reference - matched_image
        estimated_shift = np.median(shifts, axis=0)
        residuals = np.linalg.norm(shifts - estimated_shift, axis=1)
        limit = max(2.0, np.median(residuals) * 2.5)
        inliers = residuals <= limit
        if inliers.sum() >= 2:
            estimated_shift = np.median(shifts[inliers], axis=0)

        return tuple(estimated_shift)

    def _align(self, data: dict[str, np.ndarray]) -> dict[str, np.ndarray]:
        """Align all filter images to the first configured filter."""
        aligned = {}
        reference = None
        for index, filter_name in enumerate(data):
            if index == 0:
                reference = data[filter_name].copy()
                aligned[filter_name] = reference
                print(f"Base alignment filter is {filter_name}")
                continue

            dy, dx = self._find_shift(reference, data[filter_name])
            print(
                f"Alignment for filter {filter_name} is "
                f"dx = {dx}, dy = {dy}")
            aligned[filter_name] = shift(
                data[filter_name],
                shift=(dy, dx),
                order=3,
                mode="constant",
                cval=0)

        return aligned

    def _save_combined_images(self) -> None:
        """Create and save unaligned, aligned, and z-scaled RGB images."""
        output_directory = Path(self.dir_name) / "calibrated"
        filters = "".join(self.combine_filters)
        prefix = self.obj_name.upper() + "_" + filters

        images = {
            prefix + "_unaligned.pdf": self._create_rgb(
                self.data_calibrated, self._normalize
            ),
            prefix + "_unaligned_nocalib.pdf": self._create_rgb(
                self.data_uncalibrated, self._normalize
            ),
            prefix + ".pdf": self._create_rgb(
                self.data_aligned, self._normalize
            ),
            prefix + "_unaligned_zscale.pdf": self._create_rgb(
                self.data_calibrated, self._zscale
            ),
            prefix + "_unaligned_zscale_nocalib.pdf": self._create_rgb(
                self.data_uncalibrated, self._zscale
            ),
            prefix + "_zscale.pdf": self._create_rgb(
                self.data_aligned, self._zscale
            ),
        }
        print("Saving RGB images")
        for filename, image in images.items():
            plt.imsave(output_directory / filename, image)
        print("Done")


__all__ = ["ReduceRaw", "__version__"]
