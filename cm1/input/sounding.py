"""
Load ERA5 model Sounding for a user-specified time and location.

--campaign: Use campaign storage.
Otherwise use the s3fs Amazon Web Service bucket or a local cached file.
"""

import argparse
import io
import logging
import os
from pathlib import Path
from typing import TextIO, Union

import metpy.calc as mcalc
import metpy.constants
import numpy as np
import pandas as pd
import xarray as xr
from metpy.units import units
from pint import Quantity

import cm1.input.era5
from cm1.utils import TMPDIR, parse_args

# Assuming this script is located in a subdirectory of the repository
repo_base_path = os.path.abspath(os.path.join(os.path.dirname(__file__), "../.."))
soundings_path = os.path.join(repo_base_path, "soundings")


class Sounding(xr.Dataset):
    __slots__ = ()

    def __init__(self, data_or_path=None, *args, **kwargs):
        """
        Initializes the Sounding object.

        This constructor is designed to be compatible with xarray's internal
        mechanisms (like .map, .sel) which may pass a dict of variables,
        while also handling user-provided file paths or xarray.Dataset objects.
        """
        # This is the path xarray's internal methods (like .map) will take.
        # They pass a dictionary of variables as the first argument.
        if isinstance(data_or_path, dict):
            super().__init__(data_or_path, *args, **kwargs)
            return

        # This block handles user-initiated creation from a path or Dataset.
        ds = None
        if isinstance(data_or_path, (str, Path, os.PathLike)):
            ds = xr.open_dataset(data_or_path)
        else:
            # Assumes an xr.Dataset, xr.DataArray, or None
            ds = data_or_path

        if isinstance(ds, (xr.Dataset, xr.DataArray)):
            # We are creating a Sounding from another xarray object. Quantify it.
            quantified_ds = ds.metpy.quantify()
            # Explicitly pass the components to the parent constructor for robustness.
            super().__init__(quantified_ds, *args, **kwargs)
        else:
            # Handles the None case, passing it up to the parent.
            super().__init__(ds, *args, **kwargs)

    @classmethod
    def from_txt(cls, source: Union[str, os.PathLike]):
        """
        Creates a Sounding instance by parsing a special CM1 text file.

        This method can handle:
        1. A file path (as a string or Path object).
        2. The raw string content of a sounding file (e.g., from to_txt()).
        """
        # If the source is a string and contains newlines, treat it as content.
        if isinstance(source, str) and "\n" in source:
            string_stream = io.StringIO(source)
            # For content, there's no file name, so provide a default hint.
            return cls._parse_stream(string_stream, case_name_hint="from_string")

        # Otherwise, treat it as a file path.
        path = Path(source)
        if not path.is_file():
            raise FileNotFoundError(f"Sounding file not found at path: {path}")

        with path.open("r") as f:
            return cls._parse_stream(f, case_name_hint=path.stem)

    @classmethod
    def _parse_stream(cls, stream: TextIO, case_name_hint: str):
        """Helper method to parse a CM1 text file from a text stream."""
        header = stream.readline().strip()
        surface_pressure, surface_theta, surface_mixing_ratio = map(
            float, header.split()
        )

        column_names = ["Z", "theta", "qv", "U", "V"]
        df = pd.read_csv(stream, sep=r"\s+", names=column_names, engine="python")

        df = df.rename_axis("level")
        ds = df.to_xarray()

        # Add units and calculate pressure
        ds["qv"] *= units("g/kg")
        ds["qv"].attrs["long_name"] = "water vapor mixing ratio"
        ds["Q"] = mcalc.specific_humidity_from_mixing_ratio(ds["qv"])

        ds["SP"] = surface_pressure
        ds["SP"] *= units.hPa

        ds["theta"] *= units.K
        ds["Z"] *= units.m
        ds["Z"].attrs["long_name"] = "geopotential height"
        p_bot = ds.SP.copy()
        z_bot = 0.0 * units.m

        P = np.empty_like(ds.level, dtype=float)
        for i, level in enumerate(ds.level):
            T = mcalc.temperature_from_potential_temperature(
                p_bot, ds.theta.sel(level=level)
            )
            Tv = mcalc.virtual_temperature(T, ds.qv.sel(level=level))
            dz = ds.Z.sel(level=level) - z_bot
            p_bot = p_bot * np.exp(-metpy.constants.g * dz / (metpy.constants.Rd * Tv))
            assert (
                p_bot >= 0 * units.hPa
            ), f"{case_name_hint} p_bot<0 {p_bot:~} dz {dz:~} Tv {Tv:~}"
            P[i] = p_bot.data.to(ds.SP.metpy.units).m.item()
            z_bot = ds.Z.sel(level=level)

        ds["P"] = ("level", P)
        ds["P"] *= ds.SP.metpy.units
        ds["T"] = mcalc.temperature_from_potential_temperature(ds["P"], ds["theta"])
        ds["Tv"] = mcalc.virtual_temperature(ds.T, ds.qv)

        # Add surface data
        ds["surface_potential_temperature"] = surface_theta
        ds["surface_potential_temperature"] *= units.K
        ds["surface_mixing_ratio"] = surface_mixing_ratio
        ds["surface_mixing_ratio"] *= units.g / units.kg
        ds["surface_geopotential_height"] = 0.0
        ds["surface_geopotential_height"] *= units.m
        ds["U"] *= units.m / units.s
        ds["V"] *= units.m / units.s

        case_name = case_name_hint.replace("input_sounding_", "")
        sounding_obj = cls(ds)
        sounding_obj.attrs["case"] = case_name
        return sounding_obj

    @classmethod
    def get_case(cls, case: str):
        """Retrieves a predefined sounding case dataset."""
        file_path = os.path.join(soundings_path, f"input_sounding_{case}")
        return cls.from_txt(file_path)

    def plot(self, fig=None, subplot=None, **kwargs):
        """Plots the sounding using the skewt function."""
        from cm1.skewt import skewt

        return skewt(self, fig=fig, subplot=subplot, **kwargs)

    def to_txt(self) -> str:
        """Converts a Sounding into a formatted string suitable for CM1."""
        self.load()  # can't use item() on lazy dask array below
        sfc = self.level.sel(level=self.P.idxmax())
        sfc_pres = self.SP if "SP" in self else self.P.sel(level=sfc)
        self["theta"] = mcalc.potential_temperature(self.P, self.T).metpy.convert_units(
            "K"
        )
        sfc_theta_K = (
            self["surface_potential_temperature"].metpy.convert_units("K")
            if "surface_potential_temperature" in self
            else self["theta"].sel(level=sfc)
        )
        if "qv" not in self and "Q" in self:
            self["qv"] = mcalc.mixing_ratio_from_specific_humidity(
                self["Q"]
            ).metpy.convert_units("g/kg")
        sfc_qv_gkg = (
            self["surface_mixing_ratio"].metpy.convert_units("g/kg")
            if "surface_mixing_ratio" in self
            else self.qv.sel(level=sfc)
        )
        header = f"{sfc_pres.item().m_as('hPa'):.2f} {sfc_theta_K.item().m_as('K'):.2f} {sfc_qv_gkg.item().m_as('g/kg'):.2f}\n"
        df_export = self.copy(deep=False)
        for var in ["Z", "theta", "qv", "U", "V"]:
            if isinstance(df_export[var].data, Quantity):
                df_export[var].data = df_export[var].data.m
        body = (
            df_export[["Z", "theta", "qv", "U", "V"]]
            .drop_vars(["latitude", "longitude", "time"], errors="ignore")
            .to_dataframe()
            .sort_values("Z")
            .to_csv(sep=" ", header=False, index=False, float_format="%.2f")
        )
        return header + body


def era5_aws(time: pd.Timestamp, lat: Quantity, lon: Quantity, **kwargs) -> Sounding:
    """Retrieves ERA5 dataset for a specific time and location from AWS."""
    ds = cm1.input.era5.aws(time, **kwargs)
    lon = lon % (360 * units.degreeE)
    ds = ds.sel(longitude=lon, latitude=lat, method="nearest", tolerance=5 * units.deg)
    return Sounding(ds)


def era5_model_level(
    time: pd.Timestamp, lat: Quantity, lon: Quantity, **kwargs
) -> Sounding:
    """Retrieves ERA5 model-level dataset for a specific time and location."""
    ds = cm1.input.era5.model_level(time, **kwargs)
    lon = lon % (360 * units.degreeE)
    ds = ds.sel(longitude=lon, latitude=lat, method="nearest", tolerance=5 * units.deg)
    return Sounding(ds)


def era5_pressure_level(
    time: pd.Timestamp, lat: Quantity, lon: Quantity, **kwargs
) -> Sounding:
    """Retrievis ERA5 pressure-level dataset for a specific time and location."""
    ds = cm1.input.era5.pressure_level(time, **kwargs)
    lon = lon % (360 * units.degreeE)
    ds = ds.sel(longitude=lon, latitude=lat, method="nearest", tolerance=5 * units.deg)
    return Sounding(ds)


def get_ofile(args: argparse.Namespace) -> Path:
    """Generates a temporary file path for caching the dataset."""
    time_str = pd.to_datetime(args.time).strftime("%Y%m%d_%H%M%S")
    lat_str = f"{args.lat.m:~}"
    lon_str = f"{args.lon.m:~}"
    return TMPDIR / f"{time_str}.{lat_str}.{lon_str}"


# Functions for specific sounding cases now use the class method
def trier():
    return Sounding.get_case("trier")


def jordan_allmean():
    return Sounding.get_case("jordan_allmean")


def jordan_hurricane():
    return Sounding.get_case("jordan_hurricane")


def rotunno_emanuel():
    return Sounding.get_case("rotunno_emanuel")


def dunion_MT():
    return Sounding.get_case("dunion_MT")


def bryan_morrison():
    return Sounding.get_case("bryan_morrison")


def seabreeze_test():
    return Sounding.get_case("seabreeze_test")


def main() -> None:
    """Main function for loading ERA5 data and printing sounding data."""
    import pickle

    args = parse_args()
    valid_time = pd.to_datetime(args.time)
    ofile = get_ofile(args)

    if os.path.exists(ofile):
        logging.warning(f"Reading from cache: {ofile}")
        with open(ofile, "rb") as file:
            ds = pickle.load(file)
    else:
        if os.path.exists("/glade/campaign"):
            ds = era5_model_level(valid_time, args.lat, args.lon)
        else:
            logging.warning(
                "No campaign storage. Getting pressure level data from AWS."
            )
            ds = era5_aws(valid_time, args.lat, args.lon)

        with open(ofile, "wb") as file:
            logging.warning(f"Caching data to: {ofile}")
            pickle.dump(ds, file)

    print(ds.to_txt())


if __name__ == "__main__":
    main()
