"""
Load ERA5 model Sounding for a user-specified time and location.

--campaign: Use campaign storage.
Otherwise use the s3fs Amazon Web Service bucket or a local cached file.
"""

import argparse
import logging
import os
import typing
from pathlib import Path

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
    __slots__ = ("case",)

    def __init__(self, data_or_path=None, *args, case=None, **kwargs):
        """
        Initializes the Sounding object from a file path or an xarray.Dataset.
        For special text files or named cases, use the factory methods:
        Sounding.from_txt() or Sounding.get_case().
        """
        ds = None
        if isinstance(data_or_path, (str, Path)):
            ds = xr.open_dataset(data_or_path)
        else:
            ds = data_or_path

        if ds is not None:
            quantified_ds = ds.metpy.quantify()
        else:
            quantified_ds = ds

        super().__init__(quantified_ds, *args, **kwargs)
        self.case = case

    @classmethod
    def from_txt(cls, file_path: typing.Union[str, Path]):
        """Creates a Sounding instance by parsing a special CM1 text file."""
        with open(file_path, "r") as file:
            header = file.readline().strip()
            surface_pressure, surface_theta, surface_mixing_ratio = map(
                float, header.split()
            )
        column_names = ["Z", "theta", "qv", "U", "V"]
        df = pd.read_csv(
            file_path, sep=r"\s+", skiprows=1, names=column_names, engine="python"
        )
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
            ), f"{file_path} p_bot<0 {p_bot:~} dz {dz:~} Tv {Tv:~}"
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

        case_name = Path(file_path).stem.replace("input_sounding_", "")
        return cls(ds, case=case_name)

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
        sfc = self.level.sel(level=self.P.compute().idxmax())
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
    """Retrieves ERA5 pressure-level dataset for a specific time and location."""
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
