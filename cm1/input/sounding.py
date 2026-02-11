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
        if isinstance(data_or_path, dict):
            super().__init__(data_or_path, *args, **kwargs)
            return

        ds = None
        if isinstance(data_or_path, (str, Path, os.PathLike)):
            ds = xr.open_dataset(data_or_path)
        else:
            ds = data_or_path

        if isinstance(ds, (xr.Dataset, xr.DataArray)):
            quantified_ds = ds.metpy.quantify()
            super().__init__(quantified_ds, *args, **kwargs)
        else:
            super().__init__(ds, *args, **kwargs)

    def __setitem__(self, key, value):
        """
        Overridden to handle assignment of scalar Pint Quantities, which
        otherwise cause AttributeErrors in xarray because they lack a 'shape' attribute.
        """
        # Check if it looks like a Pint Quantity (has units and magnitude)
        if hasattr(value, "units") and hasattr(value, "magnitude"):
            # If the magnitude is a scalar (float/int) without a shape attribute
            if not hasattr(value.magnitude, "shape"):
                # Wrap the magnitude in a 0-d numpy array to satisfy xarray
                value = np.array(value.magnitude) * value.units

        super().__setitem__(key, value)

    @classmethod
    def from_txt(cls, source: Union[str, os.PathLike]):
        """
        Creates a Sounding instance by parsing a special CM1 text file.
        """
        if isinstance(source, str) and "\n" in source:
            string_stream = io.StringIO(source)
            return cls._parse_cm1_txt_stream(
                string_stream, case_name_hint="from_string"
            )

        path = Path(source)
        if not path.is_file():
            raise FileNotFoundError(f"Sounding file not found at path: {path}")

        with path.open("r") as f:
            return cls._parse_cm1_txt_stream(f, case_name_hint=path.stem)

    @classmethod
    def _integrate_pressure(cls, Z, theta, qv, sfc_pres, sfc_theta, sfc_mix):
        """
        Integrates the hydrostatic equation upward to calculate pressure at all levels.
        Based on the physics in George Bryan's getcape.f90.
        """
        # Constants
        g = metpy.constants.g
        Cp = metpy.constants.Cp_d
        kappa = metpy.constants.kappa  # Rd/Cp
        P0 = 1000.0 * units.hPa

        # Don't worry about swapping potential_temperature for temperature.
        # Both multiplied by same factor below to get "virtual" version.
        # (w+e)/(e(1+w))
        theta_v = mcalc.virtual_temperature(theta, qv)

        # Initialize boundary conditions with surface values
        z_prev = 0.0 * units.m
        pi_prev = (sfc_pres / P0) ** kappa

        # Surface virtual potential temperature
        theta_v_prev = mcalc.virtual_temperature(sfc_theta, sfc_mix)

        P_values = []

        for i, z_curr in enumerate(Z):
            theta_v_curr = theta_v[i]
            dz = z_curr - z_prev

            # Layer mean virtual potential temperature
            theta_v_mean = 0.5 * (theta_v_prev + theta_v_curr)

            # Change in Exner function: dpi = - (g / (Cp * theta_v_mean)) * dz
            dpi = -g * dz / (Cp * theta_v_mean)

            pi_curr = pi_prev + dpi

            # Convert back to Pressure: P = P0 * pi^(1/kappa)
            p_curr = P0 * (pi_curr ** (1.0 / kappa))
            P_values.append(p_curr.m_as(sfc_pres.units))

            # Update for next level
            pi_prev = pi_curr
            theta_v_prev = theta_v_curr
            z_prev = z_curr

        return np.array(P_values) * sfc_pres.units

    @classmethod
    def _parse_cm1_txt_stream(cls, stream: TextIO, case_name_hint: str):
        """Helper method to parse a CM1 text file from a text stream."""
        # Read the header: Surface Pressure (hPa), Surface Potential Temp (K), Surface Mixing Ratio (g/kg)
        header = stream.readline().strip()
        surface_pressure, surface_theta, surface_mixing_ratio = map(
            float, header.split()
        )

        # CM1 sounding format: Height(m), Theta(K), Qv(g/kg), U(m/s), V(m/s)
        column_names = ["Z", "theta", "qv", "U", "V"]
        df = pd.read_csv(stream, sep=r"\s+", names=column_names, engine="python")

        df = df.rename_axis("level")
        ds = df.to_xarray()

        # Apply units and derive initial variables
        ds["qv"] *= units("g/kg")
        ds["qv"].attrs["long_name"] = "water vapor mixing ratio"
        ds["Q"] = mcalc.specific_humidity_from_mixing_ratio(ds["qv"])

        # Convert scalar to numpy array so it has .shape attribute (needed by xarray)
        ds["SP"] = np.array(surface_pressure) * units.hPa
        ds["theta"] *= units.K
        ds["Z"] *= units.m
        ds["Z"].attrs["long_name"] = "geopotential height"

        # --- RECONSTRUCT PRESSURE DATA ---
        # Integrate upward from the surface using the Exner-Hydrostatic equation.
        # We pass the .data property to ensure we work with Pint Quantities directly.
        P = cls._integrate_pressure(
            ds.Z.data,
            ds.theta.data,
            ds.qv.data,
            ds.SP.data,
            surface_theta * units.K,
            surface_mixing_ratio * (units.g / units.kg),
        )

        # Attach calculated pressure and derive dependent temperatures
        ds["P"] = ("level", P)
        ds["T"] = mcalc.temperature_from_potential_temperature(ds["P"], ds["theta"])
        ds["Tv"] = mcalc.virtual_temperature(ds.T, ds.qv)

        # Final metadata and unit cleanup
        ds["surface_potential_temperature"] = np.array(surface_theta) * units.K
        ds["surface_mixing_ratio"] = np.array(surface_mixing_ratio) * (
            units.g / units.kg
        )
        ds["surface_geopotential_height"] = np.array(0.0) * units.m
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
    """Retrieves ERA5 pressure-level dataset for a specific time and location."""
    ds = cm1.input.era5.pressure_level(time, **kwargs)
    lon = lon % (360 * units.degreeE)
    ds = ds.sel(longitude=lon, latitude=lat, method="nearest", tolerance=5 * units.deg)
    return Sounding(ds)


def get_ofile(args: argparse.Namespace) -> Path:
    """Generates a temporary file path for caching the dataset."""
    time_str = pd.to_datetime(args.time).strftime("%Y%m%d_%H%M%S")
    lat_str = f"{args.lat:~}"
    lon_str = f"{args.lon:~}"
    return TMPDIR / f"{time_str}.{lat_str}.{lon_str}"


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
