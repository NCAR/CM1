"""
Load ERA5 model dataset for a user-specified time and location.

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
import xarray
from metpy.units import units
from pint import Quantity

import cm1.input.era5
from cm1.utils import TMPDIR, parse_args

# Assuming this script is located in a subdirectory of the repository
repo_base_path = os.path.abspath(os.path.join(os.path.dirname(__file__), "../.."))
soundings_path = os.path.join(repo_base_path, "soundings")


def era5_aws(time: pd.Timestamp, lat: Quantity, lon: Quantity, **kwargs):
    """
    Retrieve ERA5 dataset for a specific time and location.

    Parameters:
        time (pd.Timestamp): The time for the dataset retrieval.
        lat (float): Latitude of the desired location.
        lon (float): Longitude of the desired location.
        **kwargs: Additional arguments to pass to the data retrieval function.

    Returns:
        xarray.Dataset: ERA5 dataset for the specified time and nearest location.
    """
    ds = cm1.input.era5.aws(time, **kwargs)
    # map negative longitude to 0-360 degreeE
    lon = lon % (360 * units.degreeE)
    ds = ds.sel(
        longitude=lon,
        latitude=lat,
        method="nearest",
        tolerance=5 * units.deg,
    )

    return ds


def era5_model_level(time: pd.Timestamp, lat: Quantity, lon: Quantity, **kwargs):
    """
    Retrieve ERA5 dataset for a specific time and location.

    Parameters:
        time (pd.Timestamp): The time for the dataset retrieval.
        lat (float): Latitude of the desired location.
        lon (float): Longitude of the desired location.
        **kwargs: Additional arguments to pass to the data retrieval function.

    Returns:
        xarray.Dataset: ERA5 dataset for the specified time and nearest location.
    """
    ds = cm1.input.era5.model_level(time, **kwargs)
    # map negative longitude to 0-360 degreeE
    lon = lon % (360 * units.degreeE)
    ds = ds.sel(
        longitude=lon,
        latitude=lat,
        method="nearest",
        tolerance=5 * units.deg,
    )

    return ds


def era5_pressure_level(time: pd.Timestamp, lat: Quantity, lon: Quantity, **kwargs):
    """
    Retrieve ERA5 dataset for a specific time and location.

    Parameters:
        time (pd.Timestamp): The time for the dataset retrieval.
        lat (float): Latitude of the desired location.
        lon (float): Longitude of the desired location.
        **kwargs: Additional arguments to pass to the data retrieval function.

    Returns:
        xarray.Dataset: ERA5 dataset for the specified time and nearest location.
    """
    ds = cm1.input.era5.pressure_level(time, **kwargs)
    # map negative longitude to 0-360 degreeE
    lon = lon % (360 * units.degreeE)
    ds = ds.sel(
        longitude=lon,
        latitude=lat,
        method="nearest",
        tolerance=5 * units.deg,
    )

    return ds


def get_ofile(args: argparse.Namespace) -> Path:
    """
    Generate a temporary file path for caching the dataset.

    Parameters:
        args (Namespace): Parsed command-line arguments.

    Returns:
        pathlib.Path: File path for the temporary file.
    """
    ofile = (
        TMPDIR
        / f"{pd.to_datetime(args.time).strftime('%Y%m%d_%H%M%S')}.{args.lat:~}.{args.lon:~}.nc"
    )
    return ofile


def get_case(case: str) -> xarray.Dataset:
    """
    Retrieve a predefined sounding case dataset.

    Parameters:
        case (str): Name of the sounding case.

    Returns:
        xarray.Dataset: Dataset corresponding to the specified case.
    """
    file_path = os.path.join(soundings_path, f"input_sounding_{case}")
    ds = read_from_txt(file_path)
    ds.attrs.update({"case": case})
    return ds


# Functions for specific sounding cases
def trier():
    return get_case("trier")


def jordan_allmean():
    return get_case("jordan_allmean")


def jordan_hurricane():
    return get_case("jordan_hurricane")


def rotunno_emanuel():
    return get_case("rotunno_emanuel")


def dunion_MT():
    return get_case("dunion_MT")


def bryan_morrison():
    return get_case("bryan_morrison")


def seabreeze_test():
    return get_case("seabreeze_test")


def read_from_txt(file_path: typing.Union[str, Path]) -> xarray.Dataset:
    """
    Read a CM1 sounding file format and convert it into an xarray Dataset.

    Parameters:
        file_path (str): Path to the sounding file.

    Returns:
        xarray.Dataset: Dataset containing sounding data with appropriate units.
    """
    # Open the file and read the first line (header with surface variables)
    with open(file_path, "r") as file:
        # Read the first line
        header = file.readline().strip()
        # Split the header into surface variables
        surface_pressure, surface_theta, surface_mixing_ratio = map(
            float, header.split()
        )

    # Read the remaining columns into a pandas DataFrame
    # Specify column names and skip the first line (header)
    column_names = ["Z", "theta", "qv", "U", "V"]
    df = pd.read_csv(file_path, sep=r"\s+", skiprows=1, names=column_names)
    # Rename the index to "level"
    df = df.rename_axis("level")
    ds = df.to_xarray()

    # Add units and calculate pressure
    ds["qv"] = ds["qv"] * units.g / units.kg  # mixing ratio
    ds["qv"].attrs["long_name"] = "water vapor mixing ratio"
    ds["Q"] = mcalc.specific_humidity_from_mixing_ratio(ds["qv"])
    ds["SP"] = surface_pressure
    ds["SP"] *= units.hPa
    ds["theta"] = ds["theta"] * units.K
    ds["Z"] = ds["Z"] * units.m
    ds["Z"].attrs["long_name"] = "geopotential height"
    p_bot = ds.SP.copy()  # copy to avoid modifying the original when adding dp
    z_bot = 0.0 * units.m
    P = []  # pressure array

    for level in ds.level:
        # Add effect of water vapor. Water vapor lowers the density of air given the same temperature.
        # Find the virtual temperature, the temperature of dry air with the same density as moist air.
        T = mcalc.temperature_from_potential_temperature(
            p_bot, ds.theta.sel(level=level)
        )
        Tv = mcalc.thermo.virtual_temperature(T, ds.qv.sel(level=level))
        dz = ds.Z.sel(level=level) - z_bot
        p_bot = p_bot * np.exp(-metpy.constants.g * dz / metpy.constants.Rd / Tv)
        assert (
            p_bot >= 0 * units.hPa
        ), f"{file_path} p_bot<0 {p_bot:~} dz {dz:~} Tv {Tv:~}"

        P.append(p_bot.item().m_as(ds.SP.metpy.units))
        z_bot = ds.Z.sel(level=level)

    ds["P"] = xarray.DataArray(P, coords=[ds.level])
    ds["P"] *= ds.SP.metpy.units
    ds["T"] = mcalc.temperature_from_potential_temperature(ds["P"], ds["theta"])
    ds["Tv"] = mcalc.thermo.virtual_temperature(ds.T, ds.qv)

    # Add surface data
    ds["surface_potential_temperature"] = surface_theta
    ds["surface_potential_temperature"] *= units.K
    ds["surface_mixing_ratio"] = surface_mixing_ratio
    ds["surface_mixing_ratio"] *= units.g / units.kg
    ds["surface_geopotential_height"] = 0.0
    ds["surface_geopotential_height"] *= units.m
    ds["U"] = ds["U"] * units.m / units.s
    ds["V"] = ds["V"] * units.m / units.s
    return ds


def main() -> None:
    """
    Main function for loading ERA5 data and printing sounding data.

    If cached data exists, it is loaded from a local file. Otherwise,
    it retrieves data based on user input and caches it for later use.
    """

    import pickle

    args = parse_args()

    valid_time = pd.to_datetime(args.time)
    ofile = get_ofile(args)
    if os.path.exists(ofile):
        logging.warning(f"read {ofile}")
        with open(ofile, "rb") as file:
            ds = pickle.load(file)
    else:
        ds = cm1.input.era5.model_level(
            valid_time,
            glade=args.glade,
        )
        with open(ofile, "wb") as file:
            logging.warning(f"pickle dump {ofile}")
            pickle.dump(ds, file)

    ds = era5_model_level(valid_time, args.lat, args.lon)

    print(to_txt(ds))


def to_txt(ds: xarray.Dataset) -> str:
    """
    Convert an xarray Dataset into a formatted string suitable for CM1 or WRF input.

    Parameters:
        ds (xarray.Dataset): Dataset containing atmospheric profiles.

    Returns:
        str: Formatted string representing the sounding data.

    The CM1 input sounding file format is as follows:
    One-line header containing:   sfc pres (mb)    sfc theta (K)    sfc qv (g/kg)

      (Note1: here, "sfc" refers to near-surface atmospheric conditions.
       Technically, this should be z = 0, but in practice is obtained from the
       standard reporting height of 2 m AGL/ASL from observations)
      (Note2: land-surface temperature and/or sea-surface temperature (SST) are
       specified elsewhere: see tsk0 in namelist.input and/or tsk array in
       init_surface.F)

    Then, the following lines are:   z (m)    theta (K)   qv (g/kg)    u (m/s)    v (m/s)

      (Note3: # of levels is arbitrary)

        Index:   sfc    =  surface (technically z=0, but typically from 2 m AGL/ASL obs)
                 z      =  height AGL/ASL
                 pres   =  pressure
                 theta  =  potential temperature
                 qv     =  mixing ratio
                 u      =  west-east component of velocity
                 v      =  south-north component of velocity

    Note4:  For final line of input_sounding file, z (m) must be greater than the model top
            (which is nz * dz when stretch_z=0, or ztop when stretch_z=1,  etc)
    """
    # Find the level label with the maximum pressure.
    # Don't assign level with highest numeric value. This worked for ERA5 but not CM1 input soundings.
    sfc = ds.level.sel(level=ds.P.compute().idxmax())
    # ds.SP (surface pressure) is half-level below ds.level.max()
    sfc_pres = ds.SP if "SP" in ds else ds.P.sel(level=sfc)
    ds["theta"] = mcalc.potential_temperature(ds.P, ds.T).metpy.convert_units("K")
    sfc_theta_K = (
        ds["surface_potential_temperature"].metpy.convert_units("K")
        if "surface_potential_temperature" in ds
        else ds["theta"].sel(level=sfc)
    )
    # Convert specific humidity Q to mixing ratio qv.
    if "qv" in ds:
        logging.warning(
            "ignoring qv already in Dataset. Recompute from specific humidity Q"
        )
    ds["qv"] = mcalc.mixing_ratio_from_specific_humidity(ds["Q"]).metpy.convert_units(
        "g/kg"
    )
    sfc_qv_gkg = (
        ds["surface_mixing_ratio"].metpy.convert_units("g/kg")
        if "surface_mixing_ratio" in ds
        else ds.qv.sel(level=sfc)
    )

    s = f"{sfc_pres.compute().item().m_as('hPa')} {sfc_theta_K.values} {sfc_qv_gkg.values}\n"
    s += (
        ds[["Z", "theta", "qv", "U", "V"]]
        .drop_vars(["latitude", "longitude", "time"])
        .to_dataframe()
        .sort_values("Z")  # from surface upward
        .to_csv(sep=" ", header=False, index=False)
    )

    return s


if __name__ == "__main__":
    main()
