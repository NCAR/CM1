import logging
import os
from functools import lru_cache
from pathlib import Path
from typing import Tuple

import metpy.calc as mcalc
import numpy as np
import pandas as pd
import s3fs
import xarray as xr
from metpy.constants import Rd, g
from metpy.units import units

from cm1.utils import TMPDIR


def load_from_campaign(
    time: pd.Timestamp,
    rdaindex: str,
    level_type: str,
    varnames: list,
    start_end_str: str,
    drop_variables: list = ["utc_date"],
) -> xr.Dataset:
    """
    Load ERA5 dataset for specified time from
    campaign storage

    Parameters
    ----------
    time : pd.Timestamp
        Desired timestamp for data retrieval.
    level_type : str
        Type of ERA5 data (e.g., 'e5.oper.an.pl', 'e5.oper.an.sfc').
    varnames : list
        List of variable names to load.
    start_end_str : str
        start and end time part of filename
    drop_variables : list, optional
        List of variables to drop from the dataset (default is ["utc_date"]).

    Returns
    -------
    xr.Dataset
        Dataset containing ERA5 data for the specified time and configuration.
    """

    # model level invariant files don't have no year_month_dir.
    year_month_dir = (
        ""
        if rdaindex == "d633006" and level_type == "e5.oper.invariant"
        else time.strftime("%Y%m")
    )
    local_files = [
        Path("/glade/campaign/collections/rda/data")
        / rdaindex
        / level_type
        / year_month_dir
        / f"{level_type}.{varname}.{start_end_str}.nc"
        for varname in varnames
    ]

    ds = xr.open_mfdataset(local_files, drop_variables=drop_variables)
    logging.info(f"opened {len(local_files)} local {level_type} files")
    logging.debug(local_files)
    logging.info(f"selected {time}")
    ds = ds.sel(time=time)

    return ds


def compute_z_level(
    ds: xr.Dataset, lev: int, z_h: xr.DataArray
) -> Tuple[xr.DataArray, xr.DataArray]:
    r"""
    Compute the geopotential at a full level and the overlying half-level.

    This function calculates the geopotential \( z_f \) at a specified full level
    (given by `lev`) and updates the geopotential at the overlying half-level
    \( z_h \). The calculation is based on the virtual temperature at full level
    and pressure at half-levels.

    Parameters:
    ----------
    ds : xr.Dataset
        The dataset containing the variables required for computation,
        specifically "Tv" (virtual temperature) and "P_half" (pressure on half-levels).
    lev : int
        The level index for the desired full-level geopotential calculation.
    z_h : xr.DataArray
        The initial geopotential height at the lower half-level.

    Returns:
    -------
    Tuple[xr.DataArray, xr.DataArray]
        A tuple containing:
            - `z_h`: The updated geopotential at the overlying half-level.
            - `z_f`: The computed geopotential at the specified full level.

    References
    ----------
    ERA5: Compute pressure and geopotential on model levels, geopotential height, and geometric height.
    ECMWF Confluence: https://confluence.ecmwf.int/display/CKB/ERA5%3A+compute+pressure+and+geopotential+on+model+levels%2C+geopotential+height+and+geometric+height
    https://opensky.ucar.edu/islandora/object/%3A3444
    NCAR/TN-396+STR
    Dec 1993
    Vertical Interpolation and Truncation of Model-Coordinate Data
    Kevin E. Trenberth
    Jeffery C. Berry
    Lawrence E. Buja
    """
    # Virtual temperature at the specified level
    t_level = ds["Tv"].sel(level=lev)

    # Pressures at the half-levels above and below
    ph_lev = ds["P_half"].sel(half_level=lev)
    ph_levplusone = ds["P_half"].sel(half_level=lev + 1)
    pf_lev = ds["P"].sel(level=lev)

    if lev == 1:
        dlog_p = np.log(ph_levplusone / (0.1 * units.Pa))
        alpha = np.log(2)
    else:
        dlog_p = np.log(ph_levplusone / ph_lev)
        # TODO: understand IFS formulation for alpha. See IFS-documentation-cy47r3 eqn 2.23
        # alphaIFS = 1.0 - ((ph_lev / (ph_levplusone - ph_lev)) * dlog_p)
        # @ahijevyc formulation
        alpha = np.log(ph_levplusone / pf_lev)
        # Make sure official and @ahijevyc values are close to each other.
        # assert np.allclose(
        #    alphaIFS.load(), alpha.load(), atol=1e-2
        # ), f"{np.abs((alphaIFS-alpha)).max()}"

    # Calculate the full-level geopotential `z_f`
    # Integrate from previous (lower) half-level `z_h` to the
    z_f = z_h + (t_level * Rd * alpha)
    z_f = z_f.drop_vars("half_level")
    z_f = z_f.assign_coords(level=lev)

    # Update the half-level geopotential `z_h`
    z_h = z_h + (t_level * Rd * dlog_p)
    z_h = z_h.assign_coords(half_level=lev).drop_vars("level")

    return z_h, z_f


INVARIANT_VARNAMES = [
    "128_026_cl.ll025sc",
    "128_027_cvl.ll025sc",
    "128_028_cvh.ll025sc",
    "128_029_tvl.ll025sc",
    "128_030_tvh.ll025sc",
    "128_043_slt.ll025sc",
    "128_074_sdfor.ll025sc",
    "128_129_z.ll025sc",
    "128_160_sdor.ll025sc",
    "128_161_isor.ll025sc",
    "128_162_anor.ll025sc",
    "128_163_slor.ll025sc",
    "128_172_lsm.ll025sc",
    "228_007_dl.ll025sc",
]


def quantify_invariant(invariant: xr.Dataset) -> xr.Dataset:
    """
    Quantify invariant Dataset
    Squeeze time dimension if present
    rename Z surface geopotential
    remove units that metpy doesn't handle
    convert surface geopotential to height
    """
    invariant = invariant.squeeze()  # squeeze time, if present
    invariant = invariant.rename_vars({"Z": "surface_geopotential"})
    for var in invariant:
        u = invariant[var].attrs["units"]
        if u in ["(0-1)", "-", "index"]:
            logging.info(f"can't quantify {var} unit string '{u}'")
            invariant[var].attrs["units"] = "1"
    invariant = invariant.metpy.quantify()
    invariant["surface_geopotential_height"] = invariant["surface_geopotential"] / g

    return invariant


def model_level(
    time: pd.Timestamp,
) -> xr.Dataset:
    """
    Load native model levels ERA5 dataset for specified time
    from campaign storage. On Gaussian lat-lon grid.
    SP:grid_specification = "0.28125° x ~0.28125° from 0E to 359.71875E and 89.78488N to 89.78488S (1280 x 640 Longitude/Gaussian Latitude), ~31km at Equator" ;
        SP:original_data_representation = "N320 reduced Gaussian grid, ECMWF Meteorological Archival and Retrieval System (MARS)" ;

    Parameters
    ----------
    time : pd.Timestamp
        Desired timestamp for data retrieval.

    Returns
    -------
    xr.Dataset
        Dataset containing ERA5 data for the specified time.
    """
    # get from campaign storage
    rdaindex = "d633006"

    start_hour = time.floor("6h")
    end_hour = start_hour + pd.Timedelta(5, unit="hour")
    start_end_str = f"{start_hour.strftime('%Y%m%d%H')}_{end_hour.strftime('%Y%m%d%H')}"

    ds = load_from_campaign(
        time,
        rdaindex,
        "e5.oper.an.ml",
        [
            "0_5_0_0_0_t.regn320sc",
            "0_5_0_1_0_q.regn320sc",
            "0_5_0_2_2_u.regn320uv",
            "0_5_0_2_3_v.regn320uv",
            "0_5_0_2_8_w.regn320sc",
            "128_134_sp.regn320sc",
        ],
        start_end_str,
        ["zero", "utc_date"],
    )
    ds = ds.metpy.quantify()

    # Derive pressure from a and b coefficients
    ds["P"] = ds.a_model + ds.b_model * ds.SP
    ds["P"].attrs.update(dict(long_name="pressure"))
    ds["P"] = ds["P"].transpose(*ds.U.dims)  # keep dim order consistent with U
    ds["P_half"] = ds.a_half + ds.b_half * ds.SP
    ds["P_half"].attrs.update(dict(long_name="pressure"))

    # rdahelp says /gpfs/csfs1/collections/rda/decsdata/ds630.0/P/e5.oper.invariant/201601/
    # has same resolution as the invariant surface geopotential you refer to in d633006.

    invariant_path = Path(
        # "/glade/campaign/collections/rda/decsdata/COLD_STORAGE/d630000/P/e5.oper.invariant/201601"
        f"/glade/campaign/collections/rda/data/{rdaindex}/e5.oper.invariant"
    )
    invariant = xr.open_mfdataset(
        list(invariant_path.glob("*.nc")),
        drop_variables=["utc_date", "time"],
    )
    assert (
        invariant.latitude.size == 640
    ), "expected invariant fields on Gaussian grid like ds"
    invariant = quantify_invariant(invariant)
    ds = ds.merge(invariant)

    logging.warning("filling height using hypsometric equation")
    ds["Tv"] = mcalc.thermo.virtual_temperature(ds.T, ds.Q)
    z_h = ds.surface_geopotential.assign_coords(half_level=ds.half_level.max())
    Z = []  # geopotential on full levels
    Z_h = [z_h]  # geopotential on half levels
    # Loop from last to first level (sfc upward)
    for level in ds.level[::-1]:  # accumulate geopotential in z_h upward from sfc
        z_h, z_f = compute_z_level(ds, level, z_h)
        Z.append(z_f)
        Z_h.append(z_h)
    ds["Z_half"] = xr.concat(Z_h, dim="half_level") / g
    ds["Z_half"].attrs["long_name"] = "geopotential height"
    ds["Z"] = xr.concat(Z, dim="level") / g
    ds["Z"].attrs["long_name"] = "geopotential height"
    ds["surface_geopotential_height"] = ds["surface_geopotential"] / g
    ds["surface_geopotential_height"].attrs["long_name"] = (
        "geopotential height at surface"
    )
    # ds = ds.drop_dims("half_level") # why drop this?

    return ds


def pressure_level(
    time: pd.Timestamp,
) -> xr.Dataset:
    """
    Load ERA5 dataset for specified time and configuration.
    Comes from campaign storage.

    Parameters
    ----------
    time : pd.Timestamp
        Desired timestamp for data retrieval.

    Returns
    -------
    xr.Dataset
        Dataset containing ERA5 data for the specified time and configuration.
    """
    rdaindex = "d633000"
    start = time.floor("1d")
    end = start + pd.Timedelta(23, unit="hour")
    start_end_str = f"{start.strftime('%Y%m%d%H')}_{end.strftime('%Y%m%d%H')}"

    ds_pl = load_from_campaign(
        time,
        rdaindex,
        "e5.oper.an.pl",
        [
            "128_129_z.ll025sc",
            "128_130_t.ll025sc",
            "128_131_u.ll025uv",
            "128_132_v.ll025uv",
            "128_133_q.ll025sc",
            "128_135_w.ll025sc",
        ],
        start_end_str,
    )

    # Quantify the entire dataset to ensure all variables are pint.Quantities
    # for consistent data types throughout the merging process.
    ds_pl = ds_pl.metpy.quantify()

    # Create pressure variable from the 'level' coordinate.
    # After quantification, 'level' is already a pint.Quantity.
    ds_pl["P"] = ds_pl.level
    ds_pl["Z"] /= g

    lastdayofmonth = time + pd.offsets.MonthEnd(0)
    start_end_str = f"{time.strftime('%Y%m')}0100_{time.strftime('%Y%m')}{lastdayofmonth.strftime('%d')}23"

    ds_sfc = load_from_campaign(
        time,
        rdaindex,
        "e5.oper.an.sfc",
        [
            "128_134_sp.ll025sc",
            "128_165_10u.ll025sc",
            "128_166_10v.ll025sc",
            "128_235_skt.ll025sc",
            "128_167_2t.ll025sc",
            "128_168_2d.ll025sc",
        ],
        start_end_str,
    )

    ds_sfc = ds_sfc.metpy.quantify()
    ds_sfc["surface_potential_temperature"] = mcalc.potential_temperature(
        ds_sfc.SP,
        ds_sfc.VAR_2T,
    )
    ds_sfc["surface_mixing_ratio"] = mcalc.mixing_ratio_from_specific_humidity(
        mcalc.specific_humidity_from_dewpoint(ds_sfc.SP, ds_sfc.VAR_2D)
    )

    invariant = load_from_campaign(
        pd.to_datetime("19790101"),
        rdaindex,
        "e5.oper.invariant",
        INVARIANT_VARNAMES,
        "1979010100_1979010100",
    ).drop_vars("time")
    invariant = quantify_invariant(invariant)
    ds = ds_pl.merge(ds_sfc).merge(invariant)

    return ds


@lru_cache(maxsize=1)
def get_s3():
    return s3fs.S3FileSystem(anon=True)


def aws(time: pd.Timestamp) -> xr.Dataset:
    """
    Retrieve ERA5 data from an S3 bucket and cache it locally.

    Parameters
    ----------
    time : pd.Timestamp
        Desired timestamp for data retrieval.

    Returns
    -------
    xr.Dataset
        Dataset containing ERA5 data for the specified time, downloaded from S3.
    """
    S3_BUCKET = "nsf-ncar-era5"
    CACHE_DIR = TMPDIR / "era5_cache"
    if not os.path.exists(CACHE_DIR):
        os.makedirs(CACHE_DIR)

    def download_from_s3(var, level_type, start_end_str):
        file_name = f"e5.oper.{level_type}.{var}.{start_end_str}.nc"
        cache_file_path = CACHE_DIR / file_name
        if os.path.exists(cache_file_path):
            logging.info(f"Found cached s3 {var} {time}")
            return cache_file_path

        s3 = get_s3()

        year_month_dir = start_end_str[0:6]
        s3_file_path = S3_BUCKET + f"/e5.oper.{level_type}/{year_month_dir}/{file_name}"
        if not s3.exists(s3_file_path):
            raise FileNotFoundError(f"{s3_file_path} not found in S3 bucket")

        logging.warning(f"Downloading {s3_file_path} from S3...")
        with s3.open(s3_file_path, "rb") as f:
            xr.open_dataset(f).to_netcdf(cache_file_path)
        logging.warning(f"Downloaded and cached: {cache_file_path}")
        return cache_file_path

    start_end_str = f"{time.strftime('%Y%m%d00')}_{time.strftime('%Y%m%d23')}"
    cache_file_paths = [
        download_from_s3(var, "an.pl", start_end_str)
        for var in [
            "128_133_q.ll025sc",
            "128_130_t.ll025sc",
            "128_131_u.ll025uv",
            "128_132_v.ll025uv",
            "128_135_w.ll025sc",
            "128_129_z.ll025sc",
        ]
    ]

    ds_pl = xr.open_mfdataset(cache_file_paths).drop_vars("utc_date")
    ds_pl = ds_pl.sel(time=time)
    ds_pl = ds_pl.metpy.quantify()
    ds_pl["P"] = ds_pl.level * ds_pl.level.metpy.units
    ds_pl["Z"] /= g

    lastdayofmonth = time + pd.offsets.MonthEnd(0)
    start_end_str = f"{time.strftime('%Y%m')}0100_{time.strftime('%Y%m')}{lastdayofmonth.strftime('%d')}23"
    cache_file_paths = [
        download_from_s3(var, "an.sfc", start_end_str)
        for var in [
            "128_134_sp.ll025sc",
            "128_165_10u.ll025sc",
            "128_166_10v.ll025sc",
            "128_167_2t.ll025sc",
            "128_168_2d.ll025sc",
        ]
    ]

    ds_sfc = xr.open_mfdataset(cache_file_paths).drop_vars("utc_date")
    ds_sfc = ds_sfc.sel(time=time)
    ds_sfc = ds_sfc.metpy.quantify()

    ds_sfc["surface_potential_temperature"] = mcalc.potential_temperature(
        ds_sfc.SP,
        ds_sfc.VAR_2T,
    )
    ds_sfc["surface_mixing_ratio"] = mcalc.mixing_ratio_from_specific_humidity(
        mcalc.specific_humidity_from_dewpoint(ds_sfc.SP, ds_sfc.VAR_2D)
    )

    cache_file_paths = [
        download_from_s3(var, "invariant", "1979010100_1979010100")
        for var in INVARIANT_VARNAMES
    ]

    invariant = xr.open_mfdataset(cache_file_paths, drop_variables=["utc_date", "time"])
    invariant = quantify_invariant(invariant)

    ds = ds_pl.merge(ds_sfc).merge(invariant)

    return ds
