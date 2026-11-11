import logging
import os
from pathlib import Path
from typing import Optional, Tuple

import metpy.calc as mcalc
import numpy as np
import pandas as pd
import s3fs
import xarray as xr
from metpy.constants import Rd, g
from metpy.units import units
from pint import Quantity

from cm1.utils import TMPDIR

# ---------------------------------------------------------------------------
# Configurable paths
# ---------------------------------------------------------------------------

# Override by setting the environment variable CM1_ERA5_INVARIANT_DIR, e.g.:
#   export CM1_ERA5_INVARIANT_DIR=/my/local/copy/of/invariant
_DEFAULT_INVARIANT_DIR = "/glade/work/ahijevyc/share/ds633.6-ERA5ML/e5.oper.invariant"
ERA5_INVARIANT_DIR: str = os.environ.get(
    "CM1_ERA5_INVARIANT_DIR", _DEFAULT_INVARIANT_DIR
)


# ---------------------------------------------------------------------------
# Campaign storage loader
# ---------------------------------------------------------------------------


def load_from_campaign(
    time: pd.Timestamp,
    rdaindex: str,
    level_type: str,
    varnames: list,
    start_end_str: str,
    drop_variables: list = None,
) -> xr.Dataset:
    """
    Load an ERA5 dataset for *time* from GLADE campaign storage.

    Parameters
    ----------
    time : pd.Timestamp
        Desired timestamp.
    rdaindex : str
        RDA dataset index (e.g. ``"d633006"``).
    level_type : str
        ERA5 level type string (e.g. ``"e5.oper.an.pl"``).
    varnames : list[str]
        Variable file-name suffixes to open.
    start_end_str : str
        ``YYYYMMDDHHH_YYYYMMDDHHH`` portion of each filename.
    drop_variables : list[str], optional
        Variables to drop on open (default: ``["utc_date"]``).

    Returns
    -------
    xr.Dataset
    """
    if drop_variables is None:
        drop_variables = ["utc_date"]

    base = Path("/glade/campaign/collections/rda/data") / rdaindex / level_type
    year_month_dir = time.strftime("%Y%m")
    local_files = [
        base / year_month_dir / f"{level_type}.{v}.{start_end_str}.nc" for v in varnames
    ]

    ds = xr.open_mfdataset(
        local_files, compat="override", drop_variables=drop_variables
    )
    ds.attrs["source_files"] = [str(p) for p in local_files]
    logging.info(
        f"Opened {len(local_files)} local {level_type} files; selecting {time}"
    )
    return ds.sel(time=time)


# ---------------------------------------------------------------------------
# Geopotential helpers
# ---------------------------------------------------------------------------


def compute_z_level(
    ds: xr.Dataset, lev: int, z_h: xr.DataArray
) -> Tuple[xr.DataArray, xr.DataArray]:
    r"""
    Compute the geopotential at full level *lev* and the overlying half-level.

    Parameters
    ----------
    ds : xr.Dataset
        Must contain ``"Tv"`` (virtual temperature) and ``"P_half"``.
    lev : int
        Full-level index.
    z_h : xr.DataArray
        Geopotential at the lower half-level.

    Returns
    -------
    z_h : xr.DataArray
        Updated geopotential at the overlying half-level.
    z_f : xr.DataArray
        Geopotential at the full level.

    References
    ----------
    ERA5: Compute pressure and geopotential on model levels, geopotential height,
    and geometric height.  ECMWF Confluence.
    https://confluence.ecmwf.int/display/CKB/ERA5%3A+compute+pressure+and+geopotential+on+model+levels%2C+geopotential+height+and+geometric+height
    """
    t_level = ds["Tv"].sel(level=lev)
    ph_lev = ds["P_half"].sel(half_level=lev)
    ph_levplusone = ds["P_half"].sel(half_level=lev + 1)
    pf_lev = ds["P"].sel(level=lev)

    if lev == 1:
        dlog_p = np.log(ph_levplusone / (0.1 * units.Pa))
        alpha = np.log(2)
    else:
        dlog_p = np.log(ph_levplusone / ph_lev)
        alpha = np.log(ph_levplusone / pf_lev)

    z_f = z_h + (t_level * Rd * alpha)
    z_f = z_f.drop_vars("half_level").assign_coords(level=lev)

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
    Quantify an ERA5 invariant dataset.

    Squeezes any residual time dimension, renames the ``Z`` geopotential,
    replaces unit strings that MetPy cannot handle, and derives
    ``surface_geopotential_height``.
    """
    invariant = invariant.squeeze()
    invariant = invariant.rename_vars({"Z": "surface_geopotential"})
    for var in invariant:
        u_str = invariant[var].attrs.get("units", "")
        if u_str in {"(0-1)", "-", "index"}:
            logging.debug(
                f"Replacing unquantifiable unit '{u_str}' on {var!r} with '1'"
            )
            invariant[var].attrs["units"] = "1"
    invariant = invariant.metpy.quantify()
    invariant["surface_geopotential_height"] = invariant["surface_geopotential"] / g
    return invariant


# ---------------------------------------------------------------------------
# S3 singleton  (module-level, resettable)
# ---------------------------------------------------------------------------

_s3_instance: Optional[s3fs.S3FileSystem] = None


def get_s3(reset: bool = False) -> s3fs.S3FileSystem:
    """
    Return a (cached) anonymous S3FileSystem.

    Parameters
    ----------
    reset : bool
        If *True*, discard any cached instance and create a fresh one.
        Useful when a previous connection has gone stale.
    """
    global _s3_instance
    if reset or _s3_instance is None:
        _s3_instance = s3fs.S3FileSystem(anon=True)
    return _s3_instance


# ---------------------------------------------------------------------------
# Model-level ERA5 (campaign storage)
# ---------------------------------------------------------------------------


def model_level(time: pd.Timestamp | np.datetime64 | str) -> xr.Dataset:
    """
    Load native model-level ERA5 data from GLADE campaign storage.

    The invariant surface geopotential directory is controlled by the
    ``CM1_ERA5_INVARIANT_DIR`` environment variable (default:
    ``/glade/work/ahijevyc/share/ds633.6-ERA5ML/e5.oper.invariant``).

    Parameters
    ----------
    time : timestamp-like
        Desired valid time.

    Returns
    -------
    xr.Dataset
        Model-level dataset on the Gaussian lat-lon grid with derived
        ``P``, ``P_half``, ``Z``, ``Z_half``, ``Tv``, and
        ``surface_geopotential_height``.
    """
    time = pd.Timestamp(time)
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
        drop_variables=["zero", "utc_date"],
    )
    ds = ds.metpy.quantify()

    ds["P"] = (ds.a_model + ds.b_model * ds.SP).transpose(*ds.U.dims)
    ds["P"].attrs["long_name"] = "pressure"
    ds["P_half"] = ds.a_half + ds.b_half * ds.SP
    ds["P_half"].attrs["long_name"] = "pressure"

    # Load invariant fields from the configurable directory
    invariant_path = Path(ERA5_INVARIANT_DIR)
    if not invariant_path.exists():
        raise FileNotFoundError(
            f"ERA5 invariant directory not found: {invariant_path}\n"
            "Set the CM1_ERA5_INVARIANT_DIR environment variable to override."
        )
    invariant_files = list(invariant_path.glob("*.nc"))
    logging.debug(f"Invariant files: {invariant_files}")

    invariant = xr.open_mfdataset(invariant_files, drop_variables=["utc_date", "time"])
    invariant.attrs["source_files"] = [str(p) for p in invariant_files]
    assert invariant.latitude.size == 640, (
        "Expected invariant fields on the N320 Gaussian grid (640 latitudes)"
    )
    invariant = quantify_invariant(invariant)
    ds = ds.merge(invariant, combine_attrs="drop_conflicts")

    logging.info("Computing geopotential height via hypsometric equation")
    ds["Tv"] = mcalc.thermo.virtual_temperature(ds.T, ds.Q)
    z_h = ds.surface_geopotential.assign_coords(half_level=ds.half_level.max())
    Z, Z_h = [], [z_h]
    for level in ds.level[::-1]:
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

    return ds


# ---------------------------------------------------------------------------
# Pressure-level ERA5 (campaign storage)
# ---------------------------------------------------------------------------


def pressure_level(time: pd.Timestamp) -> xr.Dataset:
    """
    Load pressure-level ERA5 data from GLADE campaign storage.

    Parameters
    ----------
    time : pd.Timestamp
        Desired valid time.

    Returns
    -------
    xr.Dataset
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
    ds_pl = ds_pl.metpy.quantify()
    ds_pl["P"] = ds_pl.level
    ds_pl["Z"] /= g

    lastdayofmonth = time + pd.offsets.MonthEnd(0)
    sfc_start_end = (
        f"{time.strftime('%Y%m')}0100_"
        f"{time.strftime('%Y%m')}{lastdayofmonth.strftime('%d')}23"
    )

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
        sfc_start_end,
    )
    ds_sfc = ds_sfc.metpy.quantify()
    ds_sfc["surface_potential_temperature"] = mcalc.potential_temperature(
        ds_sfc.SP, ds_sfc.VAR_2T
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

    ds = ds_pl.merge(ds_sfc, compat="no_conflicts").merge(
        invariant, compat="no_conflicts"
    )

    # Consolidate source file lists from all components
    ds.attrs["history_sources"] = [
        f
        for component in [ds_pl, ds_sfc, invariant]
        for f in component.attrs.get("source_files", [])
    ]
    return ds


# ---------------------------------------------------------------------------
# AWS / S3 ERA5
# ---------------------------------------------------------------------------


def aws(time: pd.Timestamp) -> xr.Dataset:
    """
    Retrieve ERA5 pressure-level and surface data from the NSF NCAR S3 bucket.

    Files are cached locally under ``$TMPDIR/era5_cache/`` to avoid repeated
    downloads.

    Parameters
    ----------
    time : pd.Timestamp
        Desired valid time.

    Returns
    -------
    xr.Dataset
    """
    S3_BUCKET = "nsf-ncar-era5"
    CACHE_DIR = TMPDIR / "era5_cache"
    CACHE_DIR.mkdir(parents=True, exist_ok=True)

    def download_from_s3(var: str, level_type: str, start_end_str: str) -> Path:
        """Download *var* from S3 if not already cached; return local path."""
        file_name = f"e5.oper.{level_type}.{var}.{start_end_str}.nc"
        cache_path = CACHE_DIR / file_name
        if cache_path.exists():
            logging.info(f"Cache hit: {var} {time}")
            return cache_path

        s3 = get_s3()
        year_month = start_end_str[:6]
        s3_path = f"{S3_BUCKET}/e5.oper.{level_type}/{year_month}/{file_name}"
        if not s3.exists(s3_path):
            raise FileNotFoundError(f"Not found in S3: {s3_path}")

        logging.warning(f"Downloading {s3_path} …")
        with s3.open(s3_path, "rb") as f:
            xr.open_dataset(f).to_netcdf(cache_path)
        logging.warning(f"Cached to {cache_path}")
        return cache_path

    daily_str = f"{time.strftime('%Y%m%d00')}_{time.strftime('%Y%m%d23')}"
    pl_paths = [
        download_from_s3(v, "an.pl", daily_str)
        for v in [
            "128_133_q.ll025sc",
            "128_130_t.ll025sc",
            "128_131_u.ll025uv",
            "128_132_v.ll025uv",
            "128_135_w.ll025sc",
            "128_129_z.ll025sc",
        ]
    ]
    ds_pl = (
        xr.open_mfdataset(pl_paths, compat="override")
        .drop_vars("utc_date")
        .sel(time=time)
    )
    ds_pl = ds_pl.metpy.quantify()
    ds_pl["P"] = ds_pl.level * ds_pl.level.metpy.units
    ds_pl["Z"] /= g

    lastdayofmonth = time + pd.offsets.MonthEnd(0)
    monthly_str = (
        f"{time.strftime('%Y%m')}0100_"
        f"{time.strftime('%Y%m')}{lastdayofmonth.strftime('%d')}23"
    )
    sfc_paths = [
        download_from_s3(v, "an.sfc", monthly_str)
        for v in [
            "128_134_sp.ll025sc",
            "128_165_10u.ll025sc",
            "128_166_10v.ll025sc",
            "128_167_2t.ll025sc",
            "128_168_2d.ll025sc",
        ]
    ]
    ds_sfc = (
        xr.open_mfdataset(sfc_paths, compat="override")
        .drop_vars("utc_date")
        .sel(time=time)
    )
    ds_sfc = ds_sfc.metpy.quantify()
    ds_sfc["surface_potential_temperature"] = mcalc.potential_temperature(
        ds_sfc.SP, ds_sfc.VAR_2T
    )
    ds_sfc["surface_mixing_ratio"] = mcalc.mixing_ratio_from_specific_humidity(
        mcalc.specific_humidity_from_dewpoint(ds_sfc.SP, ds_sfc.VAR_2D)
    )

    invariant_paths = [
        download_from_s3(v, "invariant", "1979010100_1979010100")
        for v in INVARIANT_VARNAMES
    ]
    invariant = xr.open_mfdataset(invariant_paths, drop_variables=["utc_date", "time"])
    invariant = quantify_invariant(invariant)

    return ds_pl.merge(ds_sfc, compat="no_conflicts").merge(
        invariant, compat="no_conflicts"
    )


# ---------------------------------------------------------------------------
# Spatial selection helper
# ---------------------------------------------------------------------------


def nearest_grid_block_sel(
    dataset: xr.Dataset,
    lat: Quantity,
    lon: Quantity,
    n: int = 3,
    **kwargs,
) -> dict:
    """
    Build a ``.sel()`` dict that selects the nearest *n* × *n* grid block.

    Parameters
    ----------
    dataset : xr.Dataset
        ERA5 dataset with ``latitude`` and ``longitude`` coordinates.
    lat, lon : Quantity
        Target location (units: degrees).
    n : int
        Side length of the square block; must be a positive odd integer.

    Returns
    -------
    dict
        ``{"latitude": ..., "longitude": ...}`` suitable for ``ds.sel(**result)``.
    """
    if n < 1 or n % 2 != 1:
        raise ValueError("n must be a positive odd integer.")
    if "latitude" not in dataset.coords or "longitude" not in dataset.coords:
        raise ValueError("Dataset must contain 'latitude' and 'longitude' coordinates.")

    lon = lon % (360 * units.degree)
    center = dataset.sel(latitude=lat, longitude=lon, method="nearest", **kwargs)
    center_lat = center.latitude.values
    center_lon = center.longitude.values

    lat_vals = dataset.latitude.values
    lon_vals = dataset.longitude.values
    lat_idx = int(np.abs(lat_vals - center_lat).argmin())
    lon_idx = int(np.abs(lon_vals - center_lon).argmin())

    half = n // 2
    lat_start = max(lat_idx - half, 0)
    lat_end = min(lat_idx + half, len(lat_vals) - 1)
    lat_sel = lat_vals[lat_start : lat_end + 1]

    lon_size = len(lon_vals)
    lon_sel = np.array(
        [lon_vals[(lon_idx + i) % lon_size] for i in range(-half, half + 1)]
    )

    return {"latitude": lat_sel, "longitude": lon_sel}
