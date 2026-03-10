import argparse
import io
import logging
import os
from pathlib import Path
from typing import Optional, TextIO, Tuple, Union

from matplotlib.figure import Figure
import matplotlib.transforms as transforms
import metpy.calc as mpcalc
import metpy.constants
from metpy.interpolate import interpolate_1d
from metpy.plots import Hodograph, SkewT
import numpy as np
import pandas as pd
import xarray as xr
from metpy.units import units
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from pint import Quantity

import cm1.input.era5
from cm1.utils import TMPDIR, parse_args

# Assuming this script is located in a subdirectory of the repository
repo_base_path = os.path.abspath(os.path.join(os.path.dirname(__file__), "../.."))
soundings_path = os.path.join(repo_base_path, "soundings")


class Sounding(xr.Dataset):
    """
    An xarray-based representation of an atmospheric sounding,
    tailored for CM1 model input/output.
    """

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
        source_file = None

        if isinstance(data_or_path, (str, Path, os.PathLike)):
            source_file = str(data_or_path)
            ds = xr.open_dataset(data_or_path)
        else:
            ds = data_or_path

        if isinstance(ds, (xr.Dataset, xr.DataArray)):
            quantified_ds = ds.metpy.quantify()
            super().__init__(quantified_ds, *args, **kwargs)
            if source_file:
                self.attrs["source_file"] = source_file
        else:
            super().__init__(ds, *args, **kwargs)

    def __setitem__(self, key, value):
        """
        Overridden to handle assignment of scalar Pint Quantities, which
        otherwise cause AttributeErrors in xarray because they lack a 'shape' attribute.
        """
        if hasattr(value, "units") and hasattr(value, "magnitude"):
            if not hasattr(value.magnitude, "shape"):
                value = np.array(value.magnitude) * value.units

        super().__setitem__(key, value)

    def plot(self, fig=None, subplot=None, **kwargs):
        """Plots the sounding using the skewt function."""

        return skewt(self, fig=fig, subplot=subplot, **kwargs)

    def to_txt(self) -> str:
        """Converts a Sounding into a formatted string suitable for CM1."""
        self.load()  # can't use item() on lazy dask array below
        sfc = self.level.sel(level=self.P.idxmax())
        sfc_pres = self.SP if "SP" in self else self.P.sel(level=sfc)
        self["theta"] = mpcalc.potential_temperature(
            self.P, self.T
        ).metpy.convert_units("K")

        sfc_theta_K = (
            self["surface_potential_temperature"].metpy.convert_units("K")
            if "surface_potential_temperature" in self
            else self["theta"].sel(level=sfc)
        )

        if "qv" not in self and "Q" in self:
            self["qv"] = mpcalc.mixing_ratio_from_specific_humidity(
                self["Q"]
            ).metpy.convert_units("g/kg")

        sfc_qv_gkg = (
            self["surface_mixing_ratio"].metpy.convert_units("g/kg")
            if "surface_mixing_ratio" in self
            else self.qv.sel(level=sfc)
        )

        header = f"{sfc_pres.item().m_as('hPa'):.2f} {sfc_theta_K.item().m_as('K'):.2f} {sfc_qv_gkg.item().m_as('g/kg'):.2f}\n"

        columns = ["Z", "theta", "qv", "U", "V"]

        body = (
            self[columns]
            .sortby("Z")
            .to_dataframe()
            .to_csv(
                sep=" ", columns=columns, header=False, index=False, float_format="%.2f"
            )
        )
        return header + body


def _integrate_pressure(Z, theta, qv, sfc_pres, sfc_theta, sfc_mix):
    """
    Integrates the hydrostatic equation upward to calculate pressure at all levels.
    Based on the physics in George Bryan's getcape.f90.
    """
    # Constants
    g = metpy.constants.g
    Cp = metpy.constants.Cp_d
    kappa = metpy.constants.kappa  # Rd/Cp
    P0 = 1000.0 * units.hPa

    theta_v = mpcalc.virtual_temperature(theta, qv)

    # Initialize boundary conditions with surface values
    z_prev = 0.0 * units.m
    pi_prev = (sfc_pres / P0) ** kappa
    theta_v_prev = mpcalc.virtual_temperature(sfc_theta, sfc_mix)

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


def _parse_cm1_txt_stream(
    stream: TextIO, case_name_hint: str, source_file: str = None
) -> Sounding:
    """Helper function to parse a CM1 text file from a text stream."""
    # Read the header: Surface Pressure (hPa), Surface Potential Temp (K), Surface Mixing Ratio (g/kg)
    header = stream.readline().strip()
    surface_pressure, surface_theta, surface_mixing_ratio = map(float, header.split())

    # CM1 sounding format: Height(m), Theta(K), Qv(g/kg), U(m/s), V(m/s)
    column_names = ["Z", "theta", "qv", "U", "V"]
    df = pd.read_csv(stream, sep=r"\s+", names=column_names, engine="python")

    df = df.rename_axis("level")
    ds = df.to_xarray()

    # Apply units and derive initial variables
    ds["qv"] *= units("g/kg")
    ds["qv"].attrs["long_name"] = "water vapor mixing ratio"
    ds["Q"] = mpcalc.specific_humidity_from_mixing_ratio(ds["qv"])

    # Convert scalar to numpy array so it has .shape attribute (needed by xarray)
    ds["SP"] = np.array(surface_pressure) * units.hPa
    ds["theta"] *= units.K
    ds["Z"] *= units.m
    ds["Z"].attrs["long_name"] = "geopotential height"

    # Integrate upward from the surface using the Exner-Hydrostatic equation.
    P = _integrate_pressure(
        ds.Z.data,
        ds.theta.data,
        ds.qv.data,
        ds.SP.data,
        surface_theta * units.K,
        surface_mixing_ratio * (units.g / units.kg),
    )

    # Attach calculated pressure and derive dependent temperatures
    ds["P"] = ("level", P)
    ds["T"] = mpcalc.temperature_from_potential_temperature(ds["P"], ds["theta"])
    ds["Tv"] = mpcalc.virtual_temperature(ds.T, ds.qv)

    # Final metadata and unit cleanup
    ds["surface_potential_temperature"] = np.array(surface_theta) * units.K
    ds["surface_mixing_ratio"] = np.array(surface_mixing_ratio) * (units.g / units.kg)
    ds["surface_geopotential_height"] = np.array(0.0) * units.m
    ds["U"] *= units.m / units.s
    ds["V"] *= units.m / units.s

    case_name = case_name_hint.replace("input_sounding_", "")
    sounding_obj = Sounding(ds)
    sounding_obj.attrs["case"] = case_name
    if source_file:
        sounding_obj.attrs["source_file"] = source_file

    return sounding_obj


def from_txt(source: Union[str, os.PathLike]) -> Sounding:
    """
    Creates a Sounding instance by parsing a special CM1 text file or string.
    """
    if isinstance(source, str) and "\n" in source:
        string_stream = io.StringIO(source)
        return _parse_cm1_txt_stream(string_stream, case_name_hint="from_string")

    path = Path(source)
    if not path.is_file():
        raise FileNotFoundError(f"Sounding file not found at path: {path}")

    with path.open("r") as f:
        return _parse_cm1_txt_stream(
            f, case_name_hint=path.stem, source_file=str(path.absolute())
        )


def get_case(case: str) -> Sounding:
    """Retrieves a predefined sounding case dataset from the soundings directory."""
    file_path = os.path.join(soundings_path, f"input_sounding_{case}")
    return from_txt(file_path)


def era5_aws(
    time: pd.Timestamp | np.datetime64 | str, lat: Quantity, lon: Quantity, **kwargs
) -> Sounding:
    """Retrieves ERA5 dataset for a specific time and location from AWS."""
    valid_time = pd.Timestamp(time)
    ds = cm1.input.era5.aws(valid_time, **kwargs)
    lon = lon % (360 * units.degreeE)
    ds = ds.sel(longitude=lon, latitude=lat, method="nearest", tolerance=5 * units.deg)
    return Sounding(ds)


def era5_model_level(
    time: pd.Timestamp | np.datetime64 | str, lat: Quantity, lon: Quantity, **kwargs
) -> Sounding:
    """Retrieves ERA5 model-level dataset for a specific time and location."""
    valid_time = pd.Timestamp(time)
    ds = cm1.input.era5.model_level(valid_time, **kwargs)
    lon = lon % (360 * units.degreeE)
    ds = ds.sel(longitude=lon, latitude=lat, method="nearest", tolerance=5 * units.deg)
    return Sounding(ds)


def era5_pressure_level(
    time: pd.Timestamp | np.datetime64 | str, lat: Quantity, lon: Quantity, **kwargs
) -> Sounding:
    """Retrieves ERA5 pressure-level dataset for a specific time and location."""
    ds = cm1.input.era5.pressure_level(time, **kwargs)
    lon = lon % (360 * units.degreeE)
    ds = ds.sel(longitude=lon, latitude=lat, method="nearest", tolerance=5 * units.deg)
    return Sounding(ds)


def get_ofile(args: argparse.Namespace) -> Path:
    """Generates a temporary file path for caching the dataset."""
    time_obj = pd.Timestamp(args.time)
    time_str = time_obj.strftime("%Y%m%dT%H%M%S")
    lat_str = f"{args.lat.m:+06.2f}{args.lat.units}".replace(" ", "")
    lon_str = f"{args.lon.m:+07.2f}{args.lon.units}".replace(" ", "")
    ofile = TMPDIR / f"{time_str}.{lat_str}.{lon_str}.nc"
    return ofile


def main() -> None:
    """Main function for loading ERA5 data and printing sounding data."""
    args = parse_args()
    valid_time = pd.to_datetime(args.time)
    ofile = get_ofile(args)

    if os.path.exists(ofile):
        logging.warning(f"Reading from cache: {ofile}")
        ds = Sounding(ofile)
    else:
        if os.path.exists("/glade/campaign"):
            ds = era5_model_level(valid_time, args.lat, args.lon)
        else:
            logging.warning(
                "No campaign storage. Getting pressure level data from AWS."
            )
            ds = era5_aws(valid_time, args.lat, args.lon)

        logging.warning(f"Caching data to: {ofile}")
        ds.metpy.dequantify().to_netcdf(ofile)

    print(ds.to_txt())


def skewt(
    ds: xr.Dataset,
    fig: Optional[Figure] = None,
    subplot: Optional[Tuple[int, int, int]] = None,
    rotation: int = 40,
    ptop: Quantity = 100 * units.hPa,
    xlim: Optional[Tuple[float, float]] = (-40, 55),
) -> SkewT:
    """
    Generates a Skew-T diagram with temperature, dewpoint, wind barbs, CAPE, CIN,
    and other atmospheric features, using data from an ERA5 dataset. A hodograph
    is also plotted as an inset within the Skew-T diagram.

    Parameters:
    ----------
    ds : xr.Dataset
        ERA5 dataset containing temperature, dewpoint, wind components, and pressure
        data for a specific time and location. Expected variables include 'T'
        (temperature), 'Q' (specific humidity), 'U' (zonal wind), 'V' (meridional wind),
        'Z' (geopotential height), and 'P' (pressure).
    fig : Optional[Figure], default=None
        Matplotlib figure object to use for plotting. If None, a new figure is created.
    subplot : Optional[Tuple[int, int, int]], default=None
        Tuple specifying (nrows, ncols, index) for subplot positioning within the
        figure. This defines the grid layout and the subplot’s position within it.
    rotation : int, default=40
        Rotation angle for the Skew-T plot lines (degrees).
    ptop : Quantity, default=100 * units.hPa
        The upper limit of pressure (in hPa) for the Skew-T plot. Data below this
        pressure level will be excluded from the plot.

    Returns:
    -------
    skew : skewT
    """

    # If a figure is provided for single-panel plot, clear (0-1) axes decorations
    # If it is not a single-panel plot, you may not want to clear. See
    # input_sounding_era5.ipynb where multiple CM1 soundings are plotted.
    if fig is not None and (not subplot or subplot == (1, 1, 1)):
        fig.clear()

    # Validate required variables in dataset
    assert "surface_geopotential_height" in ds, (
        "skewt needs surface surface_geopotential_height"
    )
    assert "SP" in ds, "skewt needs surface pressure SP"
    # Load Dataset to avoid
    # plot_colormapped KeyError: 'Indexing with a boolean dask array is not allowed.
    # and allow .where function to check condition
    logging.info("load dataset and sort by descending pressure")
    ds = ds.sortby("P", ascending=False).load()
    old_nlevel = ds.level.size
    # Drop low pressure model levels to avoid ValueError: ODE Integration failed...
    # Mask high pressure levels greater than surface pressure SP
    mask = (ds.P >= 10 * units.hPa) & (ds.P < ds.SP)
    ds = ds.sel(level=mask)
    logging.info(
        f"After dropping low pressure levels, {ds.level.size}/{old_nlevel} remain"
    )

    # if args.model_levels, ds["P"] is 3D DataArray with vertical dim = model level.
    # Otherwise pressure is ds.level, a 1-D array of pressure.
    p = ds.P.data
    T = ds.T.data
    Td = mpcalc.dewpoint_from_specific_humidity(p, ds.Q.data)
    if any(Td > T):
        logging.warning("some Td > T")

    u = ds.U.data
    v = ds.V.data
    height = ds.Z.data
    sfc_hgt = ds.surface_geopotential_height.data
    agl = (height - sfc_hgt).to("m")  # So agl has .min() and .max() methods
    sp = ds.SP.data

    barb_increments = {"flag": 25, "full": 5, "half": 2.5}
    plot_barbs_units = "m/s"

    skew = SkewT(fig, subplot=subplot, rotation=rotation)
    skew.plot_dry_adiabats(lw=0.75, alpha=0.5)
    skew.plot_moist_adiabats(lw=0.75, alpha=0.25)
    skew.plot_mixing_lines(alpha=0.5)
    # 0-degC isotherm
    skew.ax.axvline(0, color="c", linestyle="--", linewidth=1.5, alpha=0.5)

    # Get parcel potential temperature and mixing ratio from
    # "surface_potential_temperature" and/or "surface_mixing_ratio" DataArray.
    # Otherwise, assume the value(s) at the first level.
    if "surface_potential_temperature" in ds:
        T_parcel = mpcalc.temperature_from_potential_temperature(
            sp, ds.surface_potential_temperature.data
        )
        logging.info(
            f"got T_parcel {T_parcel.to('degC'):~.2f} from SP {sp:~.1f} "
            f"and surface_potential_temperature {ds.surface_potential_temperature.data.to('degC'):~.2f}"
        )
    else:
        T_parcel = T[0]
        logging.info(f"got T_parcel {T_parcel:~.2f} from first level")
    if "surface_mixing_ratio" in ds:
        parcel_mixing_ratio = ds.surface_mixing_ratio.data
        logging.info(
            f"got parcel_mixing_ratio "
            f"{parcel_mixing_ratio.to('g/kg'):~.3f} "
            f"from surface_mixing_ratio"
        )
    else:
        parcel_mixing_ratio = mpcalc.mixing_ratio_from_specific_humidity(ds.Q.data[0])
        logging.info(
            f"got parcel_mixing_ratio {parcel_mixing_ratio.to('g/kg'):~.3f} "
            "from first level"
        )
    Td_parcel = mpcalc.dewpoint(mpcalc.vapor_pressure(sp, parcel_mixing_ratio))
    logging.info(
        f"p_parcel {sp:~.1f} T_parcel {T_parcel.to('degC'):~.2f} "
        f"Td_parcel {Td_parcel:~.2f} "
        f"parcel_mixing_ratio {parcel_mixing_ratio.to('g/kg'):~.3f}"
    )

    # Calculate LCL pressure and label level on SkewT.
    lcl_pressure, lcl_temperature = mpcalc.lcl(sp, T_parcel, Td_parcel)
    logging.info(f"lcl_p {lcl_pressure:~.1f} lcl_t {lcl_temperature:~.2f}")

    trans = transforms.blended_transform_factory(skew.ax.transAxes, skew.ax.transData)
    hpos = 0.82
    skew.ax.plot([hpos, hpos + 0.03], 2 * [lcl_pressure.m_as("hPa")], transform=trans)
    skew.ax.text(
        hpos + 0.035,
        lcl_pressure,
        f" LCL {lcl_pressure.to('hPa').item():~.0f}",
        transform=trans,
        horizontalalignment="left",
        verticalalignment="center",
        fontsize="x-small",
    )

    # Draw virtual temperature like temperature, but thin and dashed.
    if "Tv" in ds:
        Tv = ds["Tv"]
    else:
        logging.warning("Derive Tv from p, T, Td")
        Tv = mpcalc.virtual_temperature_from_dewpoint(p, T, Td)

    # Plot temperature and dewpoint.
    skew.plot(p, T, "r")
    skew.plot(p, Tv, "r", lw=0.5, linestyle="dashed")
    skew.plot(p, Td, "g")

    p_wLCL, T_wLCL, Td_wLCL, profT_wLCL = mpcalc.parcel_profile_with_lcl(p, T, Td)
    # parcel T and Tv virtual temperature
    profTv_wLCL = mpcalc.virtual_temperature_from_dewpoint(p_wLCL, profT_wLCL, Td_wLCL)
    skew.plot(p_wLCL, profT_wLCL, "k", linewidth=1.5, linestyle="dashed")
    skew.plot(p_wLCL, profTv_wLCL, "k", linewidth=0.5, linestyle="dashed")

    # Don't feed cape_cin virtual temperature. It assumes prof
    # is regular temperature, not virtual. It converts to virtual
    # temperature on its own.
    cape, cin = mpcalc.surface_based_cape_cin(p, T, Td)
    # Shade areas of CAPE and CIN
    Tv_wLCL = mpcalc.virtual_temperature_from_dewpoint(p_wLCL, T_wLCL, Td_wLCL)
    skew.shade_cin(p_wLCL, Tv_wLCL, profTv_wLCL)
    skew.shade_cape(p_wLCL, Tv_wLCL, profTv_wLCL)

    # Good bounds for aspect ratio
    skew.ax.set_xlim(xlim)
    skew.ax.set_ylim(None, ptop)

    title = ""
    if "time" in ds:
        title += f"{ds.time.dt.strftime('%Y-%m-%d %H:%M:%S').item()} "
    else:
        logging.warning("no 'time' variable in sounding")
    if "longitude" in ds:
        title += f"{ds.longitude.item():.3f} {ds.latitude.item():.3f}"
    if cape is not None:
        title += f"\ncape={cape:~.0f}   cin={cin:~.0f}   "
    skew.ax.set_title(title, fontsize="x-small")

    label_hgts = np.array([0, 1, 3, 6, 9, 12, 15]) * units.km
    # Label AGL intervals along the y-axis of the skewT.
    for label_hgt in label_hgts:
        if label_hgt >= agl.min() and label_hgt <= agl.max():
            (agl2p,) = interpolate_1d(label_hgt, agl, p)
        s = f"{label_hgt:~.0f}"
        if label_hgt == 0 * units.km:
            agl2p = sp.item()
            s = f"SFC {ds.surface_geopotential_height.item():~.0f}"
        skew.ax.plot([0, 0.01], 2 * [agl2p.m_as("hPa")], transform=trans, color="k")
        skew.ax.text(
            0.015,
            agl2p,
            s,
            transform=trans,
            horizontalalignment="left",
            verticalalignment="center",
            fontsize="x-small",
            clip_on=True,
        )

    skip_winds = not u.any() and not v.any()
    if skip_winds:
        return skew

    logging.info("work on winds and kinematics")
    storm_u = 0.0 * units("m/s")
    storm_v = 0.0 * units("m/s")
    right_mover, left_mover, wind_mean = mpcalc.bunkers_storm_motion(p, u, v, height)
    storm_u, storm_v = wind_mean
    srh03_pos, srh03_neg, srh03_tot = mpcalc.storm_relative_helicity(
        height, u, v, 3 * units.km, storm_u=storm_u, storm_v=storm_v
    )

    title += f"storm_u={storm_u:~.1f}   storm_v={storm_v:~.1f}"
    title += f"\n0-3km srh+={srh03_pos:~.0f}   srh-={srh03_neg:~.0f}   srh(tot)={srh03_tot:~.0f}"

    agl_colors = ["red", "red", "lime", "green", "blueviolet", "cyan"]
    # Find corresponding indices
    indices = np.digitize(agl.m_as("m"), label_hgts.m_as("m"))
    barbcolor = []
    for idx in indices:
        if idx < 0:
            barbcolor.append("#00000050")  # AGL below zero
        elif idx < len(agl_colors):
            barbcolor.append(agl_colors[idx])
        else:
            barbcolor.append("none")  # above max label

    skew.plot_barbs(
        p,
        u,
        v,
        barbcolor=barbcolor,
        length=6,
        plot_units=plot_barbs_units,
        linewidth=0.6,
        xloc=1.05,
        barb_increments=barb_increments,
    )
    skew.ax.text(
        1.07,
        0,
        f"{plot_barbs_units}",  # \n{barb_increments}",
        va="top",
        ha="right",
        transform=skew.ax.transAxes,
        fontsize="xx-small",
    )

    logging.info("Create hodograph")
    ax_hod = inset_axes(skew.ax, "40%", "40%", loc=1)
    # Initialize Hodograph without a fixed component range to allow for dynamic axes
    h = Hodograph(ax_hod)
    h.add_grid(increment=10, linewidth=0.75)

    # Calculate the data range for u and v winds to center the plot
    u_min, u_max = u.min(), u.max()
    v_min, v_max = v.min(), v.max()

    # Determine the center of the wind data
    u_center = (u_min + u_max) / 2
    v_center = (v_min + v_max) / 2

    # To make the plot square, find the largest range (u or v) and apply it to both axes
    max_range = max(u_max - u_min, v_max - v_min).to(plot_barbs_units).m
    # Add 10% padding around the maximum range
    plot_radius = (max_range / 2) * 1.1

    # Set the x and y limits of the hodograph axes, centered on the data
    ax_hod.set_xlim(u_center.m - plot_radius, u_center.m + plot_radius)
    ax_hod.set_ylim(v_center.m - plot_radius, v_center.m + plot_radius)
    # ax_hod.set_xlabel("")
    ax_hod.set_ylabel("")
    ax_hod.xaxis.label.set_fontsize("xx-small")
    ax_hod.tick_params(axis="both", which="major", labelsize="xx-small")

    # Label AGL intervals in hodograph.
    for label_hgt in label_hgts:
        ax_hod.text(
            np.interp(label_hgt, agl, u),
            np.interp(label_hgt, agl, v),
            label_hgt.to("km").m,
            fontsize=7,
        )

    h.plot_colormapped(
        u,
        v,
        agl,
        intervals=label_hgts,
        colors=agl_colors,
    )
    ax_hod.plot(storm_u, storm_v, "x")
    skew.ax.set_title(title, fontsize="x-small")

    return skew


if __name__ == "__main__":
    main()
