import io
import logging
import os
from pathlib import Path
from typing import Optional, TextIO, Tuple, Callable

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
from cm1.utils import get_ofile, parse_args

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

            # Ensure pressure is descending (Sfc -> Top)
            quantified_ds = quantified_ds.sortby("P", ascending=False)

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


def from_txt(source: str | os.PathLike) -> Sounding:
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


def ensure_quantity(val: float | Quantity, default_unit) -> Quantity:
    """
    Ensures a value is a Pint Quantity. If a raw float is passed,
    wraps it in default_unit and logs a warning.
    """
    if not hasattr(val, "units"):
        logging.warning(f"Input {val} missing units. Assuming {default_unit}.")
        return val * default_unit
    return val


def _fetch_era5_sounding(
    fetch_func: Callable,
    time: pd.Timestamp | np.datetime64 | str,
    lat: Quantity,
    lon: Quantity,
    **kwargs,
) -> Sounding:
    """
    Private helper to handle repeated ERA5 retrieval logic:
    unit validation, longitude normalization, and spatial selection.
    """
    lat = ensure_quantity(lat, units.degreeN)
    lon = ensure_quantity(lon, units.degreeE)

    valid_time = pd.Timestamp(time)
    ds = fetch_func(valid_time, **kwargs)

    # Ensure longitude is in 0-360 range for ERA5
    lon = lon % (360 * units.degreeE)

    ds = ds.sel(longitude=lon, latitude=lat, method="nearest", tolerance=5 * units.deg)
    ds = ds.sortby("P", ascending=False)

    return Sounding(ds)


def era5_aws(time, lat, lon, **kwargs) -> Sounding:
    """Retrieves ERA5 dataset for a specific time and location from AWS."""
    return _fetch_era5_sounding(cm1.input.era5.aws, time, lat, lon, **kwargs)


def era5_model_level(time, lat, lon, **kwargs) -> Sounding:
    """Retrieves ERA5 model-level dataset for a specific time and location."""
    return _fetch_era5_sounding(cm1.input.era5.model_level, time, lat, lon, **kwargs)


def era5_pressure_level(time, lat, lon, **kwargs) -> Sounding:
    """Retrieves ERA5 pressure-level dataset for a specific time and location."""
    return _fetch_era5_sounding(cm1.input.era5.pressure_level, time, lat, lon, **kwargs)


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
    rotation: int = 35,
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
    diff = Td - T

    # Identify and count exceedances
    exceedance = diff > 0
    num_affected = exceedance.sum()

    if num_affected > 0:
        max_exceedance = diff.max()

        logging.warning(
            f"Supersaturation detected in {num_affected} levels. "
            f"Max Td > T by {max_exceedance:.2f} °C. Clipping Td to T."
        )

        # "Where Td <= T, keep Td; otherwise, use T"
        Td = np.where(~exceedance, Td, T)
    else:
        logging.info("All levels are physically consistent (Td <= T).")

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
    # mixing_ratio_values from metpy default
    mixing_ratio = np.array([0.001, 0.002, 0.004, 0.007, 0.01, 0.016, 0.024, 0.032])
    lines = skew.plot_mixing_lines(mixing_ratio=mixing_ratio, alpha=0.5)
    segments = lines.get_segments()
    for i, seg in enumerate(segments):
        # seg is an array of [[temp, press], [temp, press], ...]
        # Find point with lowest pressure (top of line)
        top_idx = np.argmin(seg[:, 1])
        x_pos, y_pos = seg[top_idx]
        # kg/kg to g/kg
        label_val = f"{mixing_ratio[i] * 1000:.0f}"
        skew.ax.text(
            x_pos,
            y_pos,
            label_val,
            fontsize="xx-small",
            color="green",
            ha="center",
            va="bottom",
            alpha=0.5,
            clip_on=True,
        )

    # 0-degC isotherm
    skew.ax.axvline(0, color="c", linestyle="--", linewidth=1.5, alpha=0.5)

    # Draw virtual temperature like temperature, but thin and dashed.
    if "Tv" in ds:
        Tv = ds["Tv"].data
    else:
        logging.warning("Derive Tv from p, T, Td")
        Tv = mpcalc.virtual_temperature_from_dewpoint(p, T, Td)

    # Plot temperature and dewpoint.
    skew.plot(p, T, "r")
    skew.plot(p, Tv, "r", lw=0.5, linestyle="dashed")
    skew.plot(p, Td, "g")

    # Plot parcel T and Tv
    p_mu, t_mu, td_mu, idx_mu = mpcalc.most_unstable_parcel(p, T, Td)
    logging.info(
        f"most unstable parcel {p_mu:~.1f} {t_mu.to('degC'):~.2f} {td_mu:~.2f} i={idx_mu}"
    )
    # Calculate LCL pressure and label level on SkewT.
    lcl_p, lcl_temperature = mpcalc.lcl(p_mu, t_mu, td_mu)
    logging.info(f"lcl_p {lcl_p:~.1f} lcl_t {lcl_temperature:~.2f}")

    trans = transforms.blended_transform_factory(skew.ax.transAxes, skew.ax.transData)
    hpos = 0.82
    skew.ax.plot([hpos, hpos + 0.03], 2 * [lcl_p.m_as("hPa")], transform=trans)
    skew.ax.text(
        hpos + 0.035,
        lcl_p,
        f" LCL {lcl_p.to('hPa').item():~.0f}",
        transform=trans,
        horizontalalignment="left",
        verticalalignment="center",
        fontsize="x-small",
    )

    # Slice environment arrays to start at the MU level (the "Path" range)
    p_path = p[idx_mu:]
    Tv_env_path = Tv[idx_mu:]

    # Calculate Parcel Thermodynamic Path
    T_parcel_path = mpcalc.parcel_profile(p_path, t_mu, td_mu)

    # Moisture logic
    e_mu = mpcalc.saturation_vapor_pressure(td_mu)
    w_parcel_start = mpcalc.mixing_ratio(e_mu, p_mu)
    w_parcel_path = mpcalc.saturation_mixing_ratio(p_path, T_parcel_path)

    # Below LCL, mixing ratio is constant
    w_parcel_path[p_path >= lcl_p] = w_parcel_start
    Tv_parcel_path = mpcalc.virtual_temperature(T_parcel_path, w_parcel_path)

    # Plotting and Shading
    skew.plot(p_path, T_parcel_path, "k", ls="--")
    skew.plot(p_path, Tv_parcel_path, "k", lw=0.5, ls="--")

    # Shading (using the identical-length path arrays)
    skew.shade_cin(p_path, Tv_env_path, Tv_parcel_path)
    skew.shade_cape(p_path, Tv_env_path, Tv_parcel_path)

    skew.ax.set_xlim(xlim)
    skew.ax.set_ylim(None, ptop)

    title = ""
    if "time" in ds:
        title += f"{ds.time.dt.strftime('%Y-%m-%d %H:%M:%S').item()} "
    else:
        logging.warning("no 'time' variable in sounding")
    if "longitude" in ds:
        title += f"{ds.longitude.item():.3f} {ds.latitude.item():.3f}"
    cape, cin = mpcalc.most_unstable_cape_cin(p, T, Td)
    title += f"\nmucape={cape:~.0f}   mucin={cin:~.0f}   "
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
