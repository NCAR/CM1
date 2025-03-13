import argparse
import logging
import os
from pathlib import Path
from typing import Optional, Tuple

import matplotlib.animation as animation
import matplotlib.pyplot as plt
import matplotlib.transforms as transforms
import metpy.calc as mpcalc
import numpy as np
import pandas as pd
import xarray
from IPython.display import HTML
from matplotlib.figure import Figure
from metpy.interpolate import interpolate_1d
from metpy.plots import Hodograph, SkewT
from metpy.units import units
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from pint import Quantity

CAMPAIGNDIR = "glade/campaign/collections/rda/data"
TMPDIR = Path(os.getenv("TMPDIR"))


def parse_args() -> argparse.Namespace:
    """
    Parse command-line arguments for ERA5 data retrieval.

    Returns
    -------
    argparse.Namespace
        Parsed command-line arguments including time, longitude, latitude,
        and optional path to glade directory.
    """
    parser = argparse.ArgumentParser(description="get ERA5 sounding at time, lon, lat")
    parser.add_argument("time", help="time")
    parser.add_argument(
        "lon",
        type=lambda x: float(x) % 360 * units.degreeE,
        help="longitude in degrees East",
    )
    parser.add_argument(
        "lat",
        type=lambda x: float(x) * units.degreeN,
        help="latitude in degrees North",
    )
    parser.add_argument("--glade", default="/", help="parent of glade directory")
    args = parser.parse_args()
    logging.info(args)
    return args


def animate_cm1out_nc(
    ds: xarray.Dataset, var_name: str, height: float, dim: str = "zh", interval: int = 200, **kwargs
):
    """
    Create an animation of a user-specified variable at a given vertical level.

    Parameters:
    - ds: xarray.Dataset
    - var_name: str, Name of the variable to animate
    - height: height in km
    - dim: str, Name of height dimension
    - interval: int, Interval between frames in milliseconds (default: 200ms)
    """

    # Extract the variable at the given vertical level
    data = ds[var_name].sel({dim: height}, method="nearest")

    time = pd.to_timedelta(ds.time)

    img = data.isel(time=0).plot.imshow(origin="lower", **kwargs)

    # Animation function
    def update(frame):
        img.set_array(data.isel(time=frame))
        img.axes.set_title(f"{var_name} at {data[dim].data:.2f} km, Time: {time[frame]}")
        return [img]

    # Create animation
    ani = animation.FuncAnimation(
        img.figure, update, frames=len(time), interval=interval, blit=False
    )

    # Display in notebook
    anim_html = HTML(ani.to_jshtml())
    plt.close(img.figure)  # close figure to prevent static image display
    return anim_html


def mean_lat_lon(lats_deg, lons_deg):
    """
    Calculates the mean latitude and longitude of a set of points on a sphere.

    Args:
        lats (list): List of latitudes in degrees.
        lons (list): List of longitudes in degrees.

    Returns:
        tuple: Mean latitude and longitude in degrees.
    """

    # Convert to radians
    lats = np.radians(lats_deg)
    lons = np.radians(lons_deg)

    # Calculate Cartesian coordinates
    x = np.cos(lats) * np.cos(lons)
    y = np.cos(lats) * np.sin(lons)
    z = np.sin(lats)

    # Calculate mean Cartesian coordinates
    x_mean = np.mean(x)
    y_mean = np.mean(y)
    z_mean = np.mean(z)

    # Convert back to spherical coordinates
    lon_mean = np.arctan2(y_mean, x_mean)
    hyp = np.sqrt(x_mean**2 + y_mean**2)
    lat_mean = np.arctan2(z_mean, hyp)

    # Convert back to degrees
    lat_mean = np.degrees(lat_mean)
    lon_mean = np.degrees(lon_mean)

    # Transfer attributes of inputs to outputs lat_mean and lon_mean.
    lat_mean = lat_mean * units.degrees_N
    lon_mean = lon_mean * units.degrees_E

    return lat_mean, lon_mean


def skewt(
    ds: xarray.Dataset,
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

    # Validate required variables in dataset
    assert (
        "surface_geopotential_height" in ds
    ), "skewt needs geopotential height at the surface surface_geopotential_height"
    assert "SP" in ds, "skewt needs surface pressure SP"
    # Load Dataset to avoid
    # plot_colormapped KeyError: 'Indexing with a boolean dask array is not allowed.
    # and allow .where function to check condition
    logging.info("load dataset")
    ds = ds.load()
    old_nlevel = ds.level.size
    # Drop low pressure model levels to avoid ValueError: ODE Integration failed...
    # Side effect of Dataset.where is that DataArrays without
    # a level dimension like SP are broadcast to all levels.
    # Apply mask only to variables with the 'level' dimension
    # Mask high pressure levels greater than surface pressure SP
    ds = xarray.Dataset(
        {
            var: (
                ds[var].where((ds.P >= 10 * units.hPa) & (ds.P < ds.SP), drop=True)
                if "level" in ds[var].dims
                else ds[var]
            )
            for var in ds
        }
    )

    logging.info(
        f"After dropping low pressure levels, {ds.level.size}/{old_nlevel} remain"
    )

    # if args.model_levels, ds["P"] is 3D DataArray with vertical dim = model level.
    # Otherwise pressure is ds.level, a 1-D array of pressure.
    height = ds["Z"]
    p = ds["P"]
    T = ds["T"]
    Td = mpcalc.dewpoint_from_specific_humidity(p, ds.Q)
    if any(Td > T):
        logging.warning("some Td > T")

    barb_increments = {"flag": 25, "full": 5, "half": 2.5}
    plot_barbs_units = "m/s"
    u = ds.U.metpy.convert_units(plot_barbs_units)
    v = ds.V.metpy.convert_units(plot_barbs_units)

    skew = SkewT(fig, subplot=subplot, rotation=rotation)
    # Add the relevant special lines
    skew.plot_dry_adiabats(lw=0.75, alpha=0.5)
    skew.plot_moist_adiabats(lw=0.75, alpha=0.25)
    skew.plot_mixing_lines(alpha=0.5)
    # Slanted line on 0 isotherm
    skew.ax.axvline(0, color="c", linestyle="--", linewidth=1.5, alpha=0.5)

    # Get parcel potential temperature and mixing ratio from
    # "surface_potential_temperature" and/or "surface_mixing_ratio" DataArray.
    # Otherwise, assume the value(s) of the level with the highest pressure.
    parcel_level = ds.level.sel(level=ds.P.compute().idxmax())
    if "surface_potential_temperature" in ds:
        T_parcel = mpcalc.temperature_from_potential_temperature(
            ds.SP, ds.surface_potential_temperature
        )
        logging.warning(
            f"got T_parcel {T_parcel.metpy.convert_units('degC').item():~.2f} from SP {ds.SP.item():~.1f} "
            f"and surface_potential_temperature {ds.surface_potential_temperature.metpy.convert_units('degC').item():~.2f}"
        )
    else:
        T_parcel = ds.T.sel(level=parcel_level)
        logging.warning(
            f"got T_parcel {T_parcel.item():~.2f} from "
            f"highest pressure level {parcel_level.item()}"
        )
    if "surface_mixing_ratio" in ds:
        parcel_mixing_ratio = ds.surface_mixing_ratio
        logging.warning(
            f"got parcel_mixing_ratio "
            f"{parcel_mixing_ratio.item().to('g/kg'):~.3f} "
            f"from surface_mixing_ratio"
        )
    else:
        parcel_mixing_ratio = mpcalc.mixing_ratio_from_specific_humidity(
            ds.Q.sel(level=parcel_level)
        )
        logging.warning(
            f"got parcel_mixing_ratio {parcel_mixing_ratio.item().to('g/kg'):~.3f} "
            f"from highest pressure level {parcel_level.item()}"
        )
    Td_parcel = mpcalc.dewpoint(mpcalc.vapor_pressure(ds.SP, parcel_mixing_ratio))
    logging.warning(
        f"p_parcel {ds.SP.item():~.1f} T_parcel {T_parcel.item().to('degC'):~.2f} "
        f"Td_parcel {Td_parcel.item():~.2f} "
        f"parcel_mixing_ratio {parcel_mixing_ratio.item().to('g/kg'):~.3f}"
    )

    # Calculate LCL pressure and label level on SkewT.
    lcl_pressure, lcl_temperature = mpcalc.lcl(ds.SP, T_parcel, Td_parcel)
    logging.warning(f"lcl_p {lcl_pressure:~.1f} lcl_t {lcl_temperature:~.2f}")

    trans = transforms.blended_transform_factory(skew.ax.transAxes, skew.ax.transData)
    skew.ax.plot(
        [0.82, 0.85], 2 * [lcl_pressure.m_as("hPa")], transform=trans, color="brown"
    )
    skew.ax.text(
        1,
        lcl_pressure,
        f"LCL {lcl_pressure.to('hPa').item():~.0f}",
        transform=trans,
        horizontalalignment="right",
        verticalalignment="center",
        color="brown",
        fontsize="x-small",
    )

    # Calculate full parcel profile with LCL
    p_without_lcl = p  # remember so we can interpolate other variables later
    if min(p) <= lcl_pressure < max(p):
        # Append LCL to pressure array.
        p = np.append(p.data, lcl_pressure)
    else:
        # Averaging temperature, pressure, and mixing ratio along levels can
        # make mixing ratio above saturation, making LCL below profile.
        # Don't bother adding lcl_pressure point in that case.
        logging.warning(f"lcl outside range of p {p.min().item():~.1f} {p.max().item():~.1f}")
        p = p.data  # convert to Quantity array (with units) so we can sort it
    # Create reverse sorted array of pressure. mpcalc assumes bottom-up arrays.
    p = np.sort(p)[::-1]
    # Interpolate other variables to p array (which now includes LCL).
    T, Td, u, v, height = interpolate_1d(
        p,
        p_without_lcl,
        T,
        Td,
        u,
        v,
        height,
    )
    profT = mpcalc.parcel_profile(p, T_parcel, Td_parcel)
    prof_mixing_ratio = mpcalc.saturation_mixing_ratio(p, profT)
    prof_mixing_ratio[p >= lcl_pressure] = (
        parcel_mixing_ratio.item()
    )  # unsaturated mixing ratio (constant up to LCL)
    # parcel virtual temperature
    profTv = mpcalc.virtual_temperature(profT, prof_mixing_ratio)

    # Plot temperature and dewpoint.
    skew.plot(p, T, "r")
    skew.plot(p, Td, "g")
    # Draw virtual temperature like temperature, but thin and dashed.
    if "Tv" in ds:
        logging.warning("Ignore input Tv. Derive from T, qv(p,Td(p,Q))")

    Tv = mpcalc.virtual_temperature(T, mpcalc.saturation_mixing_ratio(p, Td))

    # Environment virtual temperature
    skew.plot(p, Tv, "r", lw=0.5, linestyle="dashed")
    skew.plot(p, profT, "k", linewidth=1.5, linestyle="dashed")
    # Parcel virtual temperature
    skew.plot(p, profTv, "k", linewidth=0.5, linestyle="dashed")

    # Don't feed cape_cin virtual temperature. It assumes prof
    # is regular temperature, not virtual. It converts to virtual
    # temperature on its own.
    cape, cin = mpcalc.cape_cin(p, T, Td, profT)
    # Shade areas of CAPE and CIN
    skew.shade_cin(p, Tv, profTv)
    skew.shade_cape(p, Tv, profTv)

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

    agl = height - ds.surface_geopotential_height.item()

    label_hgts = [0, 1, 3, 6, 9, 12, 15] * units("km")
    # Label AGL intervals along the y-axis of the skewT.
    for label_hgt in label_hgts:
        if label_hgt >= agl.min() and label_hgt <= agl.max():
            (agl2p,) = interpolate_1d(label_hgt, agl, p)
        s = f"{label_hgt:~.0f}"
        if label_hgt == 0 * units.km:
            agl2p = ds.SP.item()
            s = f"SFC {ds.surface_geopotential_height.item():~.0f}"
        skew.ax.plot([0, 0.01], 2 * [agl2p.m_as("hPa")], transform=trans, color="brown")
        skew.ax.text(
            0.01,
            agl2p,
            s,
            transform=trans,
            horizontalalignment="left",
            verticalalignment="center",
            color="brown",
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
    # Find corresponding colors
    barbcolor = []
    for h in agl:
        # Find the appropriate interval
        index = np.searchsorted(label_hgts, h, side="right") - 1
        if index == -1:
            # AGL below zero
            barbcolor.append("#00000050")
        elif index < len(agl_colors):
            barbcolor.append(agl_colors[index])
        else:
            barbcolor.append("none")

    bbz = skew.plot_barbs(
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
    h = Hodograph(ax_hod, component_range=30.0)
    h.add_grid(increment=10, linewidth=0.75)
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
