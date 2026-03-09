import argparse
import logging
import os
from pathlib import Path

import matplotlib.animation as animation
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr
from IPython.display import HTML
from metpy.units import units

TMPDIR = Path(os.getenv("TMPDIR"))


def parse_args() -> argparse.Namespace:
    """
    Parse command-line arguments for ERA5 data retrieval.

    Returns
    -------
    argparse.Namespace
        Parsed command-line arguments including time, longitude, latitude
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
    args = parser.parse_args()
    logging.info(args)
    return args


def animate_cm1out_nc(
    data: xr.DataArray,
    interval: int = 200,
    **kwargs,
):
    """
    Create an animation of a user-specified 2D field over its time dimension.

    The user is responsible for selecting any other dimensions (e.g., vertical)
    before passing the data to this function.

    Parameters:
    - data: xr.DataArray. A DataArray with a 'time' coordinate, ready to be plotted.
      Example: ds['cref'] or ds['u'].sel(zh=1.0, method='nearest')
    - interval: int, Interval between frames in milliseconds (default: 200ms).
    - **kwargs: Additional keyword arguments passed to the plot function.
    """
    if "time" not in data.dims:
        raise ValueError("Input DataArray must have a 'time' dimension.")

    time = pd.to_timedelta(data.time)

    # --- SIMPLIFIED LOGIC ---
    # The title is built from the DataArray's attributes.
    title_parts = []
    if data.name:
        title_parts.append(data.name)

    # Check for a vertical coordinate to add its value to the title.
    for coord_name in ["zh", "z", "level", "pressure"]:
        if coord_name in data.coords and data[coord_name].size == 1:
            level_val = data[coord_name].item()
            units = data[coord_name].attrs.get("units", "")
            title_parts.append(f"at {level_val:.2f} {units}")
            break

    title_parts.append("Time: {time}")
    title_template = ", ".join(title_parts)

    # --- PLOTTING LOGIC (largely unchanged) ---
    img = data.isel(time=0).plot.imshow(origin="lower", **kwargs)

    # Animation function
    def update(frame):
        img.set_array(data.isel(time=frame))

        # Use the appropriate title template.
        current_time = time[frame]
        img.axes.set_title(title_template.format(time=current_time))
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
