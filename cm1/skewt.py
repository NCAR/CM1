"""
Plot ERA5 Skew-T and hodograph for a specified time and location.

This script generates a Skew-T diagram and a hodograph for ERA5 model data at
a user-specified time and location. It includes temperature, dewpoint, and wind
data, as well as dry adiabats, moist adiabats, and mixing lines on the Skew-T diagram.
"""

import logging
import os

import matplotlib.pyplot as plt
import pandas as pd

from cm1.input.sounding import era5_aws, era5_model_level, Sounding, get_ofile
from cm1.utils import parse_args


def main() -> None:
    """
    Main function to load ERA5 data for a specified time and location,
    and plot a Skew-T diagram with hodograph.

    Parameters are read from the command-line arguments using `parse_args()`.
    """

    args = parse_args()
    valid_time = pd.to_datetime(args.time)
    ofile = get_ofile(args)
    logging.debug(f"ofile: {ofile}")
    if os.path.exists(ofile):
        logging.info(f"read cached data from: {ofile}")
        ds = Sounding(ofile)
    else:
        if os.path.exists("/glade/campaign"):
            ds = era5_model_level(
                valid_time,
                args.lat,
                args.lon,
            )
        else:
            logging.warning("No campaign storage. Fetching pres lvl data from AWS")
            ds = era5_aws(valid_time, args.lat, args.lon)

        logging.info(f"Caching data to: {ofile}")
        ds.metpy.dequantify().to_netcdf(ofile)

    ds.plot()
    plt.show()


if __name__ == "__main__":
    main()
