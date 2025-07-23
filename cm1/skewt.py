"""
Plot ERA5 Skew-T and hodograph for a specified time and location.

This script generates a Skew-T diagram and a hodograph for ERA5 model data at 
a user-specified time and location. It includes temperature, dewpoint, and wind
data, as well as dry adiabats, moist adiabats, and mixing lines on the Skew-T diagram.
"""

import logging
import os
import pickle

import matplotlib.pyplot as plt
import pandas as pd

from cm1.input.sounding import era5_aws, era5_model_level, get_ofile
from cm1.utils import parse_args, skewt


def main() -> None:
    """
    Main function to load ERA5 data for a specified time and location,
    and plot a Skew-T diagram with hodograph.

    Parameters are read from the command-line arguments using `parse_args()`.
    """

    args = parse_args()
    valid_time = pd.to_datetime(args.time)
    ofile = get_ofile(args)
    if os.path.exists(ofile):
        logging.warning(f"read {ofile}")
        with open(ofile, "rb") as file:
            ds = pickle.load(file)
    else:
        if os.path.exists("/glade/campaign"):
            ds = era5_model_level(
                valid_time,
                args.lat,
                args.lon,
            )
        else:
            logging.warning("No campaign storage. Get pressure level data from AWS")
            ds = era5_aws(valid_time, args.lat, args.lon)

    with open(ofile, "wb") as file:
        pickle.dump(ds, file)

    skewt(ds)
    plt.show()


if __name__ == "__main__":
    main()
