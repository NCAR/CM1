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

import cm1.input.era5
from cm1.input.sounding import get_ofile
from cm1.utils import CAMPAIGNDIR, parse_args, skewt


def main() -> None:
    """
    Main function to load ERA5 data for a specified time and location,
    and plot a Skew-T diagram with hodograph.

    Parameters are read from the command-line arguments using `parse_args()`.
    """

    args = parse_args()

    ofile = get_ofile(args)
    if os.path.exists(ofile):
        logging.warning(f"read {ofile}")
        with open(ofile, "rb") as file:
            ds = pickle.load(file)
    else:
        if not os.path.exists("/"+CAMPAIGNDIR):
            ds = cm1.input.era5.aws(
                pd.to_datetime(args.time)
            )
        else:
            ds = cm1.input.era5.model_level(
                pd.to_datetime(args.time),
                glade=args.glade,
            )
    with open(ofile, "wb") as file:
        pickle.dump(ds, file)

    logging.warning(f"select {args}")

    ds = ds.sel(
        longitude=args.lon,
        latitude=args.lat,
        method="nearest",
    )
    skewt(ds)
    plt.show()


if __name__ == "__main__":
    main()
