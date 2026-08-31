#!/usr/bin/env python3

import sys
import numpy as np
from netCDF4 import Dataset


def main():
    if len(sys.argv) != 2:
        print(f"Usage: {sys.argv[0]} ncdiff.nc")
        sys.exit(1)

    filename = sys.argv[1]

    with Dataset(filename, "r") as ds:
        print(f"Checking: {filename}\n")

        for name, var in ds.variables.items():

            # Skip coordinate/dimension variables
            if name in ds.dimensions:
                continue

            try:
                data = var[:]
            except Exception as err:
                print(f"Could not read {name}: {err}")
                continue

            # Remove masked values
            if np.ma.isMaskedArray(data):
                data = data.compressed()
            else:
                data = np.asarray(data).ravel()

            if data.size == 0:
                continue

            # Skip non-numeric variables
            if not np.issubdtype(data.dtype, np.number):
                continue

            # Ignore NaNs when computing statistics
            finite = data[np.isfinite(data)]

            if finite.size == 0:
                continue

            max_abs = np.max(np.abs(finite))

            if max_abs != 0.0:
                print(
                    f"{name:35s} "
                    f"min={np.min(finite): .6e}  "
                    f"max={np.max(finite): .6e}  "
                    f"maxabs={max_abs: .6e}"
                )


if __name__ == "__main__":
    main()
