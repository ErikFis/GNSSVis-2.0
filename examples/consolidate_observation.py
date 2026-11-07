#!/usr/bin/env python3
"""Example of consolidating NPZ observation data into FITS for gnssvis."""

from gnssvis.io import IOManager

# Create the I/O manager
io_mgr = IOManager()

# Consolidate data from 2025-09-10
fits_path = io_mgr.consolidate_observation_from_npz(
    data_dir="/home/erik-limeira/gnssvis/data/2025-09-10",
    output_filename="observacao_2025-09-10.fits",
    date_obs="2025-09-10",
    ra=0.0,    # Pointing RA (degrees)
    dec=-2.0   # Pointing DEC (degrees)
)

print(f"\nFITS file created: {fits_path}")
print("You can now use this file with gnssvis to plot satellites and celestial objects!")
