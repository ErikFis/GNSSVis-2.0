"""Spherical beam offsets and a smooth chromatic elliptical power response."""

import numpy as np


def beam_offsets_deg(alt_deg, az_deg, pointing_alt_deg, pointing_az_deg):
    """Spherical logarithmic map about the boresight (no azimuth-wrap jump).

    The horizontal axis follows increasing azimuth at boresight; the vertical
    axis follows increasing elevation. The norm equals great-circle distance.
    """
    a, z = np.radians(alt_deg), np.radians(az_deg)
    a0, z0 = np.radians(pointing_alt_deg), np.radians(pointing_az_deg)
    dz = z - z0
    x = np.cos(a) * np.sin(dz)
    y = np.sin(a) * np.cos(a0) - np.cos(a) * np.sin(a0) * np.cos(dz)
    den = np.sin(a) * np.sin(a0) + np.cos(a) * np.cos(a0) * np.cos(dz)
    sep = np.arctan2(np.hypot(x, y), den)
    norm = np.hypot(x, y)
    factor = np.divide(sep, norm, out=np.ones_like(sep), where=norm > 1e-14)
    return np.degrees(x * factor), np.degrees(y * factor)


def beam_gain(offset_x_deg, offset_y_deg, frequency_hz, beam):
    """Power gain; widths are full widths at half maximum at nu_ref.

    Frequencies scale both widths as nu_ref / nu (fixed aperture). No hard
    selection gate truncates emission at the half-power contour. The spherical
    tangent-plane approximation is intended for this illustrative broad beam.
    """
    ratio = np.asarray(frequency_hz) / float(beam["reference_frequency_hz"])
    angle = np.radians(float(beam.get("position_angle_deg", 0)))
    x = np.asarray(offset_x_deg) * np.cos(angle) + np.asarray(offset_y_deg) * np.sin(
        angle
    )
    y = -np.asarray(offset_x_deg) * np.sin(angle) + np.asarray(offset_y_deg) * np.cos(
        angle
    )
    q = (x / float(beam["full_width_deg"])) ** 2 + (
        y / float(beam["minor_width_deg"])
    ) ** 2
    return np.exp(-4 * np.log(2) * q * ratio**2)
