# Astronomical data repository

This is a local runtime data directory. Downloaded contents are not included
in the prepared source package and are ignored by Git.

| Directory | Contents |
| --- | --- |
| `ephemeris/` | Skyfield planetary ephemerides such as `de421.bsp` |
| `gnss/` | Current `gnss.txt` and dated `gnss_YYYYMMDD.txt` Celestrak TLEs |
| `pulsars/` | Persistent `atnf_catalog.json`, when supplied |
| `archive/YYYY/MM/` | TLE archives created with `tle_utils` |

`SkyAnalyzer` creates the needed directories on demand. DE421 is cached
long term. GNSS TLE files older than seven days are considered stale.
The current implementation downloads when no usable session cache is
available, including the first request in a new session. Dated copies retain
the downloaded elements for reproducibility. Pulsar downloads are disabled;
an existing persistent catalog can still be loaded.

See [usage](USAGE.md) for cache and archive utilities. Do not add downloaded
catalogs or ephemerides to task submissions.
