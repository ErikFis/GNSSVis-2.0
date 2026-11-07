# TLE and cache usage

## Load trajectories

```python
from astropy.time import Time
from gnssvis import SkyAnalyzer

sky = SkyAnalyzer()
times = Time(["2025-10-27T20:00:00", "2025-10-27T21:00:00"], scale="utc")
trajectories = sky.get_gnss_trajectories(times)
```

Ephemerides are downloaded during initialization when absent. The trajectory
request may download Celestrak TLEs. Network failures are reported in English.
Use observation-appropriate orbital elements when validating historical data.

## Archive and inspect TLEs

```python
from gnssvis.tle_utils import (
    archive_current_tles,
    check_tle_freshness,
    cleanup_old_archives,
    get_tle_repository_path,
    list_available_tles,
)

archive_path = archive_current_tles("gnss")
repo = get_tle_repository_path()
status = check_tle_freshness(repo / "gnss" / "gnss.txt")
available = list_available_tles("gnss")
```

`check_tle_freshness` returns `fresh`, `warning`, or `error`, using default
thresholds of seven and 30 days. Missing files return `error`.
`archive_current_tles` copies the current file to a dated archive, returning
None when the source is absent. `list_available_tles` lists dated source files.

`cleanup_old_archives(max_age_days=365)` deletes old archive files and empty
archive directories. Run it only when those archived elements are expendable.

## Clear session caches

```python
sky.clear_cache()
sky.clear_pulsar_cache()
```

These methods reset in-memory catalog references without deleting files.
All runtime data remains excluded from Git even after it is regenerated.
