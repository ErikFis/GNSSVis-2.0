"""
TLE Repository Management Utilities

Provides functions for managing the TLE repository, including:
- Archiving old TLEs
- Cleaning up outdated files
- Checking TLE age and validity
"""

from pathlib import Path
from datetime import datetime, timedelta
import shutil
from typing import List, Optional, Tuple


def tle_line_checksum(line: str) -> int:
    """Return the modulo-10 checksum of the first 68 columns of a TLE line.

    Digits count as their value, minus signs count as 1, everything else as 0.
    """
    total = 0
    for char in line[:68]:
        if char.isdigit():
            total += int(char)
        elif char == "-":
            total += 1
    return total % 10


def parse_tle_text(text: str, validate_checksum: bool = True) -> List[Tuple[str, str, str]]:
    """Parse name/line-1/line-2 TLE records from text.

    Args:
        text: TLE text in three-line (name, line 1, line 2) format.
        validate_checksum: Raise ValueError when a line checksum is wrong.

    Returns:
        List of (name, line1, line2) tuples in file order. Names have their
        surrounding whitespace stripped; element lines keep their columns.
    """
    lines = [ln.rstrip() for ln in text.splitlines() if ln.strip()]
    records: List[Tuple[str, str, str]] = []
    i = 0
    while i + 2 < len(lines):
        name, line1, line2 = lines[i], lines[i + 1], lines[i + 2]
        if line1.startswith("1 ") and line2.startswith("2 "):
            if validate_checksum:
                for line in (line1, line2):
                    if len(line) != 69 or not line[68].isdigit():
                        raise ValueError(f"Malformed TLE line for {name.strip()!r}: {line!r}")
                    if tle_line_checksum(line) != int(line[68]):
                        raise ValueError(f"TLE checksum mismatch for {name.strip()!r}: {line!r}")
                if line1[2:7] != line2[2:7]:
                    raise ValueError(f"TLE catalog numbers differ for {name.strip()!r}")
            records.append((name.strip(), line1, line2))
            i += 3
        else:
            i += 1
    return records


def read_tle_file(path, validate_checksum: bool = True) -> List[Tuple[str, str, str]]:
    """Read and parse a local three-line TLE file without network access."""
    return parse_tle_text(Path(path).read_text(encoding="utf-8"), validate_checksum)


def get_tle_repository_path() -> Path:
    """Get the path to the TLE repository."""
    # Assume this file is in src/gnssvis/utils.py
    project_root = Path(__file__).parent.parent.parent
    return project_root / "data" / "tle_repository"


def archive_current_tles(source_type: str = "gnss") -> Optional[Path]:
    """
    Archive current TLE files with today's date.
    
    Args:
        source_type: Type of TLE source ('gnss', 'active', etc.)
        
    Returns:
        Path to archived file, or None if no file to archive
    """
    repo = get_tle_repository_path()
    source_dir = repo / source_type
    archive_dir = repo / "archive"
    
    # Create archive directory structure: archive/YYYY/MM/
    today = datetime.now()
    year_month_dir = archive_dir / today.strftime("%Y") / today.strftime("%m")
    year_month_dir.mkdir(parents=True, exist_ok=True)
    
    # Find current TLE file
    current_tle = source_dir / f"{source_type}.txt"
    if not current_tle.exists():
        print(f"[ARCHIVE] No current TLE file found: {current_tle}")
        return None
    
    # Archive with date
    archive_name = f"tle_{today.strftime('%Y%m%d')}.txt"
    archive_path = year_month_dir / archive_name
    
    if archive_path.exists():
        print(f"[ARCHIVE] Already archived today: {archive_path}")
        return archive_path
    
    shutil.copy(current_tle, archive_path)
    print(f"[ARCHIVE] Archived TLE to: {archive_path}")
    return archive_path


def get_tle_age(tle_file: Path) -> Optional[int]:
    """
    Get age of TLE file in days.
    
    Args:
        tle_file: Path to TLE file
        
    Returns:
        Age in days, or None if file doesn't exist
    """
    if not tle_file.exists():
        return None
    
    file_time = datetime.fromtimestamp(tle_file.stat().st_mtime)
    age = (datetime.now() - file_time).days
    return age


def cleanup_old_archives(max_age_days: int = 365) -> int:
    """
    Remove archived TLE files older than specified days.
    
    Args:
        max_age_days: Maximum age in days before deletion
        
    Returns:
        Number of files deleted
    """
    repo = get_tle_repository_path()
    archive_dir = repo / "archive"
    
    if not archive_dir.exists():
        return 0
    
    deleted_count = 0
    cutoff_date = datetime.now() - timedelta(days=max_age_days)
    
    for tle_file in archive_dir.rglob("*.txt"):
        file_time = datetime.fromtimestamp(tle_file.stat().st_mtime)
        if file_time < cutoff_date:
            tle_file.unlink()
            deleted_count += 1
            print(f"[CLEANUP] Deleted old TLE: {tle_file}")
    
    # Remove empty directories
    for year_dir in archive_dir.iterdir():
        if year_dir.is_dir():
            for month_dir in year_dir.iterdir():
                if month_dir.is_dir() and not list(month_dir.iterdir()):
                    month_dir.rmdir()
            if not list(year_dir.iterdir()):
                year_dir.rmdir()
    
    print(f"[CLEANUP] Deleted {deleted_count} old TLE files")
    return deleted_count


def check_tle_freshness(tle_file: Path, warn_days: int = 7, error_days: int = 30) -> str:
    """
    Check if TLE file is fresh enough for accurate predictions.
    
    Args:
        tle_file: Path to TLE file
        warn_days: Age in days to trigger warning
        error_days: Age in days to trigger error
        
    Returns:
        Status string: 'fresh', 'warning', or 'error'
    """
    age = get_tle_age(tle_file)
    
    if age is None:
        return 'error'
    elif age < warn_days:
        return 'fresh'
    elif age < error_days:
        print(f"[WARNING] TLE file is {age} days old (accuracy may be reduced)")
        return 'warning'
    else:
        print(f"[ERROR] TLE file is {age} days old (accuracy significantly degraded)")
        return 'error'


def list_available_tles(source_type: str = "gnss") -> list:
    """
    List all available TLE files for a source type.
    
    Args:
        source_type: Type of TLE source
        
    Returns:
        List of tuples (date, path) for available TLEs
    """
    repo = get_tle_repository_path()
    source_dir = repo / source_type
    
    tles = []
    for tle_file in sorted(source_dir.glob(f"{source_type}_*.txt")):
        # Extract date from filename: gnss_20260204.txt
        date_str = tle_file.stem.split('_')[-1]
        try:
            date = datetime.strptime(date_str, "%Y%m%d")
            tles.append((date, tle_file))
        except ValueError:
            continue
    
    return tles


if __name__ == "__main__":
    # Example usage
    print("TLE Repository Management")
    print("=" * 50)
    
    # Archive current TLEs
    archive_current_tles("gnss")
    
    # Check freshness
    repo = get_tle_repository_path()
    gnss_tle = repo / "gnss" / "gnss.txt"
    status = check_tle_freshness(gnss_tle)
    print(f"TLE Status: {status}")
    
    # List available TLEs
    available = list_available_tles("gnss")
    print(f"\nAvailable TLEs: {len(available)}")
    for date, path in available[-5:]:  # Show last 5
        print(f"  {date.strftime('%Y-%m-%d')}: {path.name}")
