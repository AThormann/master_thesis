"""
Seasons (DJF, MAM, JJA, SON): folder names, time ranges and months, shared by all scripts in this folder.

Season (year S, type T):
    DJF  December S + January/February S+1   folder <S><S+1>_DJF, e.g. 19981999_DJF
    MAM  March-May S                         folder <S>_MAM,      e.g. 1999_MAM
    JJA  June-August S                       folder <S>_JJA
    SON  September-November S                folder <S>_SON

core range     : first and last 6-hourly step of the season itself (e.g. 1998-12-01 00 -> 1999-02-28 18)
buffered range : the core range extended by SEASON_BUFFER_DAYS (configs.py, 10 days) on both sides
                 (e.g. 1998-11-21 00 -> 1999-03-10 18). The season files (tracking field and single levels)
                 cover the buffered range; a track belongs to the season if its peak vorticity lies in the
                 core range (postprocess_on_seasons.py).
"""

import re
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "algorithms"))
from configs import SEASON_BUFFER_DAYS


SEASON_MONTHS = {"DJF": (12, 1, 2), "MAM": (3, 4, 5), "JJA": (6, 7, 8), "SON": (9, 10, 11)}
SEASON_TYPES = list(SEASON_MONTHS)
HOURS = [0, 6, 12, 18]
STEP = pd.Timedelta(hours=6)
_DIR_RE = re.compile(r"^(\d{4})(\d{4})?_(DJF|MAM|JJA|SON)$")


def season_name(year, stype):
    """(1998, 'DJF') -> '19981999_DJF', (1999, 'MAM') -> '1999_MAM'"""
    if stype not in SEASON_MONTHS:
        raise ValueError(f"unknown season type {stype}, choose from {SEASON_TYPES}")
    return f"{year}{year + 1}_DJF" if stype == "DJF" else f"{year}_{stype}"


def parse_season(name):
    """'19981999_DJF' -> (1998, 'DJF'), '1999_MAM' -> (1999, 'MAM'). Accepts a folder path too."""
    m = _DIR_RE.match(Path(name).name)
    if not m:
        raise ValueError(f"not a season folder name: {name}")
    return int(m.group(1)), m.group(3)


def core_range(year, stype):
    """(first, last) 6-hourly time step of the season itself."""
    months = SEASON_MONTHS[stype]
    last_year = year + 1 if stype == "DJF" else year
    start = pd.Timestamp(year, months[0], 1)
    end = pd.Timestamp(last_year, months[-1], 1) + pd.offsets.MonthBegin(1) - STEP
    return start, end


def buffered_range(year, stype, buffer_days=SEASON_BUFFER_DAYS):
    """Core range extended by buffer_days on both sides."""
    start, end = core_range(year, stype)
    buf = pd.Timedelta(days=buffer_days)
    return start - buf, end + buf


def season_times(year, stype, days=None, buffer_days=SEASON_BUFFER_DAYS):
    """All 6-hourly steps of the buffered range, restricted to the days of the month in `days` (None: all)."""
    start, end = buffered_range(year, stype, buffer_days)
    times = pd.date_range(start, end, freq=STEP)
    if days is not None:
        times = times[times.day.isin(list(days))]
    return times


def months_needed(year, stype, days=None, buffer_days=SEASON_BUFFER_DAYS):
    """{(year, month): sorted days} needed to build the buffered season."""
    out = {}
    for t in season_times(year, stype, days, buffer_days):
        out.setdefault((t.year, t.month), set()).add(t.day)
    return {k: sorted(v) for k, v in sorted(out.items())}


def season_list(years, types=SEASON_TYPES):
    """[(year, type)] in chronological order of the season start."""
    return sorted(((y, t) for y in years for t in types), key=lambda s: core_range(*s)[0])


def season_dirs(out_root, years=None, types=SEASON_TYPES):
    """Season folders in out_root: the given years x types, or all existing folders of these types."""
    out_root = Path(out_root)
    if years:
        return [out_root / season_name(y, t) for y, t in season_list(years, types)]
    dirs = [d for t in types for d in out_root.glob(f"*_{t}") if d.is_dir() and _DIR_RE.match(d.name)]
    return sorted(dirs, key=lambda d: core_range(*parse_season(d))[0])
