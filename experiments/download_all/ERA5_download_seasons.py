"""
Download ERA5 (6-hourly) and build one file per season (DJF, MAM, JJA, SON), for two datasets:
    pl: 850 hPa u, v and vorticity (global)          -> <out_root>/<season>/<season>_0p25.nc
    sl: 10 m u, v, mean sea level pressure and gust  -> <out_root>/<season>/<season>_sl_0p25.nc
        (limited area configs.SL_AREA, used for sampling along the tracks)

Season folders (seasons.py): 19981999_DJF (Dec 1998 + Jan/Feb 1999), 1999_MAM, 1999_JJA, 1999_SON.
Every season file covers the season plus SEASON_BUFFER_DAYS (configs.py, 10 days) on both sides,
e.g. 19981999_DJF: 1998-11-21 00:00 -> 1999-03-10 18:00. The buffer lets the tracker follow storms
across the season boundary, so a storm can be assigned to the season holding its peak vorticity
(postprocess_on_seasons.py). The sl files get the same buffer, so every point of a kept track can be sampled.
The pl files are regridded to the resolutions in regrid_era5.RESOLUTIONS (sl only with --regrid-sl).

Steps:
1. Download one request per dataset and month (only the days the seasons still need) to <out_root>/raw/.
   Monthly requests follow the ERA5 storage on the CDS (fastest) and keep each file small.
   Several requests run at the same time, as the CDS time is mostly queueing.
2. Merge the months of each season (with the buffer), check every time step and write the season file.
   Raw months are deleted as soon as every season of the run that uses them is written.
3. Regrid each season (regrid_era5.regrid_season).

Already downloaded / merged / regridded files are skipped, so a run can be restarted.

Examples:
    # test on a small subset (MAM 1998, days 1, 2 and 25 of each month)
    python ERA5_download_seasons.py --first-season 1998 --last-season 1998 --season-types MAM --days 1 2 25 --out-root /tmp/era5_test
    # everything (all four seasons)
    python ERA5_download_seasons.py --out-root /Volumes/<drive>/ERA5_pressure_levels
    # only the single levels for seasons that already have the 850 hPa files
    python ERA5_download_seasons.py --datasets sl --out-root /Volumes/<drive>/ERA5_pressure_levels
"""

import argparse
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import pandas as pd
import xarray as xr

from regrid_era5 import NATIVE_RES, RESOLUTIONS, regrid_season, season_file
from seasons import (HOURS, SEASON_BUFFER_DAYS, SEASON_TYPES, buffered_range, months_needed,
                     season_list, season_name, season_times)

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "algorithms"))
from configs import SL_AREA


DEFAULT_OUT_ROOT = Path("/Users/asgerthormann/Documents/uni/master/master_code/master_thesis/reanalysis_data/ERA5/ERA5_pressure_levels")

# key -> CDS dataset, request template (year/month/day are added per request) and raw file prefix
DATASETS = {
    "pl": dict(dataset="reanalysis-era5-pressure-levels",
               request={"product_type": ["reanalysis"],
                        "variable": ["u_component_of_wind", "v_component_of_wind", "vorticity"],
                        "time": ["00:00", "06:00", "12:00", "18:00"],
                        "pressure_level": ["850"],
                        "data_format": "netcdf",
                        "download_format": "unarchived"},
               raw="era5_pl850"),
    # no total_precipitation: at 00/06/12/18 it is only a 1-hour accumulation, and mixing it with
    # instantaneous fields makes the CDS return a zip of two files
    "sl": dict(dataset="reanalysis-era5-single-levels",
               request={"product_type": ["reanalysis"],
                        "variable": ["10m_u_component_of_wind", "10m_v_component_of_wind",
                                     "mean_sea_level_pressure", "instantaneous_10m_wind_gust"],
                        "time": ["00:00", "06:00", "12:00", "18:00"],
                        "data_format": "netcdf",
                        "download_format": "unarchived",
                        "area": SL_AREA},
               raw="era5_sl"),
}
TIME = "valid_time"


def file_times(path):
    with xr.open_dataset(path) as ds:
        return pd.DatetimeIndex(ds[TIME].values)


def season_out(out_root, year, stype, key="pl"):
    return season_file(Path(out_root) / season_name(year, stype), NATIVE_RES, key)


# DOWNLOAD

def raw_file(out_root, year, month, key="pl"):
    return Path(out_root) / "raw" / f"{DATASETS[key]['raw']}_{year}_{month:02d}.nc"


def needed_months(out_root, seasons, days, key="pl"):
    """{(year, month): days} still needed to build the seasons whose 0.25 degree file of this dataset does not exist yet."""
    months = {}
    for y, t in seasons:
        if season_out(out_root, y, t, key).exists():
            continue
        for ym, d in months_needed(y, t, days).items():
            months.setdefault(ym, set()).update(d)
    return {ym: sorted(d) for ym, d in sorted(months.items())}


def download_month(year, month, days, out_root, key="pl", retries=3):
    """
    Download the given days of one month of one dataset to raw/<prefix>_<year>_<month>.nc.
    Skips if the file exists and holds all requested days (e.g. a kept raw file from an earlier run).
    """
    target = raw_file(out_root, year, month, key)
    wanted = pd.DatetimeIndex([pd.Timestamp(year, month, d, h) for d in days for h in HOURS])
    if target.exists():
        if wanted.difference(file_times(target)).empty:
            print(f"[download] SKIP {key} {year}-{month:02d} (exists)")
            return target
        print(f"[download] {key} {year}-{month:02d}: existing raw file lacks some days, downloading again")
    import cdsapi
    target.parent.mkdir(parents=True, exist_ok=True)
    part = target.with_suffix(".nc.part")

    request = dict(DATASETS[key]["request"], year=[str(year)], month=[f"{month:02d}"],
                   day=[f"{d:02d}" for d in days])

    for attempt in range(1, retries + 1):
        print(f"[download] START {key} {year}-{month:02d} (days {days[0]}-{days[-1]}, {len(days)} days, attempt {attempt})")
        try:
            client = cdsapi.Client(quiet=True, progress=False)
            client.retrieve(DATASETS[key]["dataset"], request).download(str(part))
            part.replace(target)
            print(f"[download] OK {key} {year}-{month:02d} -> {target} ({target.stat().st_size / 1e9:.2f} GB)")
            return target
        except Exception as e:
            print(f"[download] FAILED {key} {year}-{month:02d} (attempt {attempt}): {e}")
            part.unlink(missing_ok=True)
            if attempt < retries:
                time.sleep(60 * attempt)
    raise RuntimeError(f"Download of {key} {year}-{month:02d} failed after {retries} attempts")


def download_all(tasks, out_root, workers):
    """Download all (key, (year, month), days) tasks concurrently. Returns the (key, year, month) that failed."""
    failed = []
    with ThreadPoolExecutor(max_workers=workers) as pool:
        futures = {pool.submit(download_month, y, m, d, out_root, key): (key, y, m) for key, (y, m), d in tasks}
        for f in as_completed(futures):
            try:
                f.result()
            except Exception as e:
                print(f"[download] ERROR {e}")
                failed.append(futures[f])
    return sorted(failed)


# MERGE

def check_existing(path, year, stype, days):
    """Warn if an existing season file does not have the expected (buffered) time steps, e.g. an old DJF file."""
    times = file_times(path)
    expected = season_times(year, stype, days)
    if not times.equals(expected):
        print(f"[merge] WARNING {path.name} exists but has {len(times)} steps {times[0]} -> {times[-1]}, "
              f"expected {len(expected)} steps {expected[0]} -> {expected[-1]} "
              f"(with the {SEASON_BUFFER_DAYS}-day buffer). Delete the season's files (and its tracks) to rebuild it.")


def merge_season(year, stype, days, out_root, key="pl"):
    """Combine the months of one buffered season of one dataset into one checked season file."""
    out = season_out(out_root, year, stype, key)
    name = season_name(year, stype)
    if out.exists():
        print(f"[merge] SKIP {key} {name} (exists)")
        check_existing(out, year, stype, days)
        return out
    out.parent.mkdir(parents=True, exist_ok=True)

    expected = season_times(year, stype, days)
    start, end = expected[0], expected[-1]
    parts = [xr.open_dataset(raw_file(out_root, y, m, key)) for y, m in months_needed(year, stype, days)]
    try:
        ds = xr.concat([p.sel({TIME: slice(start, end)}) for p in parts], dim=TIME).sortby(TIME)
        ds = ds.drop_vars([c for c in ("expver", "number") if c in ds.variables])
        # make sure cdo recognises the grid and does not look for the dropped variables
        for v in ds.variables:
            ds[v].encoding.pop("coordinates", None)
        ds["latitude"].attrs.setdefault("units", "degrees_north")
        ds["longitude"].attrs.setdefault("units", "degrees_east")

        # check first/last time step and that no time step is missing or duplicated
        times = pd.DatetimeIndex(ds[TIME].values)
        if not times.equals(expected):
            missing = expected.difference(times)
            raise ValueError(
                f"Season {name} ({key}): got {len(times)} steps {times[0]} -> {times[-1]}, "
                f"expected {len(expected)} steps {expected[0]} -> {expected[-1]} "
                f"({len(missing)} missing, first: {list(missing[:3])})")

        encoding = {v: {"zlib": True, "complevel": 1} for v in ds.data_vars}
        tmp = out.with_suffix(".part.nc")
        ds.to_netcdf(tmp, encoding=encoding)
    finally:
        for p in parts:
            p.close()
    tmp.rename(out)
    print(f"[merge] OK {key} {name}: {times[0]:%Y-%m-%dT%H} -> {times[-1]:%Y-%m-%dT%H} "
          f"({len(times)} steps, incl. {SEASON_BUFFER_DAYS}-day buffer) -> {out}")
    return out


def remove_unused_raw(out_root, users, done, key="pl"):
    """Delete a raw month once every season of this run that uses it is written. users: {(y, m): {seasons}}"""
    for f in sorted((Path(out_root) / "raw").glob(f"{DATASETS[key]['raw']}_*_*.nc")):
        y, m = (int(x) for x in f.stem.split("_")[-2:])
        if users.get((y, m)) and users[(y, m)] <= done:
            f.unlink()
            print(f"[cleanup] removed {f.name}")


# MAIN

def main():
    parser = argparse.ArgumentParser(description="Download ERA5 (850 hPa and single levels) and build season files.")
    parser.add_argument("--out-root", type=Path, default=DEFAULT_OUT_ROOT,
                        help="folder for the season folders (use an external drive for the full range)")
    parser.add_argument("--first-season", type=int, default=1940, help="first season year (December year for DJF)")
    parser.add_argument("--last-season", type=int, default=2025, help="last season year (December year for DJF)")
    parser.add_argument("--season-types", nargs="+", choices=SEASON_TYPES, default=SEASON_TYPES,
                        help="season types to build (default: all four)")
    parser.add_argument("--days", type=int, nargs="+", default=list(range(1, 32)),
                        help="days of the month to download (use a few for testing)")
    parser.add_argument("--workers", type=int, default=4, help="number of CDS requests at the same time")
    parser.add_argument("--threads", type=int, default=4, help="cdo -P threads for regridding")
    parser.add_argument("--datasets", nargs="+", choices=list(DATASETS), default=list(DATASETS),
                        help="pl: 850 hPa u, v, vo (global), sl: 10 m u, v, msl and gust (configs.SL_AREA)")
    parser.add_argument("--regrid", action=argparse.BooleanOptionalAction, default=True,
                        help="regrid the 850 hPa files of every season after merging")
    parser.add_argument("--regrid-sl", action="store_true",
                        help="also regrid the single-level files (cut back to configs.SL_AREA)")
    parser.add_argument("--keep-raw", action="store_true", help="do not delete the monthly raw files")
    args = parser.parse_args()

    out_root = args.out_root
    seasons = season_list(range(args.first_season, args.last_season + 1), args.season_types)
    days = sorted(set(args.days))
    out_root.mkdir(parents=True, exist_ok=True)
    first, last = seasons[0], seasons[-1]
    print(f"{len(seasons)} seasons {season_name(*first)} -> {season_name(*last)} "
          f"({buffered_range(*first)[0]:%Y-%m-%d} -> {buffered_range(*last)[1]:%Y-%m-%d} with the "
          f"{SEASON_BUFFER_DAYS}-day buffer), days {days[0]}-{days[-1]} ({len(days)}), output in {out_root}")

    # 1. download
    tasks = [(key, ym, d) for key in args.datasets for ym, d in needed_months(out_root, seasons, days, key).items()]
    print(f"[download] {len(tasks)} monthly requests")
    failed_downloads = download_all(tasks, out_root, args.workers)

    # which seasons of this run use each raw month (for the cleanup)
    users = {}
    for s in seasons:
        for ym in months_needed(*s, days):
            users.setdefault(ym, set()).add(s)

    # 2. merge (+ 3. regrid), per dataset
    regrid = {"pl": args.regrid, "sl": args.regrid_sl}
    area = {"pl": None, "sl": SL_AREA}
    done = {key: set() for key in args.datasets}
    failed = []
    for s in seasons:
        season_dir = out_root / season_name(*s)
        for key in args.datasets:
            if any((key, y, m) in failed_downloads for y, m in months_needed(*s, days)):
                failed.append(f"{key} {season_name(*s)}")
                continue
            try:
                merge_season(*s, days, out_root, key)
                done[key].add(s)
                if not args.keep_raw:
                    remove_unused_raw(out_root, users, done[key], key)
                if regrid[key]:
                    regrid_season(season_dir, threads=args.threads, kind=key, area=area[key])
            except Exception as e:
                print(f"[error] {key} season {season_name(*s)}: {e}")
                failed.append(f"{key} {season_name(*s)}")

    for key in args.datasets:
        res_list = RESOLUTIONS if regrid[key] else [NATIVE_RES]
        n_files = sum(season_file(out_root / season_name(*s), r, key).exists() for s in done[key] for r in res_list)
        print(f"[done] {key}: {len(done[key])} seasons, {n_files} files in {out_root} "
              f"(expected {len(done[key]) * len(res_list)})")
    if failed:
        print(f"[done] FAILED: {failed}  (failed downloads: {failed_downloads}) - rerun to retry")


if __name__ == "__main__":
    main()
