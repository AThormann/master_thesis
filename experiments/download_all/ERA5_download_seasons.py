"""
Download ERA5 (6-hourly) and build one file per DJF season, for two datasets:
    pl: 850 hPa u, v and vorticity (global)          -> <out_root>/<S><S+1>_DJF/<S><S+1>_DJF_0p25.nc
    sl: 10 m u, v and mean sea level pressure        -> <out_root>/<S><S+1>_DJF/<S><S+1>_DJF_sl_0p25.nc
        (limited area configs.SL_AREA, used for sampling along the tracks)

Season S = December of year S + January and February of year S+1.
The pl files are regridded to the resolutions in regrid_era5.RESOLUTIONS (sl only with --regrid-sl).

Steps:
1. Download one request per year (only the months a season still needs) to <out_root>/raw/.
   Several requests run at the same time, as the CDS time is mostly queueing.
2. Merge Dec S + Jan/Feb S+1, check the first/last timestep and write the season file.
   Yearly raw files are deleted as soon as both seasons using them are written.
3. Regrid each season (regrid_era5.regrid_season).

Already downloaded / merged / regridded files are skipped, so a run can be restarted.

Examples:
    # test on a small subset (two seasons, two days per month)
    python ERA5_download_seasons.py --first-season 1998 --last-season 1999 --days 1 2 --out-root /tmp/era5_test
    # everything
    python ERA5_download_seasons.py --out-root /Volumes/<drive>/ERA5_pressure_levels
    # only the single levels for seasons that already have the 850 hPa files
    python ERA5_download_seasons.py --datasets sl --out-root /Volumes/<drive>/ERA5_pressure_levels
"""

import argparse
import calendar
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import pandas as pd
import xarray as xr

from regrid_era5 import NATIVE_RES, RESOLUTIONS, regrid_season, season_file

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
HOURS = [0, 6, 12, 18]
TIME = "valid_time"


# DOWNLOAD

def raw_file(out_root, year, key="pl"):
    return Path(out_root) / "raw" / f"{DATASETS[key]['raw']}_{year}.nc"


def needed_months(out_root, seasons, key="pl"):
    """{year: months} still needed to build the seasons whose 0.25 degree file of this dataset does not exist yet."""
    months = {}
    for s in seasons:
        if season_file(Path(out_root) / f"{s}{s + 1}_DJF", NATIVE_RES, key).exists():
            continue
        months.setdefault(s, set()).add(12)
        months.setdefault(s + 1, set()).update([1, 2])
    return {y: sorted(m) for y, m in sorted(months.items())}


def download_year(year, months, days, out_root, key="pl", retries=3):
    """Download one year (the given months) of one dataset to raw/<prefix>_<year>.nc. Skips if it exists."""
    target = raw_file(out_root, year, key)
    if target.exists():
        print(f"[download] SKIP {key} year {year} (exists)")
        return target
    import cdsapi
    target.parent.mkdir(parents=True, exist_ok=True)
    part = target.with_suffix(".nc.part")

    request = dict(DATASETS[key]["request"],
                   year=[str(year)],
                   month=[f"{m:02d}" for m in months],
                   day=[f"{d:02d}" for d in days])

    for attempt in range(1, retries + 1):
        print(f"[download] START {key} year {year} (months {request['month']}, attempt {attempt})")
        try:
            client = cdsapi.Client(quiet=True, progress=False)
            client.retrieve(DATASETS[key]["dataset"], request).download(str(part))
            part.rename(target)
            print(f"[download] OK {key} year {year} -> {target} ({target.stat().st_size / 1e9:.2f} GB)")
            return target
        except Exception as e:
            print(f"[download] FAILED {key} year {year} (attempt {attempt}): {e}")
            part.unlink(missing_ok=True)
            if attempt < retries:
                time.sleep(60 * attempt)
    raise RuntimeError(f"Download of {key} year {year} failed after {retries} attempts")


def download_all(tasks, days, out_root, workers):
    """Download all (key, year, months) tasks concurrently. Returns the (key, year) pairs that failed."""
    failed = []
    with ThreadPoolExecutor(max_workers=workers) as pool:
        futures = {pool.submit(download_year, y, m, days, out_root, key): (key, y) for key, y, m in tasks}
        for f in as_completed(futures):
            try:
                f.result()
            except Exception as e:
                print(f"[download] ERROR {e}")
                failed.append(futures[f])
    return sorted(failed)


# MERGE

def expected_times(season, days):
    """All 6-hourly timesteps of the season, restricted to the requested days."""
    times = []
    for year, month in [(season, 12), (season + 1, 1), (season + 1, 2)]:
        n_days = calendar.monthrange(year, month)[1]
        for d in days:
            if d <= n_days:
                times += [pd.Timestamp(year, month, d, h) for h in HOURS]
    return pd.DatetimeIndex(times)


def merge_season(season, days, out_root, key="pl"):
    """Combine Dec <season> and Jan/Feb <season+1> of one dataset into one checked season file."""
    out = season_file(Path(out_root) / f"{season}{season + 1}_DJF", NATIVE_RES, key)
    if out.exists():
        print(f"[merge] SKIP {key} {season}/{season + 1} (exists)")
        return out
    out.parent.mkdir(parents=True, exist_ok=True)

    with xr.open_dataset(raw_file(out_root, season, key)) as ds_dec, \
         xr.open_dataset(raw_file(out_root, season + 1, key)) as ds_jf:
        dec = ds_dec.sel({TIME: ds_dec[TIME].dt.month == 12})
        jf = ds_jf.sel({TIME: ds_jf[TIME].dt.month.isin([1, 2])})
        ds = xr.concat([dec, jf], dim=TIME).sortby(TIME)
        ds = ds.drop_vars([c for c in ("expver", "number") if c in ds.variables])
        # make sure cdo recognises the grid and does not look for the dropped variables
        for v in ds.variables:
            ds[v].encoding.pop("coordinates", None)
        ds["latitude"].attrs.setdefault("units", "degrees_north")
        ds["longitude"].attrs.setdefault("units", "degrees_east")

        # check first/last timestep (Feb 28 or 29) and that no timestep is missing
        times = pd.DatetimeIndex(ds[TIME].values)
        expected = expected_times(season, days)
        if not times.equals(expected):
            missing = expected.difference(times)
            raise ValueError(
                f"Season {season}/{season + 1} ({key}): got {len(times)} steps {times[0]} -> {times[-1]}, "
                f"expected {len(expected)} steps {expected[0]} -> {expected[-1]} "
                f"({len(missing)} missing, first: {list(missing[:3])})")

        encoding = {v: {"zlib": True, "complevel": 1} for v in ds.data_vars}
        tmp = out.with_suffix(".part.nc")
        ds.to_netcdf(tmp, encoding=encoding)
    tmp.rename(out)
    print(f"[merge] OK {key} {season}/{season + 1}: {times[0]:%Y-%m-%dT%H} -> {times[-1]:%Y-%m-%dT%H} "
          f"({len(times)} steps) -> {out}")
    return out


def remove_unused_raw(out_root, seasons, done, key="pl"):
    """Delete raw year Y once season Y-1 (uses Jan/Feb) and season Y (uses Dec) are both written."""
    for f in sorted((Path(out_root) / "raw").glob(f"{DATASETS[key]['raw']}_*.nc")):
        year = int(f.stem.split("_")[-1])
        users = [s for s in (year - 1, year) if s in seasons]
        if all(s in done for s in users):
            f.unlink()
            print(f"[cleanup] removed {f.name}")


# MAIN

def main():
    parser = argparse.ArgumentParser(description="Download ERA5 (850 hPa and single levels) and build DJF season files.")
    parser.add_argument("--out-root", type=Path, default=DEFAULT_OUT_ROOT,
                        help="folder for the season folders (use an external drive for the full range)")
    parser.add_argument("--first-season", type=int, default=1940, help="first season (December year)")
    parser.add_argument("--last-season", type=int, default=2025, help="last season (December year)")
    parser.add_argument("--days", type=int, nargs="+", default=list(range(1, 32)),
                        help="days of the month to download (use a few for testing)")
    parser.add_argument("--workers", type=int, default=4, help="number of CDS requests at the same time")
    parser.add_argument("--threads", type=int, default=4, help="cdo -P threads for regridding")
    parser.add_argument("--datasets", nargs="+", choices=list(DATASETS), default=list(DATASETS),
                        help="pl: 850 hPa u, v, vo (global), sl: 10 m u, v and msl (configs.SL_AREA)")
    parser.add_argument("--regrid", action=argparse.BooleanOptionalAction, default=True,
                        help="regrid the 850 hPa files of every season after merging")
    parser.add_argument("--regrid-sl", action="store_true",
                        help="also regrid the single-level files (cut back to configs.SL_AREA)")
    parser.add_argument("--keep-raw", action="store_true", help="do not delete the yearly raw files")
    args = parser.parse_args()

    out_root = args.out_root
    seasons = list(range(args.first_season, args.last_season + 1))
    days = sorted(set(args.days))
    out_root.mkdir(parents=True, exist_ok=True)
    print(f"Seasons {seasons[0]}/{seasons[0] + 1} -> {seasons[-1]}/{seasons[-1] + 1}, "
          f"days {days[0]}-{days[-1]} ({len(days)}), output in {out_root}")

    # 1. download
    tasks = [(key, y, m) for key in args.datasets for y, m in needed_months(out_root, seasons, key).items()]
    failed_downloads = download_all(tasks, days, out_root, args.workers)

    # 2. merge (+ 3. regrid), per dataset
    regrid = {"pl": args.regrid, "sl": args.regrid_sl}
    area = {"pl": None, "sl": SL_AREA}
    done = {key: set() for key in args.datasets}
    failed = []
    for s in seasons:
        season_dir = out_root / f"{s}{s + 1}_DJF"
        for key in args.datasets:
            if (key, s) in failed_downloads or (key, s + 1) in failed_downloads:
                failed.append(f"{key} {s}")
                continue
            try:
                merge_season(s, days, out_root, key)
                done[key].add(s)
                if not args.keep_raw:
                    remove_unused_raw(out_root, seasons, done[key], key)
                if regrid[key]:
                    regrid_season(season_dir, threads=args.threads, kind=key, area=area[key])
            except Exception as e:
                print(f"[error] {key} season {s}/{s + 1}: {e}")
                failed.append(f"{key} {s}")

    for key in args.datasets:
        n_res = len(RESOLUTIONS) if regrid[key] else 1
        n_files = sum(season_file(out_root / f"{s}{s + 1}_DJF", r, key).exists()
                      for s in done[key] for r in (RESOLUTIONS if regrid[key] else [NATIVE_RES]))
        print(f"[done] {key}: {len(done[key])} seasons, {n_files} files in {out_root} "
              f"(expected {len(done[key]) * n_res})")
    if failed:
        print(f"[done] FAILED: {failed}  (failed downloads: {failed_downloads}) - rerun to retry")


if __name__ == "__main__":
    main()
