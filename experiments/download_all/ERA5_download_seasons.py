"""
Download ERA5 850 hPa u, v and vorticity (6-hourly) and build one file per DJF season.

Season S = December of year S + January and February of year S+1, stored as
    <out_root>/<S><S+1>_DJF/<S><S+1>_DJF_0p25.nc
after which every season is regridded to the resolutions in regrid_era5.RESOLUTIONS.

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
"""

import argparse
import calendar
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import pandas as pd
import xarray as xr

from regrid_era5 import NATIVE_RES, RESOLUTIONS, regrid_season, season_file


DEFAULT_OUT_ROOT = Path("/Users/asgerthormann/Documents/uni/master/master_code/master_thesis/reanalysis_data/ERA5/ERA5_pressure_levels")

DATASET = "reanalysis-era5-pressure-levels"
REQUEST = {
    "product_type": ["reanalysis"],
    "variable": ["u_component_of_wind", "v_component_of_wind", "vorticity"],
    "time": ["00:00", "06:00", "12:00", "18:00"],
    "pressure_level": ["850"],
    "data_format": "netcdf",
    "download_format": "unarchived",
}
HOURS = [0, 6, 12, 18]
TIME = "valid_time"


# DOWNLOAD

def raw_file(out_root, year):
    return Path(out_root) / "raw" / f"era5_pl850_{year}.nc"


def needed_months(out_root, seasons):
    """{year: months} still needed to build the seasons whose 0.25 degree file does not exist yet."""
    months = {}
    for s in seasons:
        if season_file(Path(out_root) / f"{s}{s + 1}_DJF", NATIVE_RES).exists():
            continue
        months.setdefault(s, set()).add(12)
        months.setdefault(s + 1, set()).update([1, 2])
    return {y: sorted(m) for y, m in sorted(months.items())}


def download_year(year, months, days, out_root, retries=3):
    """Download one year (the given months) to raw/era5_pl850_<year>.nc. Skips if it exists."""
    target = raw_file(out_root, year)
    if target.exists():
        print(f"[download] SKIP year {year} (exists)")
        return target
    import cdsapi
    target.parent.mkdir(parents=True, exist_ok=True)
    part = target.with_suffix(".nc.part")

    request = dict(REQUEST,
                   year=[str(year)],
                   month=[f"{m:02d}" for m in months],
                   day=[f"{d:02d}" for d in days])

    for attempt in range(1, retries + 1):
        print(f"[download] START year {year} (months {request['month']}, attempt {attempt})")
        try:
            client = cdsapi.Client(quiet=True)
            client.retrieve(DATASET, request).download(str(part))
            part.rename(target)
            print(f"[download] OK year {year} -> {target} ({target.stat().st_size / 1e9:.2f} GB)")
            return target
        except Exception as e:
            print(f"[download] FAILED year {year} (attempt {attempt}): {e}")
            part.unlink(missing_ok=True)
            if attempt < retries:
                time.sleep(60 * attempt)
    raise RuntimeError(f"Download of year {year} failed after {retries} attempts")


def download_all(year_months, days, out_root, workers):
    """Download all years concurrently. Returns the years that failed."""
    failed = []
    with ThreadPoolExecutor(max_workers=workers) as pool:
        futures = {pool.submit(download_year, y, m, days, out_root): y for y, m in year_months.items()}
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


def merge_season(season, days, out_root):
    """Combine Dec <season> and Jan/Feb <season+1> into one checked season file."""
    out = season_file(Path(out_root) / f"{season}{season + 1}_DJF", NATIVE_RES)
    if out.exists():
        print(f"[merge] SKIP {season}/{season + 1} (exists)")
        return out
    out.parent.mkdir(parents=True, exist_ok=True)

    with xr.open_dataset(raw_file(out_root, season)) as ds_dec, \
         xr.open_dataset(raw_file(out_root, season + 1)) as ds_jf:
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
                f"Season {season}/{season + 1}: got {len(times)} steps {times[0]} -> {times[-1]}, "
                f"expected {len(expected)} steps {expected[0]} -> {expected[-1]} "
                f"({len(missing)} missing, first: {list(missing[:3])})")

        encoding = {v: {"zlib": True, "complevel": 1} for v in ds.data_vars}
        tmp = out.with_suffix(".part.nc")
        ds.to_netcdf(tmp, encoding=encoding)
    tmp.rename(out)
    print(f"[merge] OK {season}/{season + 1}: {times[0]:%Y-%m-%dT%H} -> {times[-1]:%Y-%m-%dT%H} "
          f"({len(times)} steps) -> {out}")
    return out


def remove_unused_raw(out_root, seasons, done):
    """Delete raw year Y once season Y-1 (uses Jan/Feb) and season Y (uses Dec) are both written."""
    for f in sorted((Path(out_root) / "raw").glob("era5_pl850_*.nc")):
        year = int(f.stem.split("_")[-1])
        users = [s for s in (year - 1, year) if s in seasons]
        if all(s in done for s in users):
            f.unlink()
            print(f"[cleanup] removed {f.name}")


# MAIN

def main():
    parser = argparse.ArgumentParser(description="Download ERA5 850 hPa u, v, vo and build DJF season files.")
    parser.add_argument("--out-root", type=Path, default=DEFAULT_OUT_ROOT,
                        help="folder for the season folders (use an external drive for the full range)")
    parser.add_argument("--first-season", type=int, default=1940, help="first season (December year)")
    parser.add_argument("--last-season", type=int, default=2025, help="last season (December year)")
    parser.add_argument("--days", type=int, nargs="+", default=list(range(1, 32)),
                        help="days of the month to download (use a few for testing)")
    parser.add_argument("--workers", type=int, default=4, help="number of CDS requests at the same time")
    parser.add_argument("--threads", type=int, default=4, help="cdo -P threads for regridding")
    parser.add_argument("--regrid", action=argparse.BooleanOptionalAction, default=True,
                        help="regrid every season after merging")
    parser.add_argument("--keep-raw", action="store_true", help="do not delete the yearly raw files")
    args = parser.parse_args()

    out_root = args.out_root
    seasons = list(range(args.first_season, args.last_season + 1))
    days = sorted(set(args.days))
    out_root.mkdir(parents=True, exist_ok=True)
    print(f"Seasons {seasons[0]}/{seasons[0] + 1} -> {seasons[-1]}/{seasons[-1] + 1}, "
          f"days {days[0]}-{days[-1]} ({len(days)}), output in {out_root}")

    # 1. download
    year_months = needed_months(out_root, seasons)
    failed_years = download_all(year_months, days, out_root, args.workers)

    # 2. merge (+ 3. regrid)
    done, failed = set(), []
    for s in seasons:
        if s in failed_years or s + 1 in failed_years:
            failed.append(s)
            continue
        try:
            merge_season(s, days, out_root)
            done.add(s)
            if not args.keep_raw:
                remove_unused_raw(out_root, seasons, done)
            if args.regrid:
                regrid_season(out_root / f"{s}{s + 1}_DJF", threads=args.threads)
        except Exception as e:
            print(f"[error] season {s}/{s + 1}: {e}")
            failed.append(s)

    n_files = sum(len(list((out_root / f"{s}{s + 1}_DJF").glob("*.nc"))) for s in done)
    print(f"[done] {len(done)} seasons, {n_files} files in {out_root} "
          f"(expected {len(done) * (len(RESOLUTIONS) if args.regrid else 1)})")
    if failed:
        print(f"[done] FAILED seasons: {failed}  (failed downloads: {failed_years}) - rerun to retry")


if __name__ == "__main__":
    main()
