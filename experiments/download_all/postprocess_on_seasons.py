"""
Filter the raw tracks of every DJF season and resolution with POSTPROCESS_CONFIG from algorithms/configs.py.

Input:  <track_dir>/<S><S+1>_DJF_<res>_vo_raw.trackjson   (from track_on_seasons.py)
Output: <track_dir>/<S><S+1>_DJF_<res>_vo_filtered.trackjson
        (--track-dir, default <out_root>/trackjson)

Steps per file (same as experiments/validation_tracks/coarse_validation.ipynb):
1. load the raw tracks
2. time range of the season = first and last time step of the season's .nc file
   (replaces the hard-coded TIME_RANGE, handles leap years and test runs with --days)
3. clip the global tracks to DOMAIN and that time range (validation.clip_to_domain)
4. processing.postprocess_tracks with POSTPROCESS_CONFIG, time_range set per season
5. print the filter report and write the filtered tracks

Examples:
    python postprocess_on_seasons.py --out-root /Volumes/<drive>/ERA5_pressure_levels
    python postprocess_on_seasons.py --out-root /Volumes/<drive>/ERA5_pressure_levels --seasons 1998 --resolutions 3.0 --overwrite
"""

import pandas as pd
import pystormtracker as pst
import xarray as xr

from track_on_seasons import common_args, get_track_dir, season_dirs, track_file
import processing as pro
import helpers
import validation as val
from configs import DOMAIN, POSTPROCESS_CONFIG, track_variable

from regrid_era5 import RESOLUTIONS, res_tag, season_file


def season_time_range(season_dir, res):
    """(first, last) time step of the season's data file, e.g. ('1998-12-01 00:00', '1999-02-28 18:00')."""
    with xr.open_dataset(season_file(season_dir, res)) as ds:
        times = pd.DatetimeIndex(ds["valid_time"].values)
    return (f"{times[0]:%Y-%m-%d %H:%M}", f"{times[-1]:%Y-%m-%d %H:%M}")


def postprocess_season(season_dir, res, track_dir, overwrite=False):
    """Clip and filter the raw tracks of one season and resolution, write the filtered .trackjson."""
    raw = track_file(track_dir, season_dir, res, "raw")
    out = track_file(track_dir, season_dir, res, "filtered")
    if out.exists() and not overwrite:
        print(f"[post] SKIP {season_dir.name} {res_tag(res)} (exists)")
        return None
    if not raw.exists():
        print(f"[post] MISSING {raw}")
        return None

    tracks = pst.load_tracks(raw)
    time_range = season_time_range(season_dir, res)

    df = helpers.tracks_to_points(tracks)
    df = val.clip_to_domain(df, DOMAIN, time_range)
    if df.empty:
        print(f"[post] WARNING {season_dir.name} {res_tag(res)}: no track points inside DOMAIN")
        return None
    clipped = val.dataframe_to_tracks(df, variable=track_variable, scale=1)  # already in s^-1

    config = {**POSTPROCESS_CONFIG, "time_range": time_range}
    filtered, summary, _ = pro.postprocess_tracks(clipped, **config)
    print(pro.filter_report(summary).to_string())

    out.parent.mkdir(parents=True, exist_ok=True)
    tmp = out.with_name(out.name + ".part")
    filtered.write(tmp, format="json")
    tmp.replace(out)
    print(f"[post] OK {season_dir.name} {res_tag(res)} ({time_range[0]} -> {time_range[1]}): "
          f"raw {len(tracks)} -> clipped {len(clipped)} -> kept {len(filtered)} -> {out}")
    return len(filtered)


def main():
    parser = common_args("Filter the raw tracks of every season and resolution with POSTPROCESS_CONFIG.")
    parser.add_argument("--overwrite", action="store_true",
                        help="redo existing filtered files (e.g. after changing POSTPROCESS_CONFIG)")
    args = parser.parse_args()
    unknown = [r for r in args.resolutions if r not in RESOLUTIONS]
    if unknown:
        parser.error(f"unknown resolutions {unknown}, choose from {list(RESOLUTIONS)}")

    dirs = season_dirs(args.out_root, args.seasons)
    track_dir = get_track_dir(args)
    shown = {k: v for k, v in POSTPROCESS_CONFIG.items() if k != "time_range"}
    print(f"Post-processing {len(dirs)} seasons x {len(args.resolutions)} resolutions with {shown}")

    failed = []
    for d in dirs:
        for res in args.resolutions:
            try:
                postprocess_season(d, res, track_dir, args.overwrite)
            except Exception as e:
                print(f"[error] {d.name} {res_tag(res)}: {e}")
                failed.append(f"{d.name}_{res_tag(res)}")

    print(f"[done] filtered tracks in {track_dir}")
    if failed:
        print(f"[done] FAILED: {failed}")


if __name__ == "__main__":
    main()
