"""
Filter the raw tracks of every season (DJF, MAM, JJA, SON) and resolution with POSTPROCESS_CONFIG
from algorithms/configs.py.

Input:  <track_dir>/<season>/<season>_<res>_vo_raw.trackjson   (from track_on_seasons.py)
Output: <track_dir>/<season>/<season>_<res>_vo_filtered.trackjson
        (--track-dir, default <out_root>/trackjson)

Tracks are kept at full length (global tracking, no clipping). A track is kept if its
maximum vorticity lies in REGION and it has at least 4 points (24 h) in REGION
(region_criterion=("passes", "peak") in POSTPROCESS_CONFIG); the edge filter is off.
A track belongs to the season if its peak vorticity lies in the season itself (ok_season): the
season files include a buffer of SEASON_BUFFER_DAYS on both sides, and tracks peaking in the buffer
belong to the previous or next season (they are kept there).

Steps per file:
1. load the raw tracks
2. time range of the file = first and last time step of the season's .nc file (incl. buffer, used for
   the truncation flags); season range = the season without buffer (seasons.core_range of the folder name)
3. (disabled, kept as commented-out code) clip the global tracks to DOMAIN and that time range
   (validation.clip_to_domain), keeping every variable of the tracks (vo, raw_value, object_*)
4. processing.postprocess_tracks with POSTPROCESS_CONFIG, time_range and season_range set per season
5. print the filter report and write the filtered tracks

Examples:
    python postprocess_on_seasons.py --out-root /Volumes/<drive>/ERA5_pressure_levels
    python postprocess_on_seasons.py --out-root /Volumes/<drive>/ERA5_pressure_levels --seasons 1998 --resolutions 3.0 --overwrite
"""

from pathlib import Path

import numpy as np
import pandas as pd
import pystormtracker as pst
import xarray as xr

from pystormtracker.models.tracks import Tracks, TracksMetadata

from track_on_seasons import common_args, get_track_dir, season_dirs, track_file
import processing as pro
import helpers
import validation as val
from configs import DOMAIN, POSTPROCESS_CONFIG, track_variable

from regrid_era5 import RESOLUTIONS, res_tag, season_file
from seasons import SEASON_BUFFER_DAYS, buffered_range, core_range, parse_season


def season_time_range(season_dir, res):
    """
    (first, last) time step of the season's data file, e.g. ('1998-12-01 00:00', '1999-02-28 18:00').
    Uses the file of this resolution, or any other .nc file of the season if that one is missing.
    """
    path = season_file(season_dir, res)
    if not path.exists():
        path = next(iter(sorted(Path(season_dir).glob("*.nc"))), path)
    with xr.open_dataset(path) as ds:
        times = pd.DatetimeIndex(ds["valid_time"].values)
    return (f"{times[0]:%Y-%m-%d %H:%M}", f"{times[-1]:%Y-%m-%d %H:%M}")


def points_to_tracks(df, metadata):
    """
    Like validation.dataframe_to_tracks, but keeps every variable column of the per-point
    DataFrame (vo, raw_value, object_*), with the units and mode of the original tracks.
    """
    variables = [c for c in metadata.units if c in df.columns]
    df = df.dropna(subset=["lat", "lon", "time"]).sort_values(["track_id", "time"], kind="stable")
    counts = df.groupby("track_id", sort=False).size()
    meta = TracksMetadata(primary_variable=metadata.primary_variable, mode=metadata.mode,
                          units={v: metadata.units[v] for v in variables})
    return Tracks(
        ids=counts.index.to_numpy(dtype=np.int64),
        offsets=np.concatenate([[0], np.cumsum(counts.to_numpy())]),
        times=pd.DatetimeIndex(df["time"]).as_unit("ms").asi8.astype(np.int64),
        lats=df["lat"].to_numpy(float),
        lons=df["lon"].to_numpy(float),
        variables={v: df[v].to_numpy(float) for v in variables},
        metadata=meta,
    )


def postprocess_season(season_dir, res, track_dir, overwrite=False):
    """Filter the raw (full-length) tracks of one season and resolution, write the filtered .trackjson."""
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
    year, stype = parse_season(season_dir)
    season_range = tuple(f"{t:%Y-%m-%d %H:%M}" for t in core_range(year, stype))
    buf_start, buf_end = buffered_range(year, stype)
    if pd.Timestamp(time_range[0]) > buf_start or pd.Timestamp(time_range[1]) < buf_end:
        print(f"[post] WARNING {season_dir.name}: data file covers {time_range[0]} -> {time_range[1]}, "
              f"not the {SEASON_BUFFER_DAYS}-day buffer {buf_start} -> {buf_end}; tracks crossing the season "
              f"boundary may get a wrong peak (rebuild the season with ERA5_download_seasons.py)")

    # clipping to DOMAIN disabled: full-length tracks are kept (selection by peak and >= 4 points in REGION)
    # df = helpers.tracks_to_points(tracks)
    # df = val.clip_to_domain(df, DOMAIN, time_range)
    # if df.empty:
    #     print(f"[post] WARNING {season_dir.name} {res_tag(res)}: no track points inside DOMAIN")
    #     return None
    # clipped = points_to_tracks(df, tracks.metadata)  # vo already in s^-1

    config = {**POSTPROCESS_CONFIG, "time_range": time_range, "season_range": season_range}
    filtered, summary, _ = pro.postprocess_tracks(tracks, **config)
    print(pro.filter_report(summary).to_string())

    out.parent.mkdir(parents=True, exist_ok=True)
    tmp = out.with_name(out.name + ".part")
    filtered.write(tmp, format="json")
    tmp.replace(out)
    print(f"[post] OK {season_dir.name} {res_tag(res)} (peak in {season_range[0]} -> {season_range[1]}): "
          f"raw {len(tracks)} -> kept {len(filtered)} -> {out}")
    return len(filtered)


def main():
    parser = common_args("Filter the raw tracks of every season and resolution with POSTPROCESS_CONFIG.")
    parser.add_argument("--overwrite", action="store_true",
                        help="redo existing filtered files (e.g. after changing POSTPROCESS_CONFIG)")
    args = parser.parse_args()
    unknown = [r for r in args.resolutions if r not in RESOLUTIONS]
    if unknown:
        parser.error(f"unknown resolutions {unknown}, choose from {list(RESOLUTIONS)}")

    dirs = season_dirs(args.out_root, args.seasons, args.season_types)
    track_dir = get_track_dir(args)
    shown = {k: v for k, v in POSTPROCESS_CONFIG.items() if k not in ("time_range", "season_range")}
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
