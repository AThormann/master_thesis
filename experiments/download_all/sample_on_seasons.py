"""
Sample ERA5 fields along every point of the filtered tracks of every DJF season and resolution.

Input:  <track_dir>/<S><S+1>_DJF/<S><S+1>_DJF_<res>_vo_filtered.trackjson   (from postprocess_on_seasons.py)
        <out_root>/<S><S+1>_DJF/<S><S+1>_DJF_sl_<sample_res>.nc            single levels (u10, v10, msl)
        <out_root>/<S><S+1>_DJF/<S><S+1>_DJF_0p25.nc                       850 hPa (u, v, vo)
Output: <track_dir>/<S><S+1>_DJF/<S><S+1>_DJF_<res>_vo_points_sl<sample_res>.csv   one row per track point
        <track_dir>/<S><S+1>_DJF/<S><S+1>_DJF_<res>_vo_summary_sl<sample_res>.csv  one row per track

What is sampled is set by SAMPLE_CONFIG in algorithms/configs.py (statistics_on_tracks.sample_along_tracks).
--sample-res other than 0.25 regrids the single-level file first (cdo remapcon, regrid_era5.py), so the
same tracks can be sampled from coarser fields.

Examples:
    python sample_on_seasons.py --out-root /Volumes/<drive>/ERA5_pressure_levels --track-dir ../../trackjson/seasons
    python sample_on_seasons.py --out-root /tmp/era5_test --track-dir ../../trackjson/seasons --seasons 1998 --resolutions 1.0 --sample-res 1.0
"""

import time
import warnings

import pystormtracker as pst

from track_on_seasons import common_args, get_track_dir, season_dirs, track_file
import statistics_on_tracks as stats
from configs import SAMPLE_CONFIG, SL_AREA

from regrid_era5 import NATIVE_RES, RESOLUTIONS, regrid_season, res_tag, season_file


def sample_season(season_dir, res, track_dir, sample_res=NATIVE_RES, overwrite=False):
    """Sample the fields along the filtered tracks of one season and resolution and write the CSVs."""
    tag = f"sl{res_tag(sample_res)}"
    tracks_path = track_file(track_dir, season_dir, res, "filtered")
    points_path = track_file(track_dir, season_dir, res, f"points_{tag}", "csv")
    summary_path = track_file(track_dir, season_dir, res, f"summary_{tag}", "csv")
    if points_path.exists() and not overwrite:
        print(f"[sample] SKIP {season_dir.name} {res_tag(res)} {tag} (exists)")
        return None
    if not tracks_path.exists():
        print(f"[sample] MISSING {tracks_path}")
        return None

    sl_path = season_file(season_dir, sample_res, "sl")
    if not sl_path.exists() and sample_res != NATIVE_RES:
        regrid_season(season_dir, resolutions=[sample_res], kind="sl", area=SL_AREA) #get a parameter for type of remap
    pl_path = season_file(season_dir, NATIVE_RES, "pl")
    for path in (sl_path, pl_path):
        if not path.exists():
            raise FileNotFoundError(f"Missing field file {path}")

    t0 = time.time()
    tracks = pst.load_tracks(tracks_path)
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        points = stats.sample_along_tracks(tracks, sl_path, pl_path, config=SAMPLE_CONFIG, area=SL_AREA)
    for w in caught:
        if "binary incompatibility" not in str(w.message):   # harmless numpy/netCDF4 import notice
            print(f"[sample] WARNING {season_dir.name} {res_tag(res)}: {w.message}")

    # the table must hold exactly the points of the filtered tracks
    n_points = len(tracks.times)
    if len(points) != n_points or set(points["track_id"]) != set(int(i) for i in tracks.ids):
        raise ValueError(f"points table has {len(points)} rows / {points['track_id'].nunique()} tracks, "
                         f"the filtered tracks have {n_points} points / {len(tracks)} tracks")

    summary = stats.summarise_points(points)
    points_path.parent.mkdir(parents=True, exist_ok=True)
    points.to_csv(points_path, index=False)
    summary.to_csv(summary_path, index=False)
    n_out = int((~points["inside_area"]).sum())
    print(f"[sample] OK {season_dir.name} {res_tag(res)} {tag}: {len(tracks)} tracks, {len(points)} points, "
          f"{n_out} outside area ({(time.time() - t0) / 60:.1f} min) -> {points_path.name}")
    return len(points)


def main():
    parser = common_args("Sample ERA5 fields along the filtered tracks of every season and resolution.")
    parser.add_argument("--sample-res", type=float, default=SAMPLE_CONFIG.get("sample_res", NATIVE_RES),
                        help="grid of the sampled single-level field (default from SAMPLE_CONFIG, 0.25 = native)")
    parser.add_argument("--overwrite", action="store_true", help="redo existing CSV files")
    args = parser.parse_args()
    for r in [*args.resolutions, args.sample_res]:
        if r not in RESOLUTIONS:
            parser.error(f"unknown resolution {r}, choose from {list(RESOLUTIONS)}")

    dirs = season_dirs(args.out_root, args.seasons)
    track_dir = get_track_dir(args)
    print(f"Sampling {len(dirs)} seasons x {len(args.resolutions)} track resolutions, "
          f"fields on {res_tag(args.sample_res)}, {SAMPLE_CONFIG}")

    failed = []
    for d in dirs:
        for res in args.resolutions:
            try:
                sample_season(d, res, track_dir, args.sample_res, args.overwrite)
            except Exception as e:
                print(f"[error] {d.name} {res_tag(res)}: {e}")
                failed.append(f"{d.name}_{res_tag(res)}")

    print(f"[done] sampled tracks in {track_dir}")
    if failed:
        print(f"[done] FAILED: {failed}")


if __name__ == "__main__":
    main()
