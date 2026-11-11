"""
Run the Hodges tracker on every DJF season and resolution and save the raw tracks.

Input:  <out_root>/<S><S+1>_DJF/<S><S+1>_DJF_<res>.nc   (from ERA5_download_seasons.py / regrid_era5.py)
Output: <track_dir>/<S><S+1>_DJF/<S><S+1>_DJF_<res>_vo_raw.trackjson   (--track-dir, default <out_root>/trackjson)

Tracking is done on the whole globe (the spectral filter needs the full sphere) with
TRACK_CONFIGS, track_variable and object_threshold from algorithms/configs.py.
Existing output files are skipped, so a run can be restarted.

Examples:
    python track_on_seasons.py --out-root /Volumes/<drive>/ERA5_pressure_levels
    python track_on_seasons.py --out-root /Volumes/<drive>/ERA5_pressure_levels --seasons 1998 --resolutions 3.0 1.875
"""

import argparse
import sys
import time
from pathlib import Path

ALGORITHMS_DIR = Path(__file__).resolve().parents[2] / "algorithms"
sys.path.insert(0, str(ALGORITHMS_DIR))
import processing as pro
from configs import TRACK_CONFIGS, track_variable, object_threshold

from regrid_era5 import RESOLUTIONS, res_tag, season_file


DEFAULT_OUT_ROOT = Path("/Users/asgerthormann/Documents/uni/master/master_code/master_thesis/reanalysis_data/ERA5/ERA5_pressure_levels")
LEVEL = 850


# shared with postprocess_on_seasons.py

def season_dirs(out_root, seasons=None):
    """Season folders <S><S+1>_DJF in out_root (all of them if seasons is None)."""
    out_root = Path(out_root)
    if seasons:
        return [out_root / f"{s}{s + 1}_DJF" for s in seasons]
    return sorted(d for d in out_root.glob("*_DJF") if d.is_dir())


def track_file(track_dir, season_dir, res, kind="raw", suffix="trackjson"):
    """
    <track_dir>/<S><S+1>_DJF/<S><S+1>_DJF_<res>_<var>_<kind>.<suffix>
    kind = 'raw' or 'filtered' (.trackjson), or e.g. 'points_sl0p25' / 'summary_sl0p25' (.csv)
    """
    season = Path(season_dir).name
    return Path(track_dir) / season / f"{season}_{res_tag(res)}_{track_variable}_{kind}.{suffix}"


def get_track_dir(args):
    """--track-dir if given, otherwise <out_root>/trackjson"""
    return args.track_dir or args.out_root / "trackjson"


def common_args(description):
    """CLI arguments used by both track_on_seasons.py and postprocess_on_seasons.py."""
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("--out-root", type=Path, default=DEFAULT_OUT_ROOT,
                        help="folder containing the <S><S+1>_DJF season folders")
    parser.add_argument("--seasons", type=int, nargs="+", help="season start years (default: all season folders)")
    parser.add_argument("--resolutions", type=float, nargs="+", default=list(RESOLUTIONS),
                        help="grid resolutions in degrees (default: all 10)")
    parser.add_argument("--track-dir", type=Path, default=None,
                        help="folder for the .trackjson files (default: <out-root>/trackjson)")
    return parser


# TRACKING

def track_season(season_dir, res, track_dir):
    """Track one season at one resolution and write the raw .trackjson. Returns the number of tracks."""
    out = track_file(track_dir, season_dir, res, "raw")
    if out.exists():
        print(f"[track] SKIP {season_dir.name} {res_tag(res)} (exists)")
        return None
    data_path = season_file(season_dir, res)
    if not data_path.exists():
        print(f"[track] MISSING {data_path}")
        return None

    t0 = time.time()
    ds = pro.prepare_tracking_field(data_path, variable=track_variable, level=LEVEL, domain=None)
    # output_file_name=None: hodges_tracker's own output folder is not used, the file is written here instead
    tracks = pro.hodges_tracker(ds, output_file_name=None, track_var=track_variable,
                                obj_threshold=object_threshold, **TRACK_CONFIGS)
    if tracks is False:
        print(f"[track] WARNING {season_dir.name} {res_tag(res)}: no tracks found")
        return 0

    out.parent.mkdir(parents=True, exist_ok=True)
    tmp = out.with_name(out.name + ".part")
    tracks.write(tmp, format="json")
    tmp.rename(out)
    print(f"[track] OK {season_dir.name} {res_tag(res)}: {len(tracks)} tracks "
          f"({(time.time() - t0) / 60:.1f} min) -> {out}")
    return len(tracks)


def main():
    parser = common_args("Track every ERA5 DJF season and resolution with the Hodges tracker.")
    args = parser.parse_args()
    unknown = [r for r in args.resolutions if r not in RESOLUTIONS]
    if unknown:
        parser.error(f"unknown resolutions {unknown}, choose from {list(RESOLUTIONS)}")

    dirs = season_dirs(args.out_root, args.seasons)
    track_dir = get_track_dir(args)
    print(f"Tracking {len(dirs)} seasons x {len(args.resolutions)} resolutions, "
          f"variable {track_variable} at {LEVEL} hPa, {TRACK_CONFIGS}, object_threshold={object_threshold}")

    failed = []
    for d in dirs:
        for res in args.resolutions:
            try:
                track_season(d, res, track_dir)
            except Exception as e:
                print(f"[error] {d.name} {res_tag(res)}: {e}")
                failed.append(f"{d.name}_{res_tag(res)}")

    print(f"[done] tracks in {track_dir}")
    if failed:
        print(f"[done] FAILED: {failed}")


if __name__ == "__main__":
    main()
