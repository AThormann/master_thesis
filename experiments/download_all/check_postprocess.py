"""
Check the post-processing of one season and resolution (full-length tracks, region criteria
"passes" + "peak", edge filter off).

    python check_postprocess.py <raw.trackjson> [--filtered <filtered.trackjson>]

Without --filtered, the raw tracks are filtered here with POSTPROCESS_CONFIG (time_range=None).
With --filtered, the file written by postprocess_on_seasons.py is checked instead.

Checks:
- every kept track has its point of maximum vo inside REGION and at least min_points_in_region points in REGION
- the edge filter removes no tracks
- some kept tracks have points outside DOMAIN (tracks are not clipped)
- every kept track has the same number of points as the same track_id in the raw file
- all track variables of the raw file (vo, raw_value, object_*) are still present
- region_criterion="passes" (one string) gives the same tracks as n_in_region >= min_points
  combined with the other filters
"""

import argparse
import sys
from pathlib import Path

import numpy as np

ALGORITHMS_DIR = Path(__file__).resolve().parents[2] / "algorithms"
sys.path.insert(0, str(ALGORITHMS_DIR))
import pystormtracker as pst
import processing as pro
import helpers
from configs import DOMAIN, POSTPROCESS_CONFIG, REGION


def in_box(lat, lon, box):
    lon_min, lon_max, lat_min, lat_max = box
    return (lon >= lon_min) & (lon <= lon_max) & (lat >= lat_min) & (lat <= lat_max)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("raw")
    parser.add_argument("--filtered", default=None)
    args = parser.parse_args()

    raw = pst.load_tracks(args.raw)
    config = dict(POSTPROCESS_CONFIG)
    filtered, summary, _ = pro.postprocess_tracks(raw, **config)
    print(pro.filter_report(summary).to_string())
    if args.filtered:
        filtered = pst.load_tracks(args.filtered)
        print(f"checking {args.filtered}")
    n_min = config["min_points_in_region"]

    # region criteria on the kept tracks
    pts = helpers.tracks_to_points(filtered)
    var = filtered.primary_variable
    pts["in_region"] = in_box(pts["lat"], pts["lon"], REGION)
    g = pts.groupby("track_id")
    idx = g[var].idxmin() if filtered.mode == "min" else g[var].idxmax()
    peak = pts.loc[idx.values]
    assert in_box(peak["lat"], peak["lon"], REGION).all(), "a kept track peaks outside REGION"
    assert (g["in_region"].sum() >= n_min).all(), f"a kept track has < {n_min} points in REGION"

    # edge filter off
    assert summary["ok_edge"].all(), "the edge filter removed tracks"

    # full length: some points outside DOMAIN, same number of points as in the raw file
    n_outside = int((~in_box(pts["lat"], pts["lon"], DOMAIN)).sum())
    assert n_outside > 0, "no kept points outside DOMAIN: tracks look clipped"
    raw_counts = helpers.tracks_to_points(raw).groupby("track_id").size()
    kept_counts = g.size()
    assert (kept_counts == raw_counts.loc[kept_counts.index]).all(), "kept tracks differ in length from the raw ones"

    # variables
    missing = set(raw.variables) - set(filtered.variables)
    assert not missing, f"variables lost: {missing}"

    # single string criterion
    _, s_passes, _ = pro.postprocess_tracks(raw, **{**config, "region_criterion": "passes"})
    others = [c for c in s_passes.columns if c.startswith("ok_") and not c.startswith("ok_region")]
    expected = (s_passes["n_in_region"] >= n_min) & s_passes[others].all(axis=1)
    assert "ok_region_passes" in s_passes and "ok_region" not in s_passes
    assert np.array_equal(s_passes["keep"].to_numpy(), expected.to_numpy()), "'passes' alone gives other tracks"

    print(f"OK: {len(raw)} raw -> {len(filtered)} kept ('passes' alone: {int(s_passes['keep'].sum())}), "
          f"{n_outside} kept points outside DOMAIN, variables {sorted(filtered.variables)}")


if __name__ == "__main__":
    main()
