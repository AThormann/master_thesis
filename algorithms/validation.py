
import os
import warnings

import numpy as np
import xarray as xr
import netCDF4 as nc
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import pandas as pd

import pystormtracker as pst
from pystormtracker.sample import sample_tracks
from pystormtracker.hodges.rsplice import filter_rsplice
from pystormtracker.models.geo import DEG_TO_RAD, geod_dist

MISSING = 1.0e25
COLS = ["step", "lon", "lat", "vorticity", "mslp_lon", "mslp_lat", "mslp"]

def wrap_lon(lon, to="-180"):
    """Convert longitudes between 0..360 and -180..180 conventions."""
    if to == "-180":
        return ((lon + 180) % 360) - 180
    return lon % 360


def haversine_km(lat1, lon1, lat2, lon2):
    """Great-circle distance in km (works on numpy arrays)."""
    lat1, lon1, lat2, lon2 = map(np.radians, (lat1, lon1, lat2, lon2))
    a = np.sin((lat2 - lat1) / 2) ** 2 + np.cos(lat1) * np.cos(lat2) * np.sin((lon2 - lon1) / 2) ** 2
    return 2 * EARTH_RADIUS_KM * np.arcsin(np.sqrt(a))

def _in_box(lat, lon, box):
    lon_min, lon_max, lat_min, lat_max = box
    return (lon >= lon_min) & (lon <= lon_max) & (lat >= lat_min) & (lat <= lat_max)



def load_tracks(path, start_time=None, end_time=None):
    """Parse a TRACK-format file into an xarray.Dataset with one `obs` dimension.

    `track_id` is a coordinate on `obs`. Missing values (1e25) become NaN.
    """
    rows, track_ids = [], []
    current_id = None

    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith("TRACK_ID"):
                current_id = int(line.split()[1])
            elif current_id is not None and "&" in line:
                vals = line.replace("&", " ").split()
                if len(vals) == len(COLS):
                    rows.append([float(v) for v in vals])
                    track_ids.append(current_id)

    data = np.array(rows)
    data[data >= MISSING] = np.nan

    ds = xr.Dataset(
        {name: ("obs", data[:, i]) for i, name in enumerate(COLS)},
        coords={"track_id": ("obs", np.array(track_ids))}
    )
    ds["step"] = ds["step"].astype(int)
    ds["lon"] = wrap_lon(ds.lon)
    times = pd.date_range(start=start_time, end=end_time, freq="6h")
    ds['time'] = times[ds['step']-1] # times are created from the step in the data set and step 1 is equal to the start time
    ds.pop("step") #delete the step column
    ds =  ds.to_pandas()
    return ds



def summarise_tracks(dataset, region, domain, time_range=None, edge_buffer_deg=2.5):
    """
    One row per track with the quantities used by the post-processing filters.

    region : (lon_min, lon_max, lat_min, lat_max) that we are interested in
    domain : (lon_min, lon_max, lat_min, lat_max) of the field that was tracked
    time_range : (first_time, last_time) of the tracked data, used to flag
                 tracks that are cut off by the start/end of the period
    edge_buffer_deg : points closer than this to the domain edge count as "at the edge"
    """
    
   
    #var = 'vorticity'


    # distance travelled between consecutive points of the same track
    same = dataset["track_id"].eq(dataset["track_id"].shift())
    step = haversine_km(dataset["lat"].shift(), dataset["lon"].shift(), dataset["lat"], dataset["lon"])
    dataset["step_km"] = np.where(same, step, 0.0)

    lon_min, lon_max, lat_min, lat_max = domain
    inner = (lon_min + edge_buffer_deg, lon_max - edge_buffer_deg,
             lat_min + edge_buffer_deg, lat_max - edge_buffer_deg)
    dataset["in_region"] = _in_box(dataset["lat"], dataset["lon"], region)
    dataset["at_edge"] = ~_in_box(dataset["lat"], dataset["lon"], inner)

    g = dataset.groupby("track_id", sort=False)
    s = pd.DataFrame({
        "n_points": g.size(),
        "start_time": g["time"].first(),
        "end_time": g["time"].last(),
        "lat_genesis": g["lat"].first(),
        "lon_genesis": g["lon"].first(),
        "lat_lysis": g["lat"].last(),
        "lon_lysis": g["lon"].last(),
        "path_km": g["step_km"].sum(),
        "n_in_region": g["in_region"].sum(),
        "frac_at_edge": g["at_edge"].mean(),
    })
    s["duration_h"] = (s["end_time"] - s["start_time"]) / pd.Timedelta(hours=1)
    s["separation_km"] = haversine_km(s["lat_genesis"], s["lon_genesis"], s["lat_lysis"], s["lon_lysis"])
    s["genesis_in_region"] = _in_box(s["lat_genesis"], s["lon_genesis"], region)

    # intensity: vo -> maximum
    idx = g['vorticity'].idxmax()
    s[f"{'vorticity'}_peak"] = dataset.loc[idx.values, 'vorticity'].values
    s["time_peak"] = dataset.loc[idx.values, "time"].values
    s["lat_peak"] = dataset.loc[idx.values, "lat"].values
    s["lon_peak"] = dataset.loc[idx.values, "lon"].values
    s["peak_in_region"] = _in_box(s["lat_peak"], s["lon_peak"], region)

    if time_range is not None:
        t0, t1 = pd.Timestamp(time_range[0]), pd.Timestamp(time_range[1])
        s["starts_at_first_frame"] = s["start_time"] <= t0
        s["ends_at_last_frame"] = s["end_time"] >= t1
    return s.reset_index(), dataset


def postprocess_validation_tracks(dataset, region, domain, time_range=None,
                       region_criterion="passes", min_points_in_region=1,
                       rsplice_min_points=8, rsplice_max_points=None,
                       rsplice_distance_degrees=10.0, rsplice_distance_mode="endpoint",
                       rsplice_cadence=np.timedelta64(6, "h"),
                       min_peak=None, max_frac_at_edge=0.5,
                       edge_buffer_deg=2.5, drop_time_truncated=False):
    """
    Filter validation dataset tracks down to the storms we want in the region.

    Every criterion becomes a boolean column in `summary` (True = passes), and
    `keep` is the AND of all of them, so you can see why each track was removed.

    region_criterion : "passes"  -> at least `min_points_in_region` points in the region
                       "genesis" -> the first point is in the region
                       "peak"    -> the point of maximum intensity is in the region
    rsplice_*         : lifetime and displacement filter of TRACK (filter_rsplice):
                        keep tracks with min_points <= points (incl. missing frames) <= max_points
                        and displacement >= distance_degrees ("endpoint" or "travel")
    rsplice_cadence   : time step of the data, used to count missing frames (ERA5 here: 6 h)
    min_peak          : minimum peak intensity of the tracked variable (vo: s-1, msl: Pa; msl uses <=)
    max_frac_at_edge  : maximum fraction of points within `edge_buffer_deg` of the domain edge
    drop_time_truncated : also remove tracks touching the first/last time step

    Returns (filtered Tracks, summary DataFrame for all tracks, per-point DataFrame for all tracks)
    """
    summary, points = summarise_tracks(dataset, region, domain, time_range, edge_buffer_deg)


    # rs = rsplice_summary(tracks, min_points=rsplice_min_points, max_points=rsplice_max_points,
    #                      distance_degrees=rsplice_distance_degrees, distance_mode=rsplice_distance_mode,
    #                      cadence=rsplice_cadence)
    # summary = summary.merge(rs, on="track_id", how="left", validate="one_to_one")
    

    if region_criterion == "passes":
        summary["ok_region"] = summary["n_in_region"] >= min_points_in_region
    elif region_criterion == "genesis":
        summary["ok_region"] = summary["genesis_in_region"]
    elif region_criterion == "peak":
        summary["ok_region"] = summary["peak_in_region"]
    else:
        raise ValueError("region_criterion must be 'passes', 'genesis' or 'peak'")

    summary["ok_edge"] = summary["frac_at_edge"] <= max_frac_at_edge
    if min_peak is None:
        summary["ok_peak"] = True
    else:
        summary["ok_peak"] = summary[f"{'vorticity'}_peak"] >= min_peak
    if drop_time_truncated and time_range is not None:
        summary["ok_time"] = ~(summary["starts_at_first_frame"] | summary["ends_at_last_frame"])
    else:
        summary["ok_time"] = True

    ok_cols = [c for c in summary.columns if c.startswith("ok_")]
    summary["keep"] = summary[ok_cols].all(axis=1)

    # summary rows are in the same order as the tracks, so the mask lines up
    filtered = dataset.filter(summary["keep"].to_numpy())

    return filtered, summary

