import os
import numpy as np
import xarray as xr
import pandas as pd
import pystormtracker as pst
from pystormtracker.hodges.rsplice import filter_rsplice
from pystormtracker.models.geo import DEG_TO_RAD, geod_dist

EARTH_RADIUS_KM = 6371.0
REGION = (-80, 20, 30, 80)   # smaller region for tracking, limited by the domain, 10 deg smaller
DOMAIN = (-90, 30, 20, 90)   # tracking domain (region + buffer), limited by the downloaded data
#format: (lon_min, lon_max, lat_min, lat_max)




# PREPARE FIELD FOR TRACKING

def prepare_tracking_field(ds, variable="vo", level=850, domain=None):
    """
    Return a clean (time, latitude, longitude) Dataset with only the tracking variable.

    - selects the pressure level (if there is one)
    - converts longitudes to -180..180 and sorts them
    - cuts out the tracking domain (lon_min, lon_max, lat_min, lat_max)
    """
    if isinstance(ds, (str, os.PathLike)):
        ds = xr.open_dataset(ds)
    da = ds[variable]
    if "pressure_level" in da.dims:
        da = da.sel(pressure_level=level, method="nearest")
    da = da.drop_vars([c for c in ("pressure_level", "expver", "number") if c in da.coords])

    lon = da["longitude"]
    if float(lon.max()) > 180:
        da = da.assign_coords(longitude=((lon + 180) % 360) - 180)
    da = da.sortby("longitude")

    if domain is not None:
        lon_min, lon_max, lat_min, lat_max = domain
        lat_slice = slice(lat_max, lat_min) if da["latitude"][0] > da["latitude"][-1] else slice(lat_min, lat_max)
        da = da.sel(longitude=slice(lon_min, lon_max), latitude=lat_slice)
    return da.to_dataset(name=variable)



# CONDUCT TRACKING

def hodges_tracker(input_file_path='name', output_file_name='output', track_variable='vo', min_track_points=3, taper_points=10, lmin=5,lmax=42,min_object_grid_points=18,
                   object_threshold=None, exclude_boundary_extrema=False):
    """
    Run the Hodges tracker on a file path or an already prepared xarray Dataset.

    Args:
        input_file_path (str or xarray.Dataset): Path to the input file or an xarray Dataset.
        output_file_name (str): Name of the output file to save the tracks.
        track_variable (str): Variable to track ('vo' for vorticity, 'msl' for mean sea level pressure).
        min_track_points (int): Minimum number of points required for a valid track.
        taper_points (int): Number of points to taper at the ends of tracks.
        lmin (float): min wavenumber bandpass filter
        lmax (float): Maximum wavenumber bandpass filter.
        dmax (float): Maximum distance between consecutive points in a track. (dmax (degrees per time step) is only used when use_dmax_zones=False; otherwise pystormtracker uses its default latitude zones (6.5 deg in the extratropics).)
        min_object_grid_points (int): Minimum number of grid points required for an object to be considered.
        object_threshold (float or None): Threshold value for object detection. If None, automatic detection is used.
        use_dmax_zones (bool): Whether to use default latitude zones for dmax. If False, dmax is used directly.
        exclude_boundary_extrema (bool): Whether to exclude extrema at the boundaries of the domain.

    
    """
    if isinstance(input_file_path, xr.Dataset):
        data_file = input_file_path
    else:
        data_file = xr.open_dataset(input_file_path)

    
    tracker = pst.HodgesTracker(min_track_points =min_track_points, taper_points=taper_points, lmin=lmin, lmax=lmax, min_object_grid_points=min_object_grid_points,
                                exclude_boundary_extrema=exclude_boundary_extrema)

    if track_variable == 'vo':
        detect = 'max'
    elif track_variable == 'msl':
        detect = 'min'
    else:
        detect = 'auto'

    tracks = tracker.track(
        data_file,
        track_variable,
        detection_mode= detect,
        object_threshold=object_threshold,
    )

    if tracks is None:
        return False
    elif len(tracks) == 0:
        return False
    else:
        
        if output_file_name:
        
            MODULE_DIR = os.path.dirname(os.path.abspath(__file__))
            DEFAULT_OUTPUT_DIR = os.path.join(MODULE_DIR, "..", "trackjson") 
            tracks.write(os.path.join(DEFAULT_OUTPUT_DIR, f"{output_file_name}.trackjson"))

        return tracks




#POST-PROCESSING of tracks:

#Helper functions for filtering and summarising tracks:

def haversine_km(lat1, lon1, lat2, lon2):
    """Great-circle distance in km (works on numpy arrays)."""
    lat1, lon1, lat2, lon2 = map(np.radians, (lat1, lon1, lat2, lon2))
    a = np.sin((lat2 - lat1) / 2) ** 2 + np.cos(lat1) * np.cos(lat2) * np.sin((lon2 - lon1) / 2) ** 2
    return 2 * EARTH_RADIUS_KM * np.arcsin(np.sqrt(a))



def _in_box(lat, lon, box):
    lon_min, lat_min, lon_max, lat_max = box
    return (lon >= lon_min) & (lon <= lon_max) & (lat >= lat_min) & (lat <= lat_max)


def tracks_to_points(tracks):
    """One row per track point: track_id, time, lat, lon and all track variables."""
    cols = {
        "track_id": np.repeat(np.asarray(tracks.ids), np.diff(tracks.offsets)),
        "time": _to_datetime(tracks.times),
        "lat": np.asarray(tracks.lats),
        "lon": np.asarray(tracks.lons),
    }
    for name, values in tracks.variables.items():
        cols[name] = np.asarray(values)
    return pd.DataFrame(cols)


#Used  in tracks_to_points and summarise_tracks
def _to_datetime(times):
    """Track times are stored as ms offsets; also accept datetime64 arrays."""
    times = np.asarray(times)
    if np.issubdtype(times.dtype, np.datetime64):
        return pd.DatetimeIndex(times)
    return pd.to_datetime(times, unit="ms")



# Summary used in filter
def summarise_tracks(tracks, region, domain, time_range=None, edge_buffer_deg=2.5):
    """
    One row per track with the quantities used by the post-processing filters.

    region : (lon_min, lon_max, lat_min, lat_max) that we are interested in
    domain : (lon_min, lon_max, lat_min, lat_max) of the field that was tracked
    time_range : (first_time, last_time) of the tracked data, used to flag
                 tracks that are cut off by the start/end of the period
    edge_buffer_deg : points closer than this to the domain edge count as "at the edge"
    """
    
    pts = tracks_to_points(tracks)
    var = tracks.primary_variable

    # distance travelled between consecutive points of the same track
    same = pts["track_id"].eq(pts["track_id"].shift())
    step = haversine_km(pts["lat"].shift(), pts["lon"].shift(), pts["lat"], pts["lon"])
    pts["step_km"] = np.where(same, step, 0.0)

    lon_min, lon_max, lat_min, lat_max = domain
    inner = (lon_min + edge_buffer_deg, lon_max - edge_buffer_deg,
             lat_min + edge_buffer_deg, lat_max - edge_buffer_deg)
    pts["in_region"] = _in_box(pts["lat"], pts["lon"], region)
    pts["at_edge"] = ~_in_box(pts["lat"], pts["lon"], inner)

    g = pts.groupby("track_id", sort=False)
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

    # intensity: vo -> maximum, msl -> minimum
    idx = g[var].idxmin() if tracks.mode == "min" else g[var].idxmax()
    s[f"{var}_peak"] = pts.loc[idx.values, var].values
    s["time_peak"] = pts.loc[idx.values, "time"].values
    s["lat_peak"] = pts.loc[idx.values, "lat"].values
    s["lon_peak"] = pts.loc[idx.values, "lon"].values
    s["peak_in_region"] = _in_box(s["lat_peak"], s["lon_peak"], region)

    if time_range is not None:
        t0, t1 = pd.Timestamp(time_range[0]), pd.Timestamp(time_range[1])
        s["starts_at_first_frame"] = s["start_time"] <= t0
        s["ends_at_last_frame"] = s["end_time"] >= t1
    return s.reset_index(), pts


def rsplice_summary(tracks, min_points=8, max_points=None, distance_degrees=9.0,
                    distance_mode="endpoint", cadence=np.timedelta64(6, "h")):
    """
    TRACK's RSPLICE filter (pystormtracker.hodges.rsplice.filter_rsplice) as a per-track table.

    Columns:
      rsplice_points           observed points + missing frames inside the track (gaps / cadence - 1)
      rsplice_displacement_deg great-circle displacement in degrees, endpoint (genesis -> lysis)
                               or travel (sum of all steps) depending on distance_mode
      ok_rsplice               True if filter_rsplice keeps the track
    """
    cadence_ms = None if cadence is None else int(cadence / np.timedelta64(1, "ms"))
    rows = []
    for tr in tracks:
        n = len(tr)
        if cadence_ms is not None and n > 1:
            n += int(np.sum(np.diff(tr.times) // cadence_ms - 1))
        if len(tr) < 2:
            disp = 0.0
        elif distance_mode == "endpoint":
            disp = geod_dist(float(tr.lats[0]), float(tr.lons[0]), float(tr.lats[-1]), float(tr.lons[-1]))
        else:
            disp = sum(geod_dist(float(tr.lats[i]), float(tr.lons[i]), float(tr.lats[i + 1]), float(tr.lons[i + 1]))
                       for i in range(len(tr) - 1))
        rows.append((tr.track_id, n, disp / DEG_TO_RAD))
    out = pd.DataFrame(rows, columns=["track_id", "rsplice_points", "rsplice_displacement_deg"])

    # the keep/drop decision itself comes from pystormtracker, so it matches TRACK exactly
    kept = filter_rsplice(tracks, min_points=min_points, max_points=max_points,
                          distance_degrees=distance_degrees, distance_mode=distance_mode,
                          expected_cadence=cadence)
    
    out["ok_rsplice"] = out["track_id"].isin(np.asarray(kept.ids))
    return out



def postprocess_tracks(tracks, region, domain, time_range=None,
                       region_criterion="passes", min_points_in_region=1,
                       rsplice_min_points=8, rsplice_max_points=None,
                       rsplice_distance_degrees=9.0, rsplice_distance_mode="endpoint",
                       rsplice_cadence=np.timedelta64(6, "h"),
                       min_peak=None, max_frac_at_edge=0.5,
                       edge_buffer_deg=2.5, drop_time_truncated=False):
    """
    Filter raw Hodges tracks down to the storms we want in the region.

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
    summary, points = summarise_tracks(tracks, region, domain, time_range, edge_buffer_deg)
    rs = rsplice_summary(tracks, min_points=rsplice_min_points, max_points=rsplice_max_points,
                         distance_degrees=rsplice_distance_degrees, distance_mode=rsplice_distance_mode,
                         cadence=rsplice_cadence)
    summary = summary.merge(rs, on="track_id", how="left", validate="one_to_one")
    var = tracks.primary_variable

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
    elif tracks.mode == "min":
        summary["ok_peak"] = summary[f"{var}_peak"] <= min_peak
    else:
        summary["ok_peak"] = summary[f"{var}_peak"] >= min_peak
    if drop_time_truncated and time_range is not None:
        summary["ok_time"] = ~(summary["starts_at_first_frame"] | summary["ends_at_last_frame"])
    else:
        summary["ok_time"] = True

    ok_cols = [c for c in summary.columns if c.startswith("ok_")]
    summary["keep"] = summary[ok_cols].all(axis=1)

    # summary rows are in the same order as the tracks, so the mask lines up
    filtered = tracks.filter(summary["keep"].to_numpy())

    return filtered, summary, points

def filter_report(summary):
    """How many tracks fail each criterion (a track can fail several)."""
    ok_cols = [c for c in summary.columns if c.startswith("ok_")]
    report = pd.DataFrame({
        "removed_by_this_criterion": [(~summary[c]).sum() for c in ok_cols],
        "removed_only_by_this_criterion": [
            ((~summary[c]) & summary[[o for o in ok_cols if o != c]].all(axis=1)).sum() for c in ok_cols
        ],
    }, index=[c[3:] for c in ok_cols])
    print(f"raw tracks: {len(summary)}, kept: {summary['keep'].sum()}")
    return report
