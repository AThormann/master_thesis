#File containing helper functions used in both processing and validation of tracks. These functions are not specific to the tracking algorithm, but are used to process the output of the tracking algorithm and to compare tracks with reference tracks.
import numpy as np
import pandas as pd
from configs import EARTH_RADIUS_KM


# Convert longitudes between 0..360 and -180..180 conventions.
def wrap_lon(lon, to="-180"):
    """Convert longitudes between 0..360 and -180..180 conventions."""
    if to == "-180":
        return ((lon + 180) % 360) - 180
    return lon % 360

# Calculate the great-circle distance between two points on the Earth given their latitude and longitude in degrees. The result is in kilometers.   
def haversine_km(lat1, lon1, lat2, lon2):
    """Great-circle distance in km (works on numpy arrays)."""
    lat1, lon1, lat2, lon2 = map(np.radians, (lat1, lon1, lat2, lon2))
    a = np.sin((lat2 - lat1) / 2) ** 2 + np.cos(lat1) * np.cos(lat2) * np.sin((lon2 - lon1) / 2) ** 2
    return 2 * EARTH_RADIUS_KM * np.arcsin(np.sqrt(a))


# Determine if a point is inside af lat/lon box
def _in_box(lat, lon, box):
    lon_min, lat_min, lon_max, lat_max = box
    return (lon >= lon_min) & (lon <= lon_max) & (lat >= lat_min) & (lat <= lat_max)



def _to_datetime(times):
    """Track times are stored as ms offsets; also accept datetime64 arrays."""
    times = np.asarray(times)
    if np.issubdtype(times.dtype, np.datetime64):
        return pd.DatetimeIndex(times)
    return pd.to_datetime(times, unit="ms")

# convert a track object to a pandas dataframe with one row per track point, including track_id, time, lat, lon and all track variables. This is useful for post-processing and validation of tracks.
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



# summary of tracks with one row per track, including:
#  number of points, start/end time, genesis/lysis lat/lon, path length, number of points in region, fraction of points at edge, duration, separation between genesis and lysis, and peak value of the primary variable.
#  This is used in post-processing and validation of tracks.
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

