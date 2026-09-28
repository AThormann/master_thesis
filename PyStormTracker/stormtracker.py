import numpy as np
import xarray as xr
import netCDF4 as nc
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import pandas as pd

import pystormtracker as pst

def hodges_tracker(input_file_path='name', output_file_name='output', track_variable='vo', min_track_points=10, taper_points=10, lmin=5,lmax=42,dmax=9.0,min_object_grid_points=18):
    
    data_file = xr.open_dataset(input_file_path)
    tracker = pst.HodgesTracker(min_track_points =min_track_points, taper_points=taper_points, lmin=lmin, lmax=lmax, dmax=dmax, min_object_grid_points=min_object_grid_points)

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
        #object_threshold=2e-5
    )

    if tracks is None:
        return False
    elif len(tracks) == 0:
        return False
    else:
        #path = "PyStormTracker/json"
        tracks.write(f"{output_file_name}.trackjson")
        "Tracking done"
        return tracks




def plot_tracks(data, min_length=0, projection=ccrs.PlateCarree(), points=True, text=True):
    fig, ax = plt.subplots(figsize=(10, 8), subplot_kw={'projection': projection})
    ax.coastlines()
    ax.set_extent([np.min(data.lons), np.max(data.lons), np.min(data.lats), np.max(data.lats)], crs=projection)
    ax.gridlines(draw_labels=True)

    for d in data:
        if len(d) >= min_length:
            ax.plot(d.lons, d.lats, transform=projection)
            if points == True:
                ax.plot(d.lons, d.lats, 'go', transform=projection)  # Track points
            ax.plot(d.lons[0], d.lats[0], 'k', transform=projection)  # Start point
            if text == True:
                ax.text(d.lons[0], d.lats[0], str(pd.to_datetime(d.times[0], unit='ms')), transform=projection)  # Track start time
    ax.set_title(f"ETC tracks in NA between {pd.to_datetime(data[0].times[0], unit="ms")} and {pd.to_datetime(data[-1].times[-1], unit="ms")}")







#plotting method for tracks in json format
#projection

# --------------------------------------------------------------------------
# Track statistics for data analysis
# --------------------------------------------------------------------------
 
# (field in dataset, sampling method, output column name)
# "ws" is wind speed and is computed from wind_vars if it is not in the dataset.
DEFAULT_SAMPLES = (
    ("t2m", "min",  "t_min"),
    ("t2m", "max",  "t_max"),
    ("msl", "min",  "msl_min"),
    ("ws",  "max",  "ws_max"),
    ("ws",  "mean", "ws_mean"),
)
 
# column name in the output -> variable name stored on the tracks object
DEFAULT_OBJECT_COLUMNS = {
    "major_km":          "object_moment_major_axis_km",
    "minor_km":          "object_moment_minor_axis_km",
    "gridcell_area_km2": "object_gridcell_area_km2",
    "fitted_area_km2":   "object_moment_fitted_area_km2",
}
 
 
def _to_datetime(times):
    """Track times are stored as ms offsets; also accept datetime64 arrays."""
    times = np.asarray(times)
    if np.issubdtype(times.dtype, np.datetime64):
        return pd.DatetimeIndex(times)
    return pd.to_datetime(times, unit="ms")
 
 
def _check_time_overlap(track_times, ds, strict=True):
    """Make sure the track times are covered by the sampled dataset."""
    time_name = next((n for n in ("time", "valid_time") if n in ds.coords), None)
    if time_name is None:
        raise ValueError("Could not find a 'time' or 'valid_time' coordinate in the dataset.")
 
    ds_times = pd.DatetimeIndex(ds[time_name].values)
    inside = np.asarray((track_times >= ds_times.min()) & (track_times <= ds_times.max()))
    exact = np.asarray(track_times.isin(ds_times))
 
    msg = (f"tracks span {track_times.min()} -> {track_times.max()}, "
           f"dataset spans {ds_times.min()} -> {ds_times.max()}; "
           f"{inside.mean():.1%} of track points inside the dataset range, "
           f"{exact.mean():.1%} with an exact time match.")
 
    if inside.mean() < 1.0:
        if strict:
            raise ValueError("Track times not covered by dataset: " + msg)
        warnings.warn("Track times not fully covered by dataset: " + msg)
    elif exact.mean() < 1.0:
        warnings.warn("Some track times have no exact match in the dataset "
                      "(different time step?): " + msg)
 
 
def track_statistics(tracks, ds, samples=DEFAULT_SAMPLES, radius_km=500,
                     wind_vars=("u10", "v10"), level=None, strict_time=True,
                     object_columns=DEFAULT_OBJECT_COLUMNS, extra_agg=None):
    """
    Sample fields around each track point and summarise every track.
 
    Parameters
    ----------
    tracks : Tracks or str/Path
        Tracks object (from hodges_tracker) or path to a .trackjson file.
    ds : xarray.Dataset or str/Path
        Dataset with the fields to sample (e.g. t2m, msl, u10, v10).
    samples : iterable of (field, method, output_name)
        method is one of nearest, bilinear, mean, max, min.
        Use the field name "ws" for wind speed (computed from wind_vars).
    radius_km : float
        Radius used by the mean/max/min methods.
    wind_vars : (str, str)
        Names of the u and v components used to compute "ws" on the grid.
    level : float, optional
        Pressure level to select if the dataset has a 'pressure_level' dimension.
    strict_time : bool
        Raise if track times fall outside the dataset time range (else warn).
    object_columns : dict
        {output column: variable name on tracks} for the object geometry.
    extra_agg : dict, optional
        Extra per-track aggregations as {name: (column, func)}, e.g.
        {"major_km_max": ("major_km", "max")}.
 
    Returns
    -------
    points : pandas.DataFrame
        One row per track point.
    summary : pandas.DataFrame
        One row per track (extremes, lifetime, genesis/lysis position).
    """
    # --- inputs ------------------------------------------------------------
    if isinstance(tracks, (str, os.PathLike)):
        tracks = pst.load_tracks(str(tracks))
    if isinstance(ds, (str, os.PathLike)):
        ds = xr.open_dataset(ds)
    if tracks is None or len(tracks) == 0:
        raise ValueError("No tracks to sample.")
 
    samples = [tuple(s) for s in samples]
    outputs = [s[2] for s in samples]
    if len(set(outputs)) != len(outputs):
        raise ValueError(f"Output names must be unique, got {outputs}")
 
    if level is not None and "pressure_level" in ds.dims:
        ds = ds.sel(pressure_level=level, method="nearest")
        print("pressure level used:", float(ds["pressure_level"]))
 
    # --- wind speed on the grid (not from averaged u and v) -----------------
    if any(s[0] == "ws" for s in samples) and "ws" not in ds:
        u_name, v_name = wind_vars
        missing = [n for n in (u_name, v_name) if n not in ds]
        if missing:
            raise KeyError(f"Wind components {missing} not in dataset; needed for 'ws'.")
        ws = np.hypot(ds[u_name], ds[v_name])
        ws.attrs["units"] = "m s-1"
        ds = ds.assign(ws=ws)
 
    missing = sorted({s[0] for s in samples} - set(ds.data_vars))
    if missing:
        raise KeyError(f"Fields {missing} not in dataset (available: {list(ds.data_vars)}).")
 
    # --- time sanity check --------------------------------------------------
    _check_time_overlap(_to_datetime(tracks.times), ds, strict=strict_time)
 
    # --- sampling (each call returns a new Tracks object) --------------------
    trk = tracks
    for field, method, out in samples:
        trk = sample_tracks(trk, ds, variable_name=field, method=method,
                            radius_km=radius_km, output_variable_name=out)
 
    # --- per-point table ----------------------------------------------------
    V = trk.variables
    cols = {
        "track_id": np.repeat([tr.track_id for tr in trk], [len(tr) for tr in trk]),
        "time":     _to_datetime(trk.times),
        "lat":      np.asarray(trk.lats),
        "lon":      np.asarray(trk.lons),
    }
    for out in outputs:
        cols[out] = np.asarray(V[out])
    if "vo" in V:
        cols["vo"] = np.asarray(V["vo"])
    for name, key in object_columns.items():
        if key in V:
            cols[name] = np.asarray(V[key])
 
    n = len(cols["lat"])
    bad = {k: len(v) for k, v in cols.items() if len(v) != n}
    if bad:
        raise ValueError(f"Per-point arrays have inconsistent lengths (expected {n}): {bad}")
 
    points = pd.DataFrame(cols).sort_values(["track_id", "time"], kind="stable").reset_index(drop=True)
 
    nan_frac = points[outputs].isna().mean()
    if (nan_frac > 0).any():
        warnings.warn("NaNs in sampled fields (check time/space coverage):\n"
                      + nan_frac[nan_frac > 0].to_string())
 
    # --- per-track summary --------------------------------------------------
    named = {out: (out, method if method in ("min", "max") else "mean")
             for _, method, out in samples}
    if "vo" in points:
        named["vo_max"] = ("vo", "max")
    if extra_agg:
        named.update(extra_agg)
 
    g = points.groupby("track_id", sort=False)
    summary = g.agg(**named)
    summary["n_points"] = g.size()
    summary["start_time"] = g["time"].first()
    summary["end_time"] = g["time"].last()
    summary["duration_h"] = (summary["end_time"] - summary["start_time"]) / pd.Timedelta(hours=1)
    summary["lat_genesis"] = g["lat"].first()
    summary["lon_genesis"] = g["lon"].first()
    summary["lat_lysis"] = g["lat"].last()
    summary["lon_lysis"] = g["lon"].last()
    if "vo" in points and points["vo"].notna().any():
        idx = g["vo"].idxmax()
        summary["time_of_vo_max"] = pd.Series(points.loc[idx.values, "time"].values, index=idx.index)
 
    return points, summary.reset_index()
