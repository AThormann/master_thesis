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

def hodges_tracker(input_file_path='name', output_file_name='output', track_variable='vo', min_track_points=10, taper_points=10, lmin=5,lmax=42,dmax=6.5,min_object_grid_points=18,
                   object_threshold=None, use_dmax_zones=True, exclude_boundary_extrema=False):
    """
    Run the Hodges tracker on a file path or an already prepared xarray Dataset.

    dmax (degrees per time step) is only used when use_dmax_zones=False; otherwise
    pystormtracker uses its default latitude zones (6.5 deg in the extratropics).
    """
    if isinstance(input_file_path, xr.Dataset):
        data_file = input_file_path
    else:
        data_file = xr.open_dataset(input_file_path)

    extra = {} if use_dmax_zones else {"dmax_zones": np.empty((0, 5))}
    tracker = pst.HodgesTracker(min_track_points =min_track_points, taper_points=taper_points, lmin=lmin, lmax=lmax, dmax=dmax, min_object_grid_points=min_object_grid_points,
                                exclude_boundary_extrema=exclude_boundary_extrema, **extra)

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
        #path = "PyStormTracker/json"
        if output_file_name:
            tracks.write(f"{output_file_name}.trackjson")
        "Tracking done"
        return tracks


def plot_tracks(data, min_length=0, projection=ccrs.PlateCarree(), points=True, text=True, title=None, n_rows=1, n_cols=1):
    """Plot tracks on a new map, or on an existing cartopy axis given as `ax` (e.g. a subplot)."""
    fig, ax = plt.subplots(n_rows, n_cols, figsize=(10, 8), subplot_kw={'projection': projection})
    
    ax.coastlines()
    ax.set_extent([np.min(data.lons), np.max(data.lons), np.min(data.lats), np.max(data.lats)], crs=projection)
    ax.gridlines(draw_labels=True)

    for d in data:
        if len(d) >= min_length:
            # Geodetic: segments follow great circles, so tracks crossing the
            # dateline (179 -> -179) are not drawn as a line across the whole map
            ax.plot(d.lons, d.lats, transform=ccrs.Geodetic())
            if points == True:
                ax.plot(d.lons, d.lats, 'go', transform=projection)  # Track points
            ax.plot(d.lons[0], d.lats[0], 'k', transform=projection)  # Start point
            if text == True:
                ax.text(d.lons[0], d.lats[0], str(pd.to_datetime(d.times[0], unit='ms')), transform=projection)  # Track start time
    ax.set_title(f"ETC tracks in NA between {pd.to_datetime(data[0].times[0], unit="ms")} and {pd.to_datetime(data[-1].times[-1], unit="ms")}")

    return fig, ax





def load_st_report(file_path, variable="msl", units="Pa", mode="min", scale=100.0):
    """
    Read a storm report text file (e.g. Data/st_report.txt) into a Tracks object,
    so it can be used like the output of hodges_tracker / pst.load_tracks.

    The file is a list of blocks, one per start time:

        1990/01/01  00   Number of Storms Begun  2
             1   1     93.93   -5.24   1007.38    1990/01/01  00
             ...

    with the columns: storm number (within the block), point number, longitude,
    latitude, central pressure, date, hour. Storm numbers restart in every block,
    so tracks get new running ids 1, 2, 3, ... in file order.

    variable / units / mode describe the value column after multiplying it by
    `scale`. The file has central pressure in hPa; pystormtracker wants 'msl'
    in Pa, so the default scale is 100.
    """
    lons, lats, values, times, offsets = [], [], [], [], []
    n_announced = 0
    block = 0          # counts header lines, to separate storms with the same number
    current = None     # (block, storm number) of the track being read

    with open(file_path) as f:
        for line_no, line in enumerate(f, start=1):
            if not line.strip():
                continue
            parts = line.split()
            if "Number of Storms Begun" in line:
                block += 1
                n_announced += int(parts[-1])
                continue
            if len(parts) != 7:
                raise ValueError(f"{file_path}, line {line_no}: cannot read {line!r}")

            storm = (block, int(parts[0]))
            if storm != current:
                current = storm
                offsets.append(len(lons))
            lons.append(float(parts[2]))
            lats.append(float(parts[3]))
            values.append(float(parts[4]))
            times.append(f"{parts[5].replace('/', '-')}T{parts[6]}")

    n_tracks = len(offsets)
    offsets.append(len(lons))
    if n_tracks != n_announced:
        warnings.warn(f"The headers announce {n_announced} storms but {n_tracks} were read.")

    times_ms = np.array(times, dtype="datetime64[ms]").astype(np.int64)
    metadata = pst.models.tracks.TracksMetadata(primary_variable=variable, mode=mode,
                                                units={variable: units})
    return pst.Tracks(ids=np.arange(1, n_tracks + 1), offsets=np.array(offsets),
                      times=times_ms, lats=np.array(lats), lons=np.array(lons),
                      variables={variable: np.array(values) * scale}, metadata=metadata)


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


# --------------------------------------------------------------------------
# Regional tracking (North Atlantic / Northern Europe) and post-processing
# --------------------------------------------------------------------------

EARTH_RADIUS_KM = 6371.0


def prepare_tracking_field(ds, variable="vo", level=850, domain=None):
    """
    Return a clean (time, latitude, longitude) Dataset with only the tracking variable.

    - selects the pressure level (if there is one)
    - converts longitudes to -180..180 and sorts them (ERA5 is often 0..360)
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


def haversine_km(lat1, lon1, lat2, lon2):
    """Great-circle distance in km (works on numpy arrays)."""
    lat1, lon1, lat2, lon2 = map(np.radians, (lat1, lon1, lat2, lon2))
    a = np.sin((lat2 - lat1) / 2) ** 2 + np.cos(lat1) * np.cos(lat2) * np.sin((lon2 - lon1) / 2) ** 2
    return 2 * EARTH_RADIUS_KM * np.arcsin(np.sqrt(a))


def _in_box(lat, lon, box):
    lon_min, lon_max, lat_min, lat_max = box
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


def postprocess_tracks(tracks, region, domain, time_range=None,
                       region_criterion="passes", min_points_in_region=1,
                       min_duration_h=48, min_separation_km=1000,
                       min_peak=None, max_frac_at_edge=0.5,
                       edge_buffer_deg=2.5, drop_time_truncated=False):
    """
    Filter raw Hodges tracks down to the storms we want in the region.

    Every criterion becomes a boolean column in `summary` (True = passes), and
    `keep` is the AND of all of them, so you can see why each track was removed.

    region_criterion : "passes"  -> at least `min_points_in_region` points in the region
                       "genesis" -> the first point is in the region
                       "peak"    -> the point of maximum intensity is in the region
    min_duration_h    : minimum lifetime in hours
    min_separation_km : minimum distance between genesis and lysis (removes stationary features)
    min_peak          : minimum peak intensity of the tracked variable (vo: s-1, msl: Pa; msl uses <=)
    max_frac_at_edge  : maximum fraction of points within `edge_buffer_deg` of the domain edge
    drop_time_truncated : also remove tracks touching the first/last time step

    Returns (filtered Tracks, summary DataFrame for all tracks, per-point DataFrame for all tracks)
    """
    summary, points = summarise_tracks(tracks, region, domain, time_range, edge_buffer_deg)
    var = tracks.primary_variable

    if region_criterion == "passes":
        summary["ok_region"] = summary["n_in_region"] >= min_points_in_region
    elif region_criterion == "genesis":
        summary["ok_region"] = summary["genesis_in_region"]
    elif region_criterion == "peak":
        summary["ok_region"] = summary["peak_in_region"]
    else:
        raise ValueError("region_criterion must be 'passes', 'genesis' or 'peak'")

    summary["ok_duration"] = summary["duration_h"] >= min_duration_h
    summary["ok_separation"] = summary["separation_km"] >= min_separation_km
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


def plot_regional_tracks(tracks, region, domain, summary=None, ax=None, title=None):
    """
    Plot tracks with the region (red) and tracking domain (black dashed) boxes.
    If `summary` is given, removed tracks are drawn in light grey.
    """
    proj = ccrs.LambertConformal(central_longitude=-15, central_latitude=55)
    pc = ccrs.PlateCarree()
    if ax is None:
        fig, ax = plt.subplots(figsize=(11, 8), subplot_kw={"projection": proj})
    ax.set_extent([domain[0] - 2, domain[1] + 2, domain[2] - 2, domain[3]], crs=pc)
    ax.coastlines(linewidth=0.6)
    ax.gridlines(draw_labels=True, linewidth=0.3)

    for box, style in ((domain, dict(color="k", linestyle="--")), (region, dict(color="red"))):
        # densify the edges so they follow parallels/meridians in the projection
        lo = np.linspace(box[0], box[1], 100)
        la = np.linspace(box[2], box[3], 100)
        xs = np.concatenate([lo, np.full(100, box[1]), lo[::-1], np.full(100, box[0])])
        ys = np.concatenate([np.full(100, box[2]), la, np.full(100, box[3]), la[::-1]])
        ax.plot(xs, ys, lw=1.5, transform=pc, **style)

    keep = None if summary is None else dict(zip(summary["track_id"], summary["keep"]))
    for tr in tracks:
        kept = True if keep is None else keep.get(tr.track_id, True)
        if kept:
            ax.plot(tr.lons, tr.lats, "-", lw=1.4, transform=pc)
            ax.plot(tr.lons[0], tr.lats[0], "k.", ms=6, transform=pc)
        else:
            ax.plot(tr.lons, tr.lats, "-", color="0.75", lw=0.8, transform=pc, zorder=0)
    ax.set_title(title or "Tracks")
    return ax
