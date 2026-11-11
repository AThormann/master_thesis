#File containing functions to sample fields (e.g. ERA5 single levels) along tracks and to
#derive per-point quantities from the tracks. Used with the filtered tracks from processing.py.
import os
import warnings
import numpy as np
import pandas as pd
import xarray as xr
import pystormtracker as pst
from helpers import tracks_to_points, haversine_km, wrap_lon
from configs import SAMPLE_CONFIG, SL_AREA


TIME_NAMES = ("valid_time", "time")
PA_FIELDS = ("msl",)          # sampled pressure fields are converted from Pa to hPa in the output
PHYSICAL_RANGES = {           # column prefix -> (min, max); the longest matching prefix is used
    "msl": (900.0, 1070.0),   # hPa
    "ws": (0.0, 90.0),        # m s-1
    "gust": (0.0, 100.0),     # m s-1
    "vo": (-1e-3, 1e-3),      # s-1, T6-42 filtered
    "vo850_unfilt": (-5e-3, 5e-3),   # s-1, unfiltered 0.25 deg (fronts, orography)
}




# PREPARE FIELDS

def _time_name(ds):
    name = next((n for n in TIME_NAMES if n in ds.coords), None)
    if name is None:
        raise ValueError(f"No time coordinate {TIME_NAMES} in dataset (coords: {list(ds.coords)}).")
    return name


def add_derived_fields(ds):
    """
    Wind speed on the grid (never from averaged u and v): ws10 from u10/v10, ws850 from u/v.
    gust10 is the ERA5 instantaneous 10 m wind gust (stored as i10fg in the CDS netCDF files).
    """
    if "u10" in ds and "v10" in ds and "ws10" not in ds:
        ds = ds.assign(ws10=np.hypot(ds["u10"], ds["v10"]).assign_attrs(units="m s-1"))
    if "i10fg" in ds and "gust10" not in ds:
        ds = ds.assign(gust10=ds["i10fg"].assign_attrs(units="m s-1"))
    if "u" in ds and "v" in ds and "ws850" not in ds:
        ds = ds.assign(ws850=np.hypot(ds["u"], ds["v"]).assign_attrs(units="m s-1"))
    return ds


def check_times(track_times, ds):
    """Every track time must be a time step of the dataset (no silent nearest-time matching)."""
    ds_times = pd.DatetimeIndex(ds[_time_name(ds)].values)
    track_times = pd.DatetimeIndex(track_times).unique()
    missing = track_times.difference(ds_times)
    if len(missing):
        raise ValueError(f"{len(missing)} of {len(track_times)} track times are not in the dataset "
                         f"({ds_times.min()} -> {ds_times.max()}), e.g. {list(missing[:3])}")


def prepare_field(ds, fields, track_times, level=850):
    """
    Open the dataset, add derived fields, keep `fields`, select the pressure level (if any)
    and the track times, and load it into memory.
    """
    if isinstance(ds, (str, os.PathLike)):
        ds = xr.open_dataset(ds)
    if "pressure_level" in ds.dims:
        ds = ds.sel(pressure_level=level, method="nearest")
    ds = add_derived_fields(ds)
    missing = sorted(set(fields) - set(ds.data_vars))
    if missing:
        raise KeyError(f"Fields {missing} not in dataset (available: {list(ds.data_vars)}).")
    check_times(track_times, ds)
    time_name = _time_name(ds)
    times = pd.DatetimeIndex(track_times).unique().sort_values()
    # read one contiguous block (fast) and pick the track times in memory
    # (selecting scattered times directly from a compressed file is very slow)
    block = ds[list(fields)].sel({time_name: slice(times[0], times[-1])}).load()
    return block.sel({time_name: times.values})




# SAMPLING

def area_mask(lats, lons, area, radius_km):
    """
    True where the whole sampling circle (radius_km) of a point lies inside area = (north, west, south, east).
    Points too close to (or outside) the edge of a limited-area field would get truncated or edge values.
    """
    lats = np.asarray(lats, dtype=float)
    lons = wrap_lon(np.asarray(lons, dtype=float), to="-180")
    north, west, south, east = area
    dlat = radius_km / 111.0
    dlon = radius_km / (111.0 * np.maximum(np.cos(np.radians(lats)), 1e-6))
    return ((lats + dlat <= north) & (lats - dlat >= south)
            & (lons - dlon >= west) & (lons + dlon <= east))


def sample_fields(points, ds, samples, radius_km, mask=None):
    """
    Sample (field, method, output) from ds at every point of the per-point DataFrame
    (columns time, lat, lon). Same methods as pystormtracker.sample.sample_tracks, but vectorised
    per point (sample_tracks computes distances one grid point at a time, far too slow at 0.25 deg):
      nearest       value at the nearest grid point
      mean/max/min  over all grid points within radius_km (great-circle distance)
    Points where mask is False get NaN. Pressure fields (PA_FIELDS) are returned in hPa.
    Returns {output: values in the row order of points}.
    """
    for field, method, name in samples:
        if method not in ("nearest", "mean", "max", "min"):
            raise ValueError(f"Unsupported sampling method '{method}' for '{name}'")
    time_name = _time_name(ds)
    lats = np.asarray(ds["latitude"].values, dtype=float)
    lons = np.asarray(ds["longitude"].values, dtype=float)
    dlat_grid = np.abs(np.diff(lats)).max() if lats.size > 1 else 0.0
    dlon_grid = np.abs(np.diff(lons)).max() if lons.size > 1 else 0.0
    fields = {f: np.asarray(ds[f].values, dtype=float) for f in {s[0] for s in samples}}   # (time, lat, lon)
    t_index = pd.Index(pd.DatetimeIndex(ds[time_name].values)).get_indexer(pd.DatetimeIndex(points["time"]))
    if (t_index < 0).any():
        raise ValueError("Track times missing in the dataset (run check_times first)")

    n = len(points)
    out = {name: np.full(n, np.nan) for _, _, name in samples}
    p_lat, p_lon = points["lat"].to_numpy(float), points["lon"].to_numpy(float)
    for i in range(n):
        if mask is not None and not mask[i]:
            continue
        lat0, lon0 = p_lat[i], p_lon[i]
        # window around the point that surely contains the circle and the nearest grid point
        win_lat = radius_km / 111.0 + dlat_grid
        win_lon = min(radius_km / (111.0 * max(np.cos(np.radians(lat0)), 1e-6)) + dlon_grid, 180.0)
        rows = np.flatnonzero(np.abs(lats - lat0) <= win_lat)
        cols = np.flatnonzero(np.abs((lons - lon0 + 180.0) % 360.0 - 180.0) <= win_lon)
        if rows.size == 0 or cols.size == 0:
            continue
        dist = haversine_km(lat0, lon0, lats[rows][:, None], lons[cols][None, :])
        in_circle = dist <= radius_km
        nearest = np.unravel_index(np.argmin(dist), dist.shape)
        t = t_index[i]
        for field, method, name in samples:
            sub = fields[field][t][np.ix_(rows, cols)]
            if method == "nearest":
                out[name][i] = sub[nearest]
            elif in_circle.any():
                vals = sub[in_circle]
                out[name][i] = {"mean": np.nanmean, "max": np.nanmax, "min": np.nanmin}[method](vals)
    for field, _, name in samples:
        if field in PA_FIELDS:
            out[name] = out[name] / 100.0   # Pa -> hPa
    return out


def add_track_kinematics(points):
    """
    Per-point quantities derived from the track itself (first point of a track is NaN):
      speed_kmh       propagation speed from the previous point
      direction_deg   direction of motion from the previous point (0 = north, 90 = east)
      dvo_dt          change of the tracked vorticity per 6 h (s-1 / 6 h)
      dmsl_24h        change of msl_min over the previous 24 h (hPa), if msl_min was sampled
      deepening_bergeron  24 h deepening of msl_min in Bergeron units (>= 1: explosive cyclone,
                      Sanders & Gyakum 1980: 24 hPa/24 h * sin(lat)/sin(60))
    """
    p = points.sort_values(["track_id", "time"], kind="stable").reset_index(drop=True)
    g = p.groupby("track_id", sort=False)
    prev_lat, prev_lon, prev_time = g["lat"].shift(), g["lon"].shift(), g["time"].shift()
    dt_h = (p["time"] - prev_time) / pd.Timedelta(hours=1)

    p["speed_kmh"] = haversine_km(prev_lat, prev_lon, p["lat"], p["lon"]) / dt_h
    lat1, lat2 = np.radians(prev_lat), np.radians(p["lat"])
    dlon = np.radians(p["lon"] - prev_lon)
    p["direction_deg"] = np.degrees(np.arctan2(np.sin(dlon) * np.cos(lat2),
                                               np.cos(lat1) * np.sin(lat2)
                                               - np.sin(lat1) * np.cos(lat2) * np.cos(dlon))) % 360
    if "vo" in p:
        p["dvo_dt"] = (p["vo"] - g["vo"].shift()) / (dt_h / 6.0)

    if "msl_min" in p:
        # value 24 h earlier in the same track (NaN if that time is not part of the track)
        earlier = p[["track_id", "time", "lat", "msl_min"]].copy()
        earlier["time"] = earlier["time"] + pd.Timedelta(hours=24)
        earlier = earlier.rename(columns={"lat": "lat_24h", "msl_min": "msl_min_24h"})
        p = p.merge(earlier, on=["track_id", "time"], how="left")
        p["dmsl_24h"] = p["msl_min"] - p["msl_min_24h"]
        mean_lat = np.radians(0.5 * (p["lat"] + p["lat_24h"]))
        p["deepening_bergeron"] = (-p["dmsl_24h"] / 24.0) * np.sin(np.radians(60.0)) / np.sin(mean_lat)
        p = p.drop(columns=["lat_24h", "msl_min_24h"])
    return p


def check_sampled(points, columns, inside):
    """No NaN for points inside the area, and values in a physical range (warning otherwise)."""
    for col in columns:
        n_nan = int(points.loc[inside, col].isna().sum())
        if n_nan:
            raise ValueError(f"{n_nan} NaN values in '{col}' for points inside the sampling area")
        prefixes = [k for k in PHYSICAL_RANGES if col.startswith(k)]
        if prefixes:
            lo, hi = PHYSICAL_RANGES[max(prefixes, key=len)]
            bad = points[col].notna() & ~points[col].between(lo, hi)
            if bad.any():
                warnings.warn(f"{int(bad.sum())} values of '{col}' outside [{lo}, {hi}] "
                              f"(min {points[col].min():.4g}, max {points[col].max():.4g})")


def sample_along_tracks(tracks, sl_ds, pl_ds=None, config=SAMPLE_CONFIG, area=SL_AREA):
    """
    Sample fields along every point of the tracks.

    tracks : Tracks or path to a .trackjson file (e.g. the filtered tracks of one season)
    sl_ds  : single-level dataset or path (u10, v10, msl over `area`), sampled with config["samples"]
    pl_ds  : optional 850 hPa dataset or path (global, u, v, vo), sampled with config["pl_samples"]
    area   : (north, west, south, east) of sl_ds. Points whose sampling circle is not inside it
             get NaN for the single-level samples and a warning is raised.

    Returns a DataFrame with one row per track point: track_id, time, lat, lon, every variable
    stored on the tracks (vo, and raw_value / object_* if present), the sampled columns
    (msl in hPa, wind in m s-1, vorticity in s-1), the column inside_area and the kinematics
    from add_track_kinematics.
    """
    if isinstance(tracks, (str, os.PathLike)):
        tracks = pst.load_tracks(str(tracks))
    if tracks is None or len(tracks) == 0:
        raise ValueError("No tracks to sample.")
    radius = config["radius_km"]

    points = tracks_to_points(tracks)
    track_times = pd.DatetimeIndex(points["time"])

    inside = area_mask(points["lat"], points["lon"], area, radius)
    points["inside_area"] = inside
    if not inside.all():
        n_tracks = points.loc[~inside, "track_id"].nunique()
        warnings.warn(f"{int((~inside).sum())} points of {n_tracks} tracks are (within {radius} km) outside "
                      f"the single-level area {area}; their single-level samples are NaN")

    samples = list(config["samples"])
    sl = prepare_field(sl_ds, {s[0] for s in samples}, track_times)
    sampled = sample_fields(points, sl, samples, radius, mask=inside)

    pl_samples = list(config.get("pl_samples") or []) if pl_ds is not None else []
    if pl_samples:
        pl = prepare_field(pl_ds, {s[0] for s in pl_samples}, track_times)
        sampled.update(sample_fields(points, pl, pl_samples, radius))

    for name, values in sampled.items():
        points[name] = values
    check_sampled(points, [s[2] for s in samples], inside)
    check_sampled(points, [s[2] for s in pl_samples], np.ones(len(points), dtype=bool))
    return add_track_kinematics(points)


def summarise_points(points):
    """One row per track: lifetime, genesis/lysis and extremes of the per-point columns."""
    g = points.groupby("track_id", sort=False)
    s = pd.DataFrame({
        "n_points": g.size(),
        "start_time": g["time"].first(),
        "end_time": g["time"].last(),
        "lat_genesis": g["lat"].first(),
        "lon_genesis": g["lon"].first(),
        "lat_lysis": g["lat"].last(),
        "lon_lysis": g["lon"].last(),
    })
    s["duration_h"] = (s["end_time"] - s["start_time"]) / pd.Timedelta(hours=1)
    skip = {"track_id", "time", "lat", "lon", "inside_area", "direction_deg"}
    for col in points.columns:
        if col in skip or not pd.api.types.is_numeric_dtype(points[col]):
            continue
        agg = "min" if col.startswith(("msl", "dmsl")) else "max"
        s[f"{col}_{agg}"] = g[col].agg(agg)
    if "ws10_mean" in points:
        s["ws10_mean_mean"] = g["ws10_mean"].mean()
    return s.reset_index()
