
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import pandas as pd
import pystormtracker as pst

from pystormtracker.models.tracks import Tracks, TracksMetadata
from pystormtracker.metrics.compare import compare_tracks, TrackComparisonConfig

from helpers import haversine_km, _in_box, _to_datetime, wrap_lon
from configs import MISSING, COLS




def load_tracks(path, start_time, end_time):
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
    ds["lon"] = wrap_lon(ds.lon, to="-180")
    ds = ds.to_pandas()
    times = pd.date_range(start=start_time, end=end_time, freq="6h")
    smin, smax = ds["step"].min(), ds["step"].max()
    if smin < 1 or smax > len(times):
        raise ValueError(f"Track steps span {smin}..{smax}, but {start_time} to {end_time} "
                         f"has only {len(times)} time steps (valid: 1..{len(times)}).")

    ds['time'] = times[ds['step']-1] # times are created from the step in the data set and step 1 is equal to the start time
    
    return ds


def dataframe_to_tracks(df, variable="vo", units="s^-1", mode="max", scale=1.0):
    """
    Convert a per-point DataFrame (track_id, time, lat, lon, <column>) to a Tracks object.
 
    variable : name the intensity gets in the Tracks object (use the same name as your
               own tracks, e.g. "vo", so intensities can be compared)
    scale    : multiply the intensity by this. TRACK files usually store vorticity in
               units of 1e-5 s-1, so use scale=1e-5 to get s-1 like ERA5 'vo'.
    The intensity is read from the column `variable`, or 'vorticity' if that is absent.
    """
    col = variable if variable in df.columns else "vorticity"
    df = df.dropna(subset=["lat", "lon", "time"]).sort_values(["track_id", "time"], kind="stable")
    counts = df.groupby("track_id", sort=False).size()
    offsets = np.concatenate([[0], np.cumsum(counts.to_numpy())])
    times_ms = (pd.DatetimeIndex(df["time"]).as_unit("ms").asi8).astype(np.int64)
    meta = TracksMetadata(primary_variable=variable, mode=mode, units={variable: units})
    return Tracks(
        ids=counts.index.to_numpy(dtype=np.int64),
        offsets=offsets,
        times=times_ms,
        lats=df["lat"].to_numpy(float),
        lons=df["lon"].to_numpy(float),
        variables={variable: df[col].to_numpy(float) * scale},
        metadata=meta,
    )



def clip_to_domain(df, domain, time_range):
    lon_min, lon_max, lat_min, lat_max = domain
    t0, t1 = map(pd.Timestamp, time_range)
    m = (df["lon"].between(lon_min, lon_max) & df["lat"].between(lat_min, lat_max)
         & df["time"].between(t0, t1))
    out = df[m].sort_values(["track_id", "time"]).copy()
    gap = out.groupby("track_id")["time"].diff() > pd.Timedelta(hours=6)
    out["segment"] = gap.groupby(out["track_id"]).cumsum()
    out["track_id"] = out["track_id"] * 100 + out["segment"]   # unique id per piece
    return out


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


def compare_to_validation(reference, candidate, matching="global_assignment",
                          max_mean_separation_deg=2.0, min_overlap_fraction=0.6,
                          variable=None):
    """
    Run compare_tracks and return (result, per-match DataFrame).
 
    matching : "nearest"           -> TRACK's own rule (a candidate can be reused)
               "mutual_nearest"    -> reciprocal nearest pairs only
               "global_assignment" -> one-to-one; gives tp/fp/fn, precision, recall, f1
    variable : intensity variable to compare along matched pairs (units must be identical
               in both Tracks objects), or None to skip intensity statistics.
    """
    cfg = TrackComparisonConfig(matching=matching,
                                max_mean_separation_deg=max_mean_separation_deg,
                                min_overlap_fraction=min_overlap_fraction,
                                variable=variable)
    res = compare_tracks(reference, candidate, config=cfg)
    rows = []
    for m in res.matches:
        row = {
            "reference_id": m.reference_id, "candidate_id": m.candidate_id,
            "overlap_fraction": m.overlap_fraction, "overlap_count": m.overlap_count,
            "mean_sep_km": m.mean_separation_km, "max_sep_km": m.maximum_separation_km,
            "ref_points": m.reference.point_count, "cand_points": m.candidate.point_count,
            "ref_peak": m.reference.peak_intensity, "cand_peak": m.candidate.peak_intensity,
            "n_eligible_candidates": m.eligible_candidate_count,
        }
        if m.intensity_difference is not None:
            row.update(intensity_bias=m.intensity_difference.bias,
                       intensity_rmse=m.intensity_difference.rmse)
        rows.append(row)
    return res, pd.DataFrame(rows)




