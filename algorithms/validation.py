import numpy as np
import xarray as xr

def wrap_lon(lon, to="-180"):
    """Convert longitudes between 0..360 and -180..180 conventions."""
    if to == "-180":
        return ((lon + 180) % 360) - 180
    return lon % 360

def load_tracks(path):
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

    ds =  ds.to_pandas()
    return ds