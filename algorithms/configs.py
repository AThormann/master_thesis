import numpy as np
import pandas as pd


EARTH_RADIUS_KM = 6371.0
MISSING = 1.0e25
#format (lon_min, lon_max, lat_min, lat_max) or (west, east, south, north)
REGION = (-80, 20, 30, 80)   # smaller region for tracking, limited by the domain, 10 deg smaller

DOMAIN = (-90, 30, 20, 90)   # no longer limits the tracks (tracks are full length); only used for plotting and the (disabled) edge filter
COLS = ["step", "lon", "lat", "vorticity", "mslp_lon", "mslp_lat", "mslp"]
TIME_RANGE = None  # set to (start_time, end_time) if you want to filter tracks by time


# for processing raw tracks in TRACK object format 
POSTPROCESS_CONFIG = dict(region=REGION,
                          domain=DOMAIN,              # only used by the edge filter, which is now off
                          time_range=TIME_RANGE,
                          region_criterion=("passes", "peak"),  # both must hold
                          min_points_in_region=4,     # "passes": at least 4 points (24 h) in REGION
                                                      # "peak": max vo must lie in REGION
                          rsplice_min_points=8,       # 48 hours
                          rsplice_max_points=None,
                          rsplice_distance_degrees=10.0,  # reproducing Priestley
                          rsplice_distance_mode="endpoint",
                          rsplice_cadence=np.timedelta64(6, "h"),
                          min_peak=None,
                          max_frac_at_edge=None,      # no edge filter: tracking is global, tracks are full length
                          edge_buffer_deg=2.5,
                          drop_time_truncated=False,
)

TRACK_CONFIGS = dict(min_track_points=3, 
                     taper_points=0, 
                     lmin=6,
                     lmax=42, 
                     spectral_taper = 0.1, 
                     min_object_grid_points=1, 
                     exclude_boundary_extrema=False)

track_variable = 'vo'
object_threshold=1e-5

# Seasons (DJF, MAM, JJA, SON): the season files (tracking field + single levels) cover the season
# plus this many days on either side. A track belongs to the season if its peak vorticity lies in
# the season itself (experiments/download_all/seasons.py, postprocess_on_seasons.py).
SEASON_BUFFER_DAYS = 10

# ERA5 single levels: area of the CDS request (north, west, south, east).
# Larger than DOMAIN so that sampling within radius_km of any track point stays inside the field.
SL_AREA = [90, -120, 10, 40]

# Sampling of fields along the tracks (statistics_on_tracks.py)
SAMPLE_CONFIG = dict(radius_km=500,
                     # (field, method, output column) from the single-level file; method: nearest, bilinear, mean, max, min
                     samples=[("msl",  "min",     "msl_min"),      # central pressure (Pa)
                              ("msl",  "nearest", "msl_center"),   # pressure at the vorticity centre (Pa)
                              ("ws10", "max",     "ws10_max"),     # ws10 = hypot(u10, v10) on the grid (m s-1)
                              ("ws10", "mean",    "ws10_mean"),
                              ("gust10", "max",   "gust10_max"),
                              ("gust10", "mean",  "gust10_mean")],  # gust10 = max wind gust on the grid (m s-1)
                     # same, from the 0.25 degree pressure-level file (850 hPa)
                     pl_samples=[("vo",    "max", "vo850_unfilt_max"),  # unfiltered vorticity (s-1)
                                 ("ws850", "max", "ws850_max")],        # ws850 = hypot(u, v) on the grid (m s-1)
                     # outputs (max/min samples) whose grid location is also stored (<output>_lat/_lon) and
                     # summarised per track (location of the maximum 10 m wind)
                     locate=["ws10_max"],
                     sample_res=0.25,  # grid of the sampled fields, 0.25 (native) or one of regrid_era5.RESOLUTIONS
)
