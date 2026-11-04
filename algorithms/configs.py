import numpy as np
import pandas as pd


EARTH_RADIUS_KM = 6371.0
MISSING = 1.0e25
#format (lon_min, lon_max, lat_min, lat_max) or (west, east, south, north)
REGION = (-80, 20, 30, 80)   # smaller region for tracking, limited by the domain, 10 deg smaller

DOMAIN = (-90, 30, 20, 90)   # tracking domain (region + buffer), limited by the downloaded data
COLS = ["step", "lon", "lat", "vorticity", "mslp_lon", "mslp_lat", "mslp"]
TIME_RANGE = None  # set to (start_time, end_time) if you want to filter tracks by time


# for processing raw tracks in TRACK object format 
POSTPROCESS_CONFIG = dict(region=REGION,
                          domain=DOMAIN,
                          time_range=TIME_RANGE,
                          region_criterion="passes",
                          min_points_in_region=1,
                          rsplice_min_points=8, #42 hours
                          rsplice_max_points=None,
                          rsplice_distance_degrees=10.0, #Reproducing priestley
                          rsplice_distance_mode="endpoint",
                          rsplice_cadence=np.timedelta64(6, "h"),
                          min_peak=None,
                          max_frac_at_edge=0.5,
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