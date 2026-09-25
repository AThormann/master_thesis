import numpy as np
import xarray as xr
import netCDF4 as nc
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
        tracks.write(f"{output_file_name}.trackjson")
        return "Tracking done"



