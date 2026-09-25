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
        tracks.write(f"{output_file_name}.trackjson")
        return "Tracking done"




def plot_tracks(data, min_length=0, projection=ccrs.PlateCarree()):
    fig, ax = plt.subplots(figsize=(10, 8), subplot_kw={'projection': projection})
    ax.coastlines()
    ax.set_extent([np.min(data.lons), np.max(data.lons), np.min(data.lats), np.max(data.lats)], crs=projection)
    ax.gridlines(draw_labels=True)

    for d in data:
        if len(d) >= min_length:
            ax.plot(d.lons, d.lats, transform=projection)
            ax.plot(d.lons[0], d.lats[0], 'go', transform=projection)  # Start point
            ax.text(d.lons[0], d.lats[0], str(pd.to_datetime(d.times[0], unit='ms')), transform=projection)  # Track start time
    ax.set_title(f"ETC tracks in NA between {pd.to_datetime(data[0].times[0], unit="ms")} and {pd.to_datetime(data[-1].times[-1], unit="ms")}")







#plotting method for tracks in json format
#projection


