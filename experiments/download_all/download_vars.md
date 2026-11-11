## Downloading ERA5 single levels

# Goal
 I want to use the tracks found on the ERA5 850 hPa pressure levels to extract other variables along the tracks to create distributions. 

# Steps
 1. Rewrite statistics_on_tracks.py to be compatible with processing.py, helpers.py and configs.py. Remove the default variables for t2m min and max. Add precipitation and vorticity. Vorticity should come from the track itself. check if there is both a scaled vorticity variable and a raw vorticity connected to the track. There should be. Check/rewrite statistics_on_tracks.py shuch that I get variables for all track points and not just summaries.
 2. Add the file to the work flow of ERA5 download seasons. Add a function for an API call to download the correct field from single levels in the same fashion as the vorticity from the pressure levels.
 3. Organize the files such that the output trackjson and the statistics of the tracks are in the same directory with clear names indicating the season it belongs to.
 4. Write appropiate checks that ensures that variables are sampled and that the times are correct

# API call
 This call downloads one year of data. This is what should be used in step 2:
 The API call is shown for 1940 but should function like ERA5_download_seasons.py

 API call:
 import cdsapi
 dataset = "reanalysis-era5-single-levels"
 request = {
    "product_type": ["reanalysis"],
    "variable": [
        "10m_u_component_of_wind",
        "10m_v_component_of_wind",
        "mean_sea_level_pressure",
        "total_precipitation"
    ],
    "year": ["1940"],
    "month": ["01", "02", "12"],
    "day": [
        "01", "02", "03",
        "04", "05", "06",
        "07", "08", "09",
        "10", "11", "12",
        "13", "14", "15",
        "16", "17", "18",
        "19", "20", "21",
        "22", "23", "24",
        "25", "26", "27",
        "28", "29", "30",
        "31"
    ],
    "time": [
        "00:00", "06:00", "12:00",
        "18:00"
    ],
    "data_format": "netcdf",
    "download_format": "unarchived",
    "area": [90, -120, 20, 40]
 }
 client = cdsapi.Client()
 client.retrieve(dataset, request).download()


# Notes:
The area is limited as the tracks all pass through the box indicated in the configs.py I'm only interested in variables connected to the tracks in this region so I've set a slightly larger box that should cover the tracks. If a track is outside the region a warning should be raised and the track should have NaN values.
Give me suggestions on how to best conduct the sampling of variables and what extra (dervied) variables I could extract from the fields. 
Point out what steps might be missing.
