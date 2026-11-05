## Downloading and regridding data

I want a process in a .py file that downloads yearly ERA5 data, concatenates the correct times and saves a .nc file. 

Rules:
1. One .py file to complete the download
2. one method or .py file to do the regridding

# ERA5 reanalysis data
From the Copernicus website one can download data with an API request. I want to download data the pressure levels data set. We are interested in 850 hPa relative vorticity and U- and V- components of the wind at the same pressure level.
It is important to note that the request size has a limit so I would need to loop over the API request to download the data sets.
It is possible to choose hours, days, months and years when downloading.
We want 6 hour resolution meaning times (00, 06, 12, 18), on all days in the months. december, january and february. We want these months for all years.

This is the API request:
import cdsapi

dataset = "reanalysis-era5-pressure-levels"
request = {
    "product_type": ["reanalysis"],
    "variable": [
        "u_component_of_wind",
        "v_component_of_wind",
        "vorticity"
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
    "pressure_level": ["850"],
    "data_format": "netcdf",
    "download_format": "unarchived"
}

client = cdsapi.Client()
client.retrieve(dataset, request).download()

 

# Structuring of the datasets
We are interested in the winter season meaning that we need to concatenate two years of data each containing january, febuary and december within the year. I want an output datafile containing december of the earliest year and january and febuary of the latest year.
We want this for all years so requesting and downloading two years and deleting the remainder is a waste of downloading.
We propose downloading all years first then concatenating and removing the old yearly files.
Thus we will have a folder for each winter season since 1940. 1940 wil contain december of 1940 and january and febuary of 1941.
The upper edge will be 2025 as there is no data for 2026.

The folders should be located at this path: /Users/asgerthormann/Documents/uni/master/master_code/master_thesis/reanalysis_data/ERA5/ERA5_pressure_levels


# Regridding

Within each directory each containing a season we want to run the command "cdo remapcon" to create different grid resolutions of the datafiles. We want 10 different resolutions (meaning 9 additional as the original counts) from 0.25 degrees (original) to 3 degrees in. It would make sense to use resolutions that other reanalysis datasets and climate models use so they could be a proxy of these.
After the algorithm each directory should contain 10 .nc files with names describing what season it is and what grid resolution they are.


# Testing fase
I want to be able to control the functionality by initially testing the script on a smaller subset by being able to choose the timerange when downloading.

# Checking progress
I want some checks to see if everything runs correctly. 
1. Print statement for each begun download request
2. Print statement for every succesful download
3. print for every succesful merging of datasets, where it is checked that the first timestep is december 1st and last timestep is febuary of the following year (which can have different dates because of leap year).
4. Print statement for succesful regridding for each new grid made.
5. Done statement


# functionality and speed
You are allowed to suggest methods to optimize the downloading and regridding functionality.