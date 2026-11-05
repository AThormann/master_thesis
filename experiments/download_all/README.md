# ERA5 DJF seasons: download, regrid, track, post-process

Scripts in this folder:

| Step | Script | Does |
|---|---|---|
| 1 | `ERA5_download_seasons.py` | downloads ERA5 850 hPa `u`, `v`, `vo` (6-hourly), builds one file per DJF season, then regrids (calls step 2) |
| 2 | `regrid_era5.py` | `regrid_season()`: conservative regridding (`cdo remapcon`) to 10 resolutions. Also runs on its own |
| 3 | `track_on_seasons.py` | Hodges tracker (`processing.hodges_tracker`) on every season × resolution → raw `.trackjson` |
| 4 | `postprocess_on_seasons.py` | clips to `DOMAIN` and filters with `processing.postprocess_tracks` + `POSTPROCESS_CONFIG` → filtered `.trackjson` |

Every script skips outputs that already exist, so after a crash or interruption you just run the same command again.

## Output layout

```
<out_root>/
├── raw/                                  yearly downloads (temporary, deleted after merging)
├── _regrid_weights/                      cdo grid files + remapping weights (computed once)
├── 19401941_DJF/
│   ├── 19401941_DJF_0p25.nc              original ERA5 grid (Dec 1940 + Jan/Feb 1941)
│   ├── 19401941_DJF_0p5.nc
│   ├── ...                               10 files in total, see the resolution table
│   └── 19401941_DJF_3p0.nc
├── ...
├── 20252026_DJF/
└── trackjson/
    ├── 19401941_DJF_0p25_vo_raw.trackjson
    ├── 19401941_DJF_0p25_vo_filtered.trackjson
    └── ...                               raw + filtered for every season and resolution
```

Season `S` is December of year `S` plus January and February of year `S+1`. The full range is 1940/41 to 2025/26 (86 seasons).

## One-time setup

1. Use the `master_env` environment (`mamba activate master_env`). It has xarray, netCDF4, cdo and pystormtracker.
2. Install the CDS client:
   ```bash
   mamba install -n master_env -c conda-forge cdsapi
   ```
3. Log in at https://cds.climate.copernicus.eu, copy your API token from your profile page and create `~/.cdsapirc`:
   ```
   url: https://cds.climate.copernicus.eu/api
   key: <your-token>
   ```
4. Accept the licence on the dataset page (`reanalysis-era5-pressure-levels` → Download tab → Terms of use).

**Disk space:** all seasons at all resolutions take roughly 250–300 GB. Point `--out-root` at an external drive.

## 1. `ERA5_download_seasons.py`

| Argument | Default | Meaning |
|---|---|---|
| `--out-root` | `reanalysis_data/ERA5/ERA5_pressure_levels` | folder for the season folders (use the external drive) |
| `--first-season` | `1940` | first season (its December year) |
| `--last-season` | `2025` | last season (2025 → Dec 2025 + Jan/Feb 2026) |
| `--days` | `1 … 31` | days of each month to download. Use a few for testing |
| `--workers` | `4` | CDS requests running at the same time |
| `--threads` | `4` | `cdo -P` threads for regridding |
| `--regrid` / `--no-regrid` | `--regrid` | regrid every season right after merging |
| `--keep-raw` | off | keep the yearly raw files instead of deleting them |

**What it does:**
- One request per year, containing only the months a season still needs: the first year gets only December, the last year only Jan/Feb.
- Merges Dec `S` with Jan/Feb `S+1`.
- Checks the merge: the first step must be Dec 1 00:00, the last must be Feb 28/29 18:00, and no step may be missing.
- Deletes raw years once they are no longer needed.

**Progress prints:** `[download] START / OK`, `[merge] OK`, `[cleanup]`, `[regrid] OK` and `[done]`.

```bash
# small test: two seasons, two days per month
python ERA5_download_seasons.py --first-season 1998 --last-season 1999 --days 1 2 --out-root /tmp/era5_test

# everything, in the background with a log
nohup python ERA5_download_seasons.py --out-root /Volumes/<drive>/ERA5_pressure_levels > download.log 2>&1 &
```

## 2. `regrid_era5.py`

Runs automatically from step 1. Run it yourself only to regrid again or to add resolutions.

| Argument | Default | Meaning |
|---|---|---|
| `out_root` (positional) | – | folder containing the season folders |
| `--seasons` | all season folders | season start years, e.g. `--seasons 1998 1999` |
| `--threads` | `4` | `cdo -P` threads |

The resolutions live in the `RESOLUTIONS` dictionary at the top of the file. All grids are global regular lon-lat grids that include the poles, like ERA5.

| Resolution | File tag | Grid (lon × lat) | Proxy for |
|---|---|---|---|
| 0.25° | `0p25` | 1440 × 721 | ERA5 native (original) |
| 0.5° | `0p5` | 720 × 361 | MERRA-2 / CFSR class, HighResMIP |
| 0.75° | `0p75` | 480 × 241 | ERA-Interim native |
| 1.0° | `1p0` | 360 × 181 | CMIP6 HR models (e.g. MPI-ESM1-2-HR), GFS |
| 1.125° | `1p125` | 320 × 161 | BCC-CSM2-MR (T106) |
| 1.25° | `1p25` | 288 × 145 | JRA-55 |
| 1.5° | `1p5` | 240 × 121 | ERA-Interim / ERA-40 distribution grid |
| 1.875° | `1p875` | 192 × 97 | MPI-ESM1-2-LR (T63), NorESM2-LR latitude |
| 2.5° | `2p5` | 144 × 73 | NCEP/NCAR R1/R2, NorESM2-LR longitude |
| 3.0° | `3p0` | 120 × 61 | coarse CMIP-class upper bound (~T42) |

The remapping weights are computed once per resolution (`cdo gencon`) and reused for every season. This gives the same result as `cdo remapcon`, only faster.

```bash
python regrid_era5.py /Volumes/<drive>/ERA5_pressure_levels --seasons 1998
```

## 3. `track_on_seasons.py`

Runs the Hodges tracker on every season and resolution and writes `<out_root>/trackjson/<S><S+1>_DJF_<res>_vo_raw.trackjson`.

| Argument | Default | Meaning |
|---|---|---|
| `--out-root` | `reanalysis_data/ERA5/ERA5_pressure_levels` | folder containing the season folders |
| `--seasons` | all season folders | season start years |
| `--resolutions` | all 10 | resolutions in degrees, e.g. `--resolutions 3.0 1.875` |
| `--track-dir` | `<out-root>/trackjson` | folder for the `.trackjson` files |

**Steps per file:**
1. `processing.prepare_tracking_field(file, variable=track_variable, level=850, domain=None)`. The whole globe is tracked, because the spectral filter needs the full sphere.
2. `processing.hodges_tracker(ds, output_file_name=None, track_var=track_variable, obj_threshold=object_threshold, **TRACK_CONFIGS)`. The tracker settings come from `algorithms/configs.py`, so change them there. The file is written by the script, not by `hodges_tracker`'s own output folder.
3. Prints `[track] OK <season> <res>: N tracks (x min)`.

**Run time:** a 3° or 1.875° season takes about half a minute. 0.25° takes much longer, so test on the coarse grids first.

```bash
python track_on_seasons.py --out-root /Volumes/<drive>/ERA5_pressure_levels --seasons 1998 --resolutions 3.0 1.875
nohup python track_on_seasons.py --out-root /Volumes/<drive>/ERA5_pressure_levels > track.log 2>&1 &
```

## 4. `postprocess_on_seasons.py`

Filters the raw tracks and writes `<out_root>/trackjson/<S><S+1>_DJF_<res>_vo_filtered.trackjson`.

| Argument | Default | Meaning |
|---|---|---|
| `--out-root` | `reanalysis_data/ERA5/ERA5_pressure_levels` | folder containing the season folders and `trackjson/` |
| `--seasons` | all season folders | season start years |
| `--resolutions` | all 10 | resolutions in degrees |
| `--track-dir` | `<out-root>/trackjson` | folder with the raw tracks; the filtered tracks are written here too |
| `--overwrite` | off | redo existing filtered files, e.g. after changing `POSTPROCESS_CONFIG` |

**Steps per file** (same as `experiments/validation_tracks/coarse_validation.ipynb`):
1. `pst.load_tracks(raw file)`.
2. **Time range of the season:** the first and last time step of the season's `.nc` file, e.g. `('1998-12-01 00:00', '1999-02-28 18:00')`. This replaces the hard-coded `TIME_RANGE` from the notebooks, and it handles leap years and test runs with `--days`.
3. **Clip to the domain:** `helpers.tracks_to_points`, then `validation.clip_to_domain(df, DOMAIN, time_range)`, then `validation.dataframe_to_tracks(df, variable="vo", scale=1)`. The tracks are global, so without this the points outside `DOMAIN` would count as "at the edge".
4. **Filter:** `processing.postprocess_tracks(tracks, **POSTPROCESS_CONFIG)` with `time_range` set per season. All other settings come from `POSTPROCESS_CONFIG` in `configs.py`. With a real time range, `drop_time_truncated=True` in `configs.py` now also works.
5. **Report:** prints `processing.filter_report(summary)` and `[post] OK <season> <res> (<time range>): raw N -> clipped N -> kept N`.

```bash
python postprocess_on_seasons.py --out-root /Volumes/<drive>/ERA5_pressure_levels --seasons 1998 --resolutions 3.0 1.875
python postprocess_on_seasons.py --out-root /Volumes/<drive>/ERA5_pressure_levels --overwrite   # after changing POSTPROCESS_CONFIG
```

## Order of a full run

```bash
python ERA5_download_seasons.py   --out-root /Volumes/<drive>/ERA5_pressure_levels
python track_on_seasons.py        --out-root /Volumes/<drive>/ERA5_pressure_levels
python postprocess_on_seasons.py  --out-root /Volumes/<drive>/ERA5_pressure_levels
```

Load the results in a notebook with `pst.load_tracks("<out_root>/trackjson/19981999_DJF_1p875_vo_filtered.trackjson")`.
