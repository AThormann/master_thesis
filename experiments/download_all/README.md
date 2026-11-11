# ERA5 DJF seasons: download, regrid, track, post-process, sample

Scripts in this folder:

| Step | Script | Does |
|---|---|---|
| 1 | `ERA5_download_seasons.py` | downloads ERA5 (6-hourly) 850 hPa `u`, `v`, `vo` (`pl`) and single levels `u10`, `v10`, `msl` (`sl`), builds one file per DJF season, then regrids the 850 hPa files (calls step 2) |
| 2 | `regrid_era5.py` | `regrid_season()`: conservative regridding (`cdo remapcon`) to the resolutions in `RESOLUTIONS`. Also runs on its own |
| 3 | `track_on_seasons.py` | Hodges tracker (`processing.hodges_tracker`) on every season × resolution → raw `.trackjson` |
| 4 | `postprocess_on_seasons.py` | clips to `DOMAIN` and filters with `processing.postprocess_tracks` + `POSTPROCESS_CONFIG` → filtered `.trackjson` |
| 5 | `sample_on_seasons.py` | samples fields along every point of the filtered tracks (`statistics_on_tracks.py` + `SAMPLE_CONFIG`) → points and summary `.csv` |

Every script skips outputs that already exist, so after a crash or interruption you just run the same command again.

## Output layout

```
<out_root>/                                    data (external drive)
├── raw/                                       yearly downloads (temporary, deleted after merging)
├── _regrid_weights/                           cdo grid files + remapping weights (sl/ for single levels)
├── 19401941_DJF/
│   ├── 19401941_DJF_0p25.nc                   850 hPa, original ERA5 grid (Dec 1940 + Jan/Feb 1941)
│   ├── 19401941_DJF_1p0.nc                    850 hPa, one file per resolution in RESOLUTIONS
│   ├── ...
│   ├── 19401941_DJF_sl_0p25.nc                single levels, limited area SL_AREA
│   └── 19401941_DJF_sl_1p0.nc                 single levels regridded (only with --regrid-sl / --sample-res)
└── ...

<track_dir>/                                   tracks + statistics (e.g. master_thesis/trackjson/seasons)
├── 19401941_DJF/
│   ├── 19401941_DJF_0p25_vo_raw.trackjson
│   ├── 19401941_DJF_0p25_vo_filtered.trackjson
│   ├── 19401941_DJF_0p25_vo_points_sl0p25.csv     one row per track point, fields sampled on the 0.25° grid
│   ├── 19401941_DJF_0p25_vo_summary_sl0p25.csv    one row per track
│   └── ...                                        the same for every track resolution
└── ...
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
4. Accept the licence on both dataset pages (`reanalysis-era5-pressure-levels` and `reanalysis-era5-single-levels` → Download tab → Terms of use).

**Disk space:** the 850 hPa files at all resolutions take a few GB per season. Point `--out-root` at an external drive.

## Settings in `algorithms/configs.py`

| Name | Used by | Meaning |
|---|---|---|
| `TRACK_CONFIGS`, `track_variable`, `object_threshold` | step 3 | tracker settings |
| `REGION`, `DOMAIN`, `POSTPROCESS_CONFIG` | step 4 | clipping and filtering |
| `SL_AREA` | steps 1, 2, 5 | `[north, west, south, east]` of the single-level download. Must contain `DOMAIN` plus `radius_km` |
| `SAMPLE_CONFIG` | step 5 | `radius_km`, `samples` (from the single levels), `pl_samples` (from the 0.25° 850 hPa file), `sample_res` |

## 1. `ERA5_download_seasons.py`

| Argument | Default | Meaning |
|---|---|---|
| `--out-root` | `reanalysis_data/ERA5/ERA5_pressure_levels` | folder for the season folders (use the external drive) |
| `--first-season` | `1940` | first season (its December year) |
| `--last-season` | `2025` | last season (2025 → Dec 2025 + Jan/Feb 2026) |
| `--days` | `1 … 31` | days of each month to download. Use a few for testing |
| `--datasets` | `pl sl` | `pl`: 850 hPa `u`, `v`, `vo` (global); `sl`: `u10`, `v10`, `msl` (area `SL_AREA`) |
| `--workers` | `4` | CDS requests running at the same time |
| `--threads` | `4` | `cdo -P` threads for regridding |
| `--regrid` / `--no-regrid` | `--regrid` | regrid the 850 hPa files right after merging |
| `--regrid-sl` | off | also regrid the single-level files (cut back to `SL_AREA`) |
| `--keep-raw` | off | keep the yearly raw files instead of deleting them |

**What it does:**
- One request per dataset and year, containing only the months a season still needs: the first year gets only December, the last year only Jan/Feb. A dataset whose season file already exists is not downloaded again, so `--datasets sl` adds the single levels to existing seasons.
- Merges Dec `S` with Jan/Feb `S+1`.
- Checks the merge: the first step must be Dec 1 00:00, the last must be Feb 28/29 18:00, and no step may be missing.
- Deletes raw years once they are no longer needed.
- Precipitation is not downloaded. At 00/06/12/18 ERA5 `total_precipitation` is only a 1-hour accumulation (6-hour totals need all 24 hours). Mixing it with instantaneous fields also makes the CDS return a zip of two files.

**Progress prints:** `[download] START / OK <dataset> year`, `[merge] OK <dataset>`, `[cleanup]`, `[regrid] OK` and `[done]`.

```bash
# small test: two seasons, two days per month
python ERA5_download_seasons.py --first-season 1998 --last-season 1999 --days 1 2 --out-root /tmp/era5_test

# add the single levels to seasons that already have the 850 hPa files
python ERA5_download_seasons.py --datasets sl --out-root /Volumes/<drive>/ERA5_pressure_levels

# everything, in the background with a log
nohup python ERA5_download_seasons.py --out-root /Volumes/<drive>/ERA5_pressure_levels > download.log 2>&1 &
```

## 2. `regrid_era5.py`

Runs automatically from step 1 (and from step 5 for `--sample-res`). Run it yourself only to regrid again or to add resolutions.

| Argument | Default | Meaning |
|---|---|---|
| `out_root` (positional) | – | folder containing the season folders |
| `--seasons` | all season folders | season start years, e.g. `--seasons 1998 1999` |
| `--threads` | `4` | `cdo -P` threads |
| `--kind` | `pl` | `pl`: 850 hPa files (global grids); `sl`: single-level files (same grids, cut back to `SL_AREA`) |

The resolutions live in the `RESOLUTIONS` dictionary at the top of the file; commented-out entries are skipped. All grids are regular lon-lat grids that include the poles, like ERA5.

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
python regrid_era5.py /Volumes/<drive>/ERA5_pressure_levels --kind sl
```

## 3. `track_on_seasons.py`

Runs the Hodges tracker on every season and resolution and writes `<track_dir>/<S><S+1>_DJF/<S><S+1>_DJF_<res>_vo_raw.trackjson`.

| Argument | Default | Meaning |
|---|---|---|
| `--out-root` | `reanalysis_data/ERA5/ERA5_pressure_levels` | folder containing the season folders |
| `--seasons` | all season folders | season start years |
| `--resolutions` | all in `RESOLUTIONS` | resolutions in degrees, e.g. `--resolutions 3.0 1.875` |
| `--track-dir` | `<out-root>/trackjson` | folder for the season folders with tracks and statistics |

**Steps per file:**
1. `processing.prepare_tracking_field(file, variable=track_variable, level=850, domain=None)`. The whole globe is tracked, because the spectral filter needs the full sphere.
2. `processing.hodges_tracker(ds, output_file_name=None, track_var=track_variable, obj_threshold=object_threshold, **TRACK_CONFIGS)`. The tracker settings come from `algorithms/configs.py`, so change them there. The file is written by the script, not by `hodges_tracker`'s own output folder.
3. Prints `[track] OK <season> <res>: N tracks (x min)`.

**Run time:** a coarse season takes about half a minute, 0.25° about 6 minutes.

```bash
python track_on_seasons.py --out-root /Volumes/<drive>/ERA5_pressure_levels --track-dir ../../trackjson/seasons --seasons 1998 --resolutions 3.0 1.875
nohup python track_on_seasons.py --out-root /Volumes/<drive>/ERA5_pressure_levels --track-dir ../../trackjson/seasons > track.log 2>&1 &
```

## 4. `postprocess_on_seasons.py`

Filters the raw tracks and writes `<track_dir>/<S><S+1>_DJF/<S><S+1>_DJF_<res>_vo_filtered.trackjson`.

| Argument | Default | Meaning |
|---|---|---|
| `--out-root` | `reanalysis_data/ERA5/ERA5_pressure_levels` | folder containing the season folders (used for each season's time range) |
| `--seasons` | all season folders | season start years |
| `--resolutions` | all in `RESOLUTIONS` | resolutions in degrees |
| `--track-dir` | `<out-root>/trackjson` | folder with the season folders of tracks |
| `--overwrite` | off | redo existing filtered files, e.g. after changing `POSTPROCESS_CONFIG` |

**Steps per file** (same as `experiments/validation_tracks/coarse_validation.ipynb`):
1. `pst.load_tracks(raw file)`.
2. **Time range of the season:** the first and last time step of the season's `.nc` file, e.g. `('1998-12-01 00:00', '1999-02-28 18:00')`. This replaces the hard-coded `TIME_RANGE` from the notebooks, and it handles leap years and test runs with `--days`.
3. **Clip to the domain:**
   - `helpers.tracks_to_points`, then `validation.clip_to_domain(df, DOMAIN, time_range)`.
   - The points are then turned back into tracks while keeping **every** track variable (`vo`, `raw_value`, `object_*`). `validation.dataframe_to_tracks` would keep only `vo`.
   - The tracks are global, so without this the points outside `DOMAIN` would count as "at the edge".
4. **Filter:** `processing.postprocess_tracks(tracks, **POSTPROCESS_CONFIG)` with `time_range` set per season. All other settings come from `POSTPROCESS_CONFIG` in `configs.py`. With a real time range, `drop_time_truncated=True` in `configs.py` now also works.
5. **Report:** prints `processing.filter_report(summary)` and `[post] OK <season> <res> (<time range>): raw N -> clipped N -> kept N`.

```bash
python postprocess_on_seasons.py --out-root /Volumes/<drive>/ERA5_pressure_levels --track-dir ../../trackjson/seasons
python postprocess_on_seasons.py --out-root /Volumes/<drive>/ERA5_pressure_levels --track-dir ../../trackjson/seasons --overwrite   # after changing POSTPROCESS_CONFIG
```

## 5. `sample_on_seasons.py`

Samples fields along every point of the filtered tracks with `algorithms/statistics_on_tracks.py` and writes, per season and track resolution:
- `…_vo_points_sl<res>.csv`: one row per track point
- `…_vo_summary_sl<res>.csv`: one row per track

| Argument | Default | Meaning |
|---|---|---|
| `--out-root` | `reanalysis_data/ERA5/ERA5_pressure_levels` | folder containing the season folders with the `.nc` files |
| `--seasons` | all season folders | season start years |
| `--resolutions` | all in `RESOLUTIONS` | track resolutions to sample |
| `--track-dir` | `<out-root>/trackjson` | folder with the season folders of tracks; the CSVs are written there too |
| `--sample-res` | `SAMPLE_CONFIG["sample_res"]` (0.25) | grid of the sampled single-level field. Other values regrid `…_sl_0p25.nc` first (`cdo remapcon`) |
| `--overwrite` | off | redo existing CSV files |

**Columns of the points CSV:**

| Column | Unit | Source |
|---|---|---|
| `track_id`, `time`, `lat`, `lon` | | filtered track |
| `vo` | s⁻¹ | tracked (T6–42 filtered) vorticity |
| `raw_value` | s⁻¹ | stored by the tracker; currently identical to `vo` |
| `object_*` | km, km², ° | size and shape of the detected vorticity object |
| `inside_area` | | False if the 500 km circle is not inside `SL_AREA`. Those points have NaN single-level samples |
| `msl_min`, `msl_center` | hPa | minimum within `radius_km`, value at the track point |
| `ws10_max`, `ws10_mean` | m s⁻¹ | 10 m wind speed (computed on the grid) within `radius_km` |
| `vo850_unfilt_max` | s⁻¹ | maximum unfiltered 850 hPa vorticity within `radius_km` (0.25° file) |
| `ws850_max` | m s⁻¹ | maximum 850 hPa wind speed within `radius_km` (0.25° file) |
| `speed_kmh`, `direction_deg` | km h⁻¹, ° | motion since the previous point (0° = north) |
| `dvo_dt` | s⁻¹ per 6 h | change of `vo` |
| `dmsl_24h` | hPa | change of `msl_min` over the previous 24 h |
| `deepening_bergeron` | Bergeron | 24 h deepening normalised to 60°N; ≥ 1 is an explosive cyclone |

**Checks** (an error stops that file; warnings are printed with `[sample] WARNING`):
- Every track time must be a time step of the field files. `sample_tracks` would otherwise silently take the nearest time.
- Points outside `SL_AREA` (including `radius_km`) give a warning and NaN, instead of the edge value `sample_tracks` would return.
- No NaN is allowed for points inside the area.
- The row count must equal the number of points in the filtered tracks, with the same `track_id`s.
- Values outside physical ranges give a warning: `msl` 900–1070 hPa, wind 0–90 m s⁻¹, |vo| < 1e-3 s⁻¹.

```bash
python sample_on_seasons.py --out-root /Volumes/<drive>/ERA5_pressure_levels --track-dir ../../trackjson/seasons
python sample_on_seasons.py --out-root /Volumes/<drive>/ERA5_pressure_levels --track-dir ../../trackjson/seasons --sample-res 1.0
```

## Order of a full run

```bash
python ERA5_download_seasons.py   --out-root /Volumes/<drive>/ERA5_pressure_levels
python track_on_seasons.py        --out-root /Volumes/<drive>/ERA5_pressure_levels --track-dir ../../trackjson/seasons
python postprocess_on_seasons.py  --out-root /Volumes/<drive>/ERA5_pressure_levels --track-dir ../../trackjson/seasons
python sample_on_seasons.py       --out-root /Volumes/<drive>/ERA5_pressure_levels --track-dir ../../trackjson/seasons
```

Load the results in a notebook with `pst.load_tracks(".../19981999_DJF/19981999_DJF_1p875_vo_filtered.trackjson")` and `pd.read_csv(".../19981999_DJF/19981999_DJF_1p875_vo_points_sl0p25.csv", parse_dates=["time"])`.
