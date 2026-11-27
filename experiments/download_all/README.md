# ERA5 seasons (DJF, MAM, JJA, SON): download, regrid, track, post-process, sample

Scripts in this folder:

| Step | Script | Does |
|---|---|---|
| 1 | `ERA5_download_seasons.py` | downloads ERA5 (6-hourly, monthly requests) 850 hPa `u`, `v`, `vo` (`pl`) and single levels `u10`, `v10`, `msl`, gust (`sl`), builds one file per season with a 10-day buffer on both sides, then regrids the 850 hPa files (calls step 2) |
| 2 | `regrid_era5.py` | `regrid_season()`: conservative regridding (`cdo remapcon`) to the resolutions in `RESOLUTIONS`. Also runs on its own |
| 3 | `track_on_seasons.py` | Hodges tracker (`processing.hodges_tracker`) on every season × resolution → raw `.trackjson` |
| 4 | `postprocess_on_seasons.py` | filters the full-length tracks with `processing.postprocess_tracks` + `POSTPROCESS_CONFIG`, keeps the tracks whose peak vorticity is in the season → filtered `.trackjson` |
| 5 | `sample_on_seasons.py` | samples fields along every point of the filtered tracks (`statistics_on_tracks.py` + `SAMPLE_CONFIG`) → points and summary `.csv` |
| – | `seasons.py` | shared helpers: season folder names, season / buffered time ranges, months needed |

Every script skips outputs that already exist, so after a crash or interruption you just run the same command again.

## Output layout

```
<out_root>/                                    data (external drive)
├── raw/                                       monthly downloads (temporary, deleted after merging)
├── _regrid_weights/                           cdo grid files + remapping weights (sl/ for single levels)
├── 19401941_DJF/
│   ├── 19401941_DJF_0p25.nc                   850 hPa, original ERA5 grid (21 Nov 1940 - 10 Mar 1941: DJF + 10-day buffer)
│   ├── 19401941_DJF_1p0.nc                    850 hPa, one file per resolution in RESOLUTIONS
│   ├── ...
│   ├── 19401941_DJF_sl_0p25.nc                single levels, limited area SL_AREA
│   └── 19401941_DJF_sl_1p0.nc                 single levels regridded (only with --regrid-sl / --sample-res)
├── 1941_MAM/ ...                               same layout for MAM, JJA, SON
└── ...

<track_dir>/                                   tracks + statistics (e.g. master_thesis/trackjson/seasons)
├── 19401941_DJF/
│   ├── 19401941_DJF_0p25_vo_raw.trackjson
│   ├── 19401941_DJF_0p25_vo_filtered.trackjson
│   ├── 19401941_DJF_0p25_vo_points_sl0p25.csv     one row per track point, fields sampled on the 0.25° grid
│   ├── 19401941_DJF_0p25_vo_summary_sl0p25.csv    one row per track
│   └── ...                                        the same for every track resolution
├── 1941_MAM/ ...                               same layout for MAM, JJA, SON
└── ...
```

## Seasons and the 10-day buffer

| Type | Months | Folder (year `S`) |
|---|---|---|
| DJF | Dec `S` + Jan/Feb `S+1` | `<S><S+1>_DJF`, e.g. `19981999_DJF` |
| MAM | Mar–May `S` | `<S>_MAM`, e.g. `1999_MAM` |
| JJA | Jun–Aug `S` | `<S>_JJA` |
| SON | Sep–Nov `S` | `<S>_SON` |

The other season types use the same file layout, e.g. `1999_MAM/1999_MAM_0p25.nc` and `1999_MAM/1999_MAM_1p875_vo_summary_sl0p25.csv`. `name.split("_")[2]` is the resolution for every season type.

- **Buffer:** every season file (`pl` and `sl`) covers the season plus `SEASON_BUFFER_DAYS` (10, in `configs.py`) on both sides. For example, `19981999_DJF` runs from 1998-11-21 00:00 to 1999-03-10 18:00. So the tracker follows storms across the season boundary, and every point of a kept track can be sampled.
- **Season of a storm:** a storm belongs to the season in which its peak (tracked, T6–42 filtered) vorticity lies. Season boundaries exclude the buffer. A storm that peaks in the buffer belongs to the previous or next season, and is kept when that season is processed. A storm crossing a boundary is tracked in both season files but kept in only one.
- **Full range:** 1940 to 2025 (DJF 2025 = Dec 2025 + Jan/Feb 2026).
- **Old DJF data:** DJF files built before this change have no buffer. Step 1 warns about them (`[merge] WARNING … Delete the season's files`), and step 4 warns that tracks crossing the boundary may get a wrong peak. Delete the season folder (data and tracks) to rebuild it with the buffer.

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
| `REGION` | step 4 | selection box: a track is kept if its maximum `vo` lies in `REGION` and it has at least 4 points in `REGION` |
| `DOMAIN` | plotting | no longer limits the tracks (no clipping); only used by the edge filter, which is off |
| `POSTPROCESS_CONFIG` | step 4 | filter settings: `region_criterion=("passes", "peak")`, `min_points_in_region=4`, RSPLICE, `max_frac_at_edge=None` (edge filter off) |
| `SL_AREA` | steps 1, 2, 5 | `[north, west, south, east]` of the single-level download. Must contain `DOMAIN` plus `radius_km` |
| `SAMPLE_CONFIG` | step 5 | `radius_km`, `samples` (from the single levels), `pl_samples` (from the 0.25° 850 hPa file), `sample_res` |

## 1. `ERA5_download_seasons.py`

| Argument | Default | Meaning |
|---|---|---|
| `--out-root` | `reanalysis_data/ERA5/ERA5_pressure_levels` | folder for the season folders (use the external drive) |
| `--first-season` | `1940` | first season year (December year for DJF) |
| `--last-season` | `2025` | last season year (DJF 2025 → Dec 2025 + Jan/Feb 2026) |
| `--season-types` | `DJF MAM JJA SON` | season types to build |
| `--days` | `1 … 31` | days of each month to download. Use a few for testing |
| `--datasets` | `pl sl` | `pl`: 850 hPa `u`, `v`, `vo` (global); `sl`: `u10`, `v10`, `msl` (area `SL_AREA`) |
| `--workers` | `4` | CDS requests running at the same time |
| `--threads` | `4` | `cdo -P` threads for regridding |
| `--regrid` / `--no-regrid` | `--regrid` | regrid the 850 hPa files right after merging |
| `--regrid-sl` | off | also regrid the single-level files (cut back to `SL_AREA`) |
| `--keep-raw` | off | keep the monthly raw files instead of deleting them |

**What it does:**
- One request per dataset and month (`raw/era5_pl850_YYYY_MM.nc`, `raw/era5_sl_YYYY_MM.nc`). Each request contains only the days the seasons of the run still need: a buffer month on its own needs only 10 days.
  - Monthly requests follow ERA5's storage on the CDS, and are faster and smaller than one global year (about 18 GB).
  - A raw file that exists is reused only if it holds all requested days.
  - A dataset whose season file already exists is not downloaded again, so `--datasets sl` adds the single levels to existing seasons.
- Merges the months of each season, cut to the buffered range.
- Checks the merge: the time steps must be exactly the 6-hourly steps from (season start − 10 days) to (season end + 10 days), with no step missing or duplicated. February 28/29 is handled.
- Deletes a raw month once every season of the run that uses it is written. For example, May is used by both MAM (all of it) and the JJA buffer.
- Precipitation is not downloaded. At 00/06/12/18 ERA5 `total_precipitation` is only a 1-hour accumulation (6-hour totals need all 24 hours). Mixing it with instantaneous fields also makes the CDS return a zip of two files.

**Progress prints:** `[download] START / OK <dataset> year`, `[merge] OK <dataset>`, `[cleanup]`, `[regrid] OK` and `[done]`.

```bash
# small test: MAM 1998, three days per month (the buffer days 25, 1, 2 fall in Feb / Jun)
python ERA5_download_seasons.py --first-season 1998 --last-season 1998 --season-types MAM --days 1 2 25 --out-root /tmp/era5_test

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
| `--seasons` | all season folders | season years, e.g. `--seasons 1998 1999` |
| `--season-types` | `DJF MAM JJA SON` | season types |
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

Runs the Hodges tracker on every season and resolution, on the buffered season file. Writes `<track_dir>/<season>/<season>_<res>_vo_raw.trackjson`.

| Argument | Default | Meaning |
|---|---|---|
| `--out-root` | `reanalysis_data/ERA5/ERA5_pressure_levels` | folder containing the season folders |
| `--seasons` | all season folders | season years (December year for DJF) |
| `--season-types` | `DJF MAM JJA SON` | season types (also for steps 4 and 5) |
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

Filters the raw tracks and writes `<track_dir>/<season>/<season>_<res>_vo_filtered.trackjson`.

| Argument | Default | Meaning |
|---|---|---|
| `--out-root` | `reanalysis_data/ERA5/ERA5_pressure_levels` | folder containing the season folders (used for each season's time range) |
| `--seasons` | all season folders | season start years |
| `--resolutions` | all in `RESOLUTIONS` | resolutions in degrees |
| `--track-dir` | `<out-root>/trackjson` | folder with the season folders of tracks |
| `--overwrite` | off | redo existing filtered files, e.g. after changing `POSTPROCESS_CONFIG` |

Tracks are **not clipped** to `DOMAIN`: the global tracks are kept at full length. A track is kept only if both region criteria hold (`region_criterion=("passes", "peak")`):
- `passes`: at least 4 points (24 h) inside `REGION`,
- `peak`: the point of maximum `vo` (the tracked T6–42 filtered vorticity, over the whole track) lies inside `REGION`.

A track is also kept only if its peak vorticity lies in the season itself, not in the buffer (`season`, column `ok_season`, `postprocess_tracks(season_range=...)`).

The edge filter is off (`max_frac_at_edge=None`): tracking is global, so there is no domain edge. The other filters are unchanged: RSPLICE (≥ 8 points / 48 h, ≥ 10° endpoint displacement), `min_peak=None`, `drop_time_truncated=False`.

**Steps per file:**
1. `pst.load_tracks(raw file)`.
2. **Time ranges:**
   - `time_range` is the first and last time step of the season's `.nc` file, including the buffer, e.g. `('1998-11-21 00:00', '1999-03-10 18:00')`. It is used for the truncation flags.
   - `season_range` is the season without the buffer, taken from the folder name (`seasons.core_range`), e.g. `('1998-12-01 00:00', '1999-02-28 18:00')`.
   - A warning is printed if the file does not cover the buffered range (old DJF files).
3. **Filter:** `processing.postprocess_tracks(tracks, **POSTPROCESS_CONFIG)` with `time_range` and `season_range` set per season. All other settings come from `POSTPROCESS_CONFIG` in `configs.py`. With a real time range, `drop_time_truncated=True` in `configs.py` now also works. (Clipping to `DOMAIN` is kept as commented-out code in `postprocess_season`, so it can be switched back on.)
4. **Report:** prints `processing.filter_report(summary)`, with one row per criterion (`region_passes`, `region_peak`, `edge`, `peak`, `time`, `rsplice`, `season`), and `[post] OK <season> <res> (<time range>): raw N -> kept N -> <file>`.

**Known consequences:**
- **Fewer tracks than with `passes` alone.** Tracks with ≥ 4 points in `REGION` that peak outside it are removed, and so are tracks that peak in `REGION` but spend less than 24 h there. The filter report shows how many tracks each criterion removes, and how many it removes on its own.
- **Truncated tracks.** With the 10-day buffer, a track is cut off only at the edge of the buffer. Its peak can then be wrong only if the storm lives more than 10 days across the season boundary.
- **More NaN single-level samples in step 5.** Full-length tracks often leave `SL_AREA`, so more points get `inside_area=False` and the `[sample] WARNING` appears more often. This is expected. The peak point is always in `REGION`, so values at peak time are almost always sampled. The exception is the far north-east corner of `REGION` (north of about 77°N and east of about 15°E): there 500 km spans more than the 20° of longitude between `REGION` and the east edge of `SL_AREA` (40°E). The check that points inside the area have no NaN stays as it is. Per-track summaries of single-level fields (e.g. `msl_min_min`) only cover the part of the track inside `SL_AREA`. The 850 hPa samples come from the global 0.25° file and cover the whole track.

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
| `ws10_max_lat`, `ws10_max_lon` | ° | grid point of the 10 m wind maximum within `radius_km` (outputs in `SAMPLE_CONFIG["locate"]`; NaN outside the area; longitude −180…180) |

**Columns of the summary CSV** (one row per track), in addition to lifetime, genesis/lysis and the per-column extremes (`<col>_max`, `msl…_min`):

| Column | Meaning |
|---|---|
| `vo_peak_time`, `vo_peak_lat`, `vo_peak_lon` | track point of the peak tracked vorticity (the point that decides the season) |
| `ws10_max_time`, `ws10_max_track_lat`, `ws10_max_track_lon` | track centre at the time of the track's maximum `ws10_max` |
| `ws10_max_lat`, `ws10_max_lon` | grid point of that wind maximum (NaN if the track has no single-level samples) |

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
