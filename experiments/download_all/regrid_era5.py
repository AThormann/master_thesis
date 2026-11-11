"""
Conservative regridding of the merged ERA5 DJF season files with CDO.

Every season folder <S><S+1>_DJF/ contains <S><S+1>_DJF_0p25.nc (pressure levels, kind="pl")
and optionally <S><S+1>_DJF_sl_0p25.nc (single levels, kind="sl"), both on the original ERA5 grid.
regrid_season() adds one file per resolution in RESOLUTIONS for the chosen kind.
Single-level files cover a limited area: they are regridded to the same global grids and
then cut back to that area (sellonlatbox), so the grid points match the pressure-level grids.

The weights for each target grid are computed once (cdo gencon) and reused for
every season (cdo remap), which gives the same result as cdo remapcon.

Usage (standalone):
    python regrid_era5.py <out_root>                 # all season folders in out_root
    python regrid_era5.py <out_root> --seasons 1998  # only 19981999_DJF
    python regrid_era5.py <out_root> --kind sl       # single-level files instead of 850 hPa
"""

import argparse
import shutil
import subprocess
import sys
from pathlib import Path

# cdo on PATH, otherwise the one installed next to this python (env not activated)
CDO = shutil.which("cdo") or str(Path(sys.executable).parent / "cdo")


# Resolution (degrees) -> reanalysis / model the grid is a proxy for
RESOLUTIONS = {
    0.25:  "ERA5 native (original, not regridded)",
    #0.5:   "MERRA-2 / CFSR class, HighResMIP",
    #0.75:  "ERA-Interim native",
    1.0:   "CMIP6 HR models (e.g. MPI-ESM1-2-HR), GFS",
    #1.125: "BCC-CSM2-MR (T106)",
    1.25:  "JRA-55",
    #1.5:   "ERA-Interim / ERA-40 distribution grid",
    1.875: "MPI-ESM1-2-LR (T63), NorESM2-LR latitude",
    #2.5:   "NCEP/NCAR R1/R2, NorESM2-LR longitude",
    3.0:   "coarse CMIP-class upper bound (~T42)",
}
NATIVE_RES = 0.25


def res_tag(res):
    """0.25 -> '0p25', 1.0 -> '1p0', 1.875 -> '1p875'"""
    tag = f"{res:g}"
    if "." not in tag:
        tag += ".0"
    return tag.replace(".", "p")


def season_file(season_dir, res, kind="pl"):
    """
    Path of the season file at a given resolution.
    kind="pl": 19981999_DJF/19981999_DJF_1p875.nc (850 hPa), kind="sl": 19981999_DJF/19981999_DJF_sl_1p875.nc
    """
    season_dir = Path(season_dir)
    prefix = "" if kind == "pl" else f"{kind}_"
    return season_dir / f"{season_dir.name}_{prefix}{res_tag(res)}.nc"


def grid_size(res):
    """(nlon, nlat) of a global pole-inclusive regular grid."""
    return round(360 / res), round(180 / res) + 1


def write_grid_file(res, weights_dir):
    """
    CDO grid description of a global regular lon-lat grid.
    Same layout as ERA5: longitudes from 0, latitudes from 90 to -90 (poles included).
    """
    nlon, nlat = grid_size(res)
    path = Path(weights_dir) / f"grid_{res_tag(res)}.txt"
    path.write_text(
        "gridtype = lonlat\n"
        "xname    = longitude\n"
        "yname    = latitude\n"
        f"xsize    = {nlon}\n"
        f"ysize    = {nlat}\n"
        "xfirst   = 0\n"
        f"xinc     = {res}\n"
        "yfirst   = 90\n"
        f"yinc     = {-res}\n"
    )
    return path


def run_cdo(args, quiet=False):
    """Run cdo. quiet=True hides stderr unless cdo fails (gencon prints harmless HDF5 noise)."""
    result = subprocess.run([CDO, "-s", *args], stderr=subprocess.PIPE if quiet else None, text=True)
    if result.returncode != 0:
        if quiet:
            print(result.stderr)
        raise subprocess.CalledProcessError(result.returncode, result.args)


def get_weights(res, source_file, weights_dir):
    """Grid file and conservative remapping weights for one resolution (computed only once)."""
    weights_dir = Path(weights_dir)
    weights_dir.mkdir(parents=True, exist_ok=True)
    grid = write_grid_file(res, weights_dir)
    weights = weights_dir / f"weights_con_{res_tag(res)}.nc"
    if not weights.exists():
        try:
            run_cdo([f"gencon,{grid}", str(source_file), str(weights)], quiet=True)
        except BaseException:
            weights.unlink(missing_ok=True)
            raise
    return grid, weights


def regrid_season(season_dir, resolutions=RESOLUTIONS, weights_dir=None, threads=4, kind="pl", area=None):
    """
    Regrid the 0.25 degree season file of one kind ("pl" or "sl") in season_dir to every resolution.

    area : (north, west, south, east) to cut the result back to (used for the limited-area
           single-level files, e.g. configs.SL_AREA). None keeps the global grid.
    Existing output files are skipped. Returns the list of this kind's files in season_dir.
    """
    season_dir = Path(season_dir)
    source = season_file(season_dir, NATIVE_RES, kind)
    if not source.exists():
        raise FileNotFoundError(f"Missing original season file: {source}")
    if weights_dir is None:
        weights_dir = season_dir.parent / "_regrid_weights" / ("" if kind == "pl" else kind)
    cut = []
    if area is not None:
        north, west, south, east = area
        cut = [f"-sellonlatbox,{west},{east},{south},{north}"]

    for res in resolutions:
        if res == NATIVE_RES:
            continue
        out = season_file(season_dir, res, kind)
        nlon, nlat = grid_size(res)
        if out.exists():
            print(f"[regrid] SKIP {out.name} (exists)")
            continue
        grid, weights = get_weights(res, source, weights_dir)
        tmp = out.with_suffix(".part.nc")
        run_cdo(["-P", str(threads), "-f", "nc4", "-z", "zip_1",
                 *cut, f"-remap,{grid},{weights}", str(source), str(tmp)])
        tmp.rename(out)
        print(f"[regrid] OK {season_dir.name} -> {out.name} ({nlon}x{nlat}{' cut to area' if cut else ''})")

    files = [season_file(season_dir, res, kind) for res in resolutions]
    missing = [f.name for f in files if not f.exists()]
    if missing:
        print(f"[regrid] WARNING {season_dir.name}: missing {missing}")
    return files


def main():
    parser = argparse.ArgumentParser(description="Regrid ERA5 DJF season folders with cdo remapcon.")
    parser.add_argument("out_root", type=Path, help="folder containing the <S><S+1>_DJF season folders")
    parser.add_argument("--seasons", type=int, nargs="+", help="season start years (default: all folders)")
    parser.add_argument("--threads", type=int, default=4, help="cdo -P threads")
    parser.add_argument("--kind", choices=["pl", "sl"], default="pl",
                        help="pl: 850 hPa files (global), sl: single-level files (cut to configs.SL_AREA)")
    args = parser.parse_args()

    if args.seasons:
        dirs = [args.out_root / f"{s}{s + 1}_DJF" for s in args.seasons]
    else:
        dirs = sorted(d for d in args.out_root.glob("*_DJF") if d.is_dir())

    area = None
    if args.kind == "sl":
        sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "algorithms"))
        from configs import SL_AREA as area
    for d in dirs:
        regrid_season(d, threads=args.threads, kind=args.kind, area=area)
    print(f"[done] regridded {len(dirs)} seasons in {args.out_root}")


if __name__ == "__main__":
    main()
