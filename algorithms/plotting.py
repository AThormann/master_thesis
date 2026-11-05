"""
plotting.py — one consistent figure style for distributions and maps.

Usage
-----
    from plotting import setup_plotting, plot_style, map_axes, add_colorbar

    colors = setup_plotting()            # apply globally, once
    with plot_style(fontsize=12):        # or temporarily, for one figure
        ...
"""
from contextlib import contextmanager

import matplotlib as mpl
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from mpl_toolkits.axes_grid1 import make_axes_locatable

PALETTE = sns.color_palette("colorblind")  # 8 soft, distinguishable colours
#PALETTE = PALETTE.pop(-2)  # remove yellow, which is hard to see on white


def _rc(ticksize=14, fontsize=20, figsize=(16, 9), dpi=150, facecolor="#e8e8e8"):
    """Build the rcParams dict. Everything stylize_ticks() did lives here now."""
    rc = {
        # Figure / saving
        "figure.figsize": figsize,
        "figure.dpi": dpi,
        "savefig.dpi": dpi,
        "savefig.bbox": "tight",
        # Fonts
        "font.family": "serif",
        "font.serif": ["DejaVu Serif"],
        "mathtext.fontset": "stix",
        "font.size": fontsize,
        "legend.fontsize": ticksize,
        # Axes
        "axes.facecolor": facecolor,
        "axes.linewidth": 2.25,
        "axes.formatter.limits": (-1, 3),
        "axes.grid": True,
        "axes.axisbelow": True,  # grid behind data
        "axes.prop_cycle": mpl.cycler(color=PALETTE),
        # Grid
        "grid.linestyle": "--",
        "grid.linewidth": 1.5,
        # Legend
        "legend.facecolor": "white",
        # Ticks on all four sides
        "xtick.top": True,
        "ytick.right": True,
    }
    for ax in ("x", "y"):
        rc.update({
            f"{ax}tick.labelsize": ticksize,
            f"{ax}tick.direction": "in",
            f"{ax}tick.minor.visible": True,
            f"{ax}tick.major.size": 12,
            f"{ax}tick.major.width": 2.5,
            f"{ax}tick.minor.size": 6,
            f"{ax}tick.minor.width": 1.5,
            # default pad is 3.5 pt; +3 pt replaces the old ScaledTranslation offset
            f"{ax}tick.major.pad": 6.5,
        })
    return rc


def setup_plotting(show_palette=False, **kwargs):
    """Apply the style globally. Returns the colour palette."""
    plt.rcParams.update(_rc(**kwargs))
    if show_palette:
        fig, ax = plt.subplots(figsize=(len(PALETTE), 1))
        ax.imshow([list(range(len(PALETTE)))], cmap=mpl.colors.ListedColormap(PALETTE))
        ax.set_axis_off()
        plt.show()
    return PALETTE


@contextmanager
def plot_style(**kwargs):
    """Apply the style only inside a `with` block."""
    with mpl.rc_context(_rc(**kwargs)):
        yield PALETTE


def map_axes(ax, xlabel="Longitude", ylabel="Latitude", grid=False):
    """Map-friendly axes: equal aspect, no dashed grid over the field."""
    ax.set_aspect("equal")
    ax.grid(grid)
    ax.set(xlabel=xlabel, ylabel=ylabel)
    return ax


def add_colorbar(mappable, ax, label=None, location="right", size="4%", pad=None,
                 labelsize=None, ticksize=None, **kwargs):
    """
    Add a colorbar.

    ax        : one axes -> colorbar matched exactly to that axes.
                list/array of axes -> ONE shared colorbar spanning all of them.
    location  : "right", "left", "bottom" or "top".
    pad       : gap between axes and colorbar (inches for a shared colorbar).
    labelsize : font size of the colorbar label (default: rcParams font.size).
    ticksize  : font size of the colorbar tick labels (default: rcParams tick size).
    Works with equal-aspect axes and cartopy GeoAxes.
    """
    horizontal = location in ("top", "bottom")
    shared = np.ndim(ax) > 0

    if shared:
        # Place one colorbar by hand next to the block of axes. This avoids
        # layout engines, which can collapse cartopy axes to zero size.
        axes = list(np.ravel(ax))
        fig = axes[0].figure
        shrink = kwargs.pop("shrink", 0.8)
        thick = kwargs.pop("thickness", 0.18)               # inches
        if pad is None:
            pad = 0.55 if location == "bottom" else 0.2    # inches
        for a in axes:
            a.apply_aspect()                                # real (aspect-locked) positions
        boxes = [a.get_position() for a in axes]
        x0, x1 = min(b.x0 for b in boxes), max(b.x1 for b in boxes)
        y0, y1 = min(b.y0 for b in boxes), max(b.y1 for b in boxes)
        W, H = fig.get_size_inches()
        if horizontal:
            w = shrink * (x1 - x0)
            y = y0 - (pad + thick) / H if location == "bottom" else y1 + pad / H
            rect = [x0 + (x1 - x0 - w) / 2, y, w, thick / H]
        else:
            h = shrink * (y1 - y0)
            x = x1 + pad / W if location == "right" else x0 - (pad + thick) / W
            rect = [x, y0 + (y1 - y0 - h) / 2, thick / W, h]
        cax = fig.add_axes(rect)
        cbar = fig.colorbar(mappable, cax=cax, location=location, **kwargs)
    else:
        # axes_class=plt.Axes: otherwise the divider tries to clone a GeoAxes
        # (which needs a projection) and crashes.
        cax = make_axes_locatable(ax).append_axes(
            location, size=size, pad=0.15 if pad is None else pad, axes_class=plt.Axes)
        cbar = ax.figure.colorbar(mappable, cax=cax, location=location, **kwargs)

    cbar.outline.set_linewidth(plt.rcParams["axes.linewidth"])
    if horizontal:
        cbar.ax.tick_params(which="both", top=False, bottom=(location == "bottom"))
    else:
        cbar.ax.tick_params(which="both", left=False, right=(location == "right"))
    if ticksize is not None:
        cbar.ax.tick_params(labelsize=ticksize)
    if label:
        cbar.set_label(label, size=labelsize)
    return cbar


def _style_geo_axes(ax, extent, coastlines, gridlines, resolution):
    import cartopy.crs as ccrs

    if extent is not None:
        ax.set_extent(extent, crs=ccrs.PlateCarree())
    else:
        ax.set_global()
    if coastlines:
        ax.coastlines(resolution=resolution, linewidth=1.0, color="k")
    if gridlines:
        gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                          linestyle=plt.rcParams["grid.linestyle"],
                          linewidth=0.8, color="grey", alpha=0.7)
        # Only label the outer edges of a grid of maps
        ss = ax.get_subplotspec()
        gl.top_labels = gl.right_labels = False
        gl.bottom_labels = ss is None or ss.is_last_row()
        gl.left_labels = ss is None or ss.is_first_col()
        gl.xlabel_style = gl.ylabel_style = {"size": plt.rcParams["xtick.labelsize"]}
        _fix_tightbbox(ax, gl)
        # Titles are normally placed above the axes' contents; the broken label
        # boxes would push them to y=inf (invisible), so pin them just above the map.
        ax._autotitlepos = False
    ax.spines["geo"].set_linewidth(plt.rcParams["axes.linewidth"])


def _fix_tightbbox(ax, gl):
    """
    Work around a cartopy < 0.26 bug: with gridline labels, a GeoAxes can
    report a NaN/inf bounding box. Jupyter (and savefig with bbox="tight")
    then crops the map away and only the colorbar is left.
    """
    from matplotlib.transforms import Bbox

    original = ax.get_tightbbox

    def safe_tightbbox(renderer=None, *args, **kwargs):
        bb = original(renderer, *args, **kwargs)
        if bb is not None and np.all(np.isfinite(bb.extents)):
            return bb
        renderer = renderer or ax.figure.canvas.get_renderer()
        boxes = [ax.get_window_extent(renderer)]
        for label in getattr(gl, "label_artists", []):
            if label.get_visible():
                b = label.get_window_extent(renderer)
                if np.all(np.isfinite(b.extents)):
                    boxes.append(b)
        for t in (ax.title, ax._left_title, ax._right_title):
            if t.get_text():
                b = t.get_window_extent(renderer)
                if np.all(np.isfinite(b.extents)):
                    boxes.append(b)
        return Bbox.union(boxes)

    ax.get_tightbbox = safe_tightbbox


def geo_map(projection=None, extent=None, nrows=1, ncols=1, width=10, ax=None,
            coastlines=True, gridlines=True, resolution="110m", **subplot_kw):
    """
    Create (or style) one or several cartopy map axes.

    projection : cartopy CRS for the *map*, default PlateCarree().
                 Your data's CRS goes in `transform=` when plotting.
    extent     : [lon_min, lon_max, lat_min, lat_max] in degrees, or None for global.
    nrows/ncols: grid of maps sharing the same projection and extent.
    width      : figure width in inches; the height is chosen to fit the maps
                 tightly. Don't call tight_layout() on the result.
    ax         : existing GeoAxes (or array of them) to style instead.
    Returns (fig, ax) like plt.subplots: a single axes, or an array for a grid.
    """
    import cartopy.crs as ccrs  # imported here so the module works without cartopy

    if ax is not None:
        for a in np.ravel(ax):
            _style_geo_axes(a, extent, coastlines, gridlines, resolution)
        return np.ravel(ax)[0].figure, ax

    projection = projection or ccrs.PlateCarree()
    fig, axes = plt.subplots(nrows, ncols, squeeze=False,
                             subplot_kw={"projection": projection, **subplot_kw})
    for a in axes.flat:
        _style_geo_axes(a, extent, coastlines, gridlines, resolution)

    # Fixed layout in inches: figure height follows the map aspect ratio,
    # so there is no dead space. Bottom margin leaves room for a colorbar.
    left, right, top, bottom = 0.7, 0.3, 0.55, 1.4
    vgap, hgap = 0.6, 0.35                       # room for titles / lat labels
    x0, x1, y0, y1 = axes[0, 0].get_extent()
    map_w = (width - left - right - (ncols - 1) * hgap) / ncols
    map_h = map_w * (y1 - y0) / (x1 - x0)
    height = top + bottom + nrows * map_h + (nrows - 1) * vgap
    fig.set_size_inches(width, height)
    fig.subplots_adjust(left=left / width, right=1 - right / width,
                        top=1 - top / height, bottom=bottom / height,
                        hspace=vgap / map_h, wspace=hgap / map_w)

    return fig, (axes[0, 0] if nrows == ncols == 1 else axes.squeeze())


if __name__ == "__main__":
    import numpy as np

    colors = setup_plotting()
    rng = np.random.default_rng(0)
    fig, (ax1, ax2) = plt.subplots(1, 2)

    # Distribution
    for mu, s in [(0, 1), (1.5, 0.7)]:
        ax1.hist(rng.normal(mu, s, 5000), bins=60, density=True,
                 histtype="stepfilled", alpha=0.6, lw=2,
                 label=rf"$\mu={mu},\ \sigma={s}$")
    ax1.set(xlabel=r"$x$", ylabel="Density")
    ax1.legend()

    # Map (NaNs show the grey background, handy for masked regions)
    x = np.linspace(-3, 3, 300)
    X, Y = np.meshgrid(x, x)
    Z = np.sin(X) * np.cos(Y)
    Z[X**2 + Y**2 > 8] = np.nan
    im = ax2.pcolormesh(X, Y, Z, cmap="RdBu_r", shading="auto")
    map_axes(ax2)
    add_colorbar(im, ax2, label=r"$\phi$")

    fig.tight_layout()
    fig.savefig("demo.png")
