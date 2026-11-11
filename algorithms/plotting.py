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
from mpl_toolkits.axes_grid1 import make_axes_locatable

PALETTE = mpl.colormaps["Set2"].colors  # 8 soft, distinguishable colours


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


def _domain_path(extent, n=200):
    """Closed lon/lat path around `extent`, densified so it curves when projected."""
    import matplotlib.path as mpath

    lon0, lon1, lat0, lat1 = extent
    # A frame edge exactly at a pole collapses to one point, which older cartopy
    # cannot handle; stopping just short of it looks identical.
    lat0, lat1 = max(lat0, -89.0), min(lat1, 89.0)
    lon = np.linspace(lon0, lon1, n)
    lat = np.linspace(lat0, lat1, n)
    verts = np.vstack([
        np.column_stack([lon, np.full(n, lat0)]),         # south edge
        np.column_stack([np.full(n, lon1), lat]),         # east edge
        np.column_stack([lon[::-1], np.full(n, lat1)]),   # north edge
        np.column_stack([np.full(n, lon0), lat[::-1]]),   # west edge
    ])
    return mpath.Path(verts, closed=True)


def _nice_ticks(lo, hi, nbins=6):
    """Round gridline positions strictly inside (lo, hi)."""
    import matplotlib.ticker as mticker

    ticks = mticker.MaxNLocator(nbins=nbins, steps=[1, 2, 2.5, 3, 5, 10]).tick_values(lo, hi)
    return [float(t) for t in ticks if lo < t < hi]


def _fmt_lat(v):
    return f"{abs(v):g}°" + ("N" if v > 0 else "S" if v < 0 else "")


def _fmt_lon(v):
    v = (v + 180) % 360 - 180
    return f"{abs(v):g}°" + ("E" if 0 < v < 180 else "W" if v < 0 else "")


def _edge_labels(ax, extent, xlocs, ylocs, left=True, bottom=True, offset=6):
    """Latitude labels along the western edge, longitude labels along the southern edge."""
    import cartopy.crs as ccrs
    from matplotlib.transforms import offset_copy

    lon0, lon1, lat0, lat1 = extent
    data = ccrs.PlateCarree()._as_mpl_transform(ax)
    size = plt.rcParams["xtick.labelsize"]
    if left:
        tr = offset_copy(data, fig=ax.figure, x=-offset, units="points")
        for lat in ylocs:
            ax.text(lon0, lat, _fmt_lat(lat), transform=tr, ha="right", va="center",
                    size=size, clip_on=False)
    if bottom:
        tr = offset_copy(data, fig=ax.figure, y=-offset, units="points")
        for lon in xlocs:
            ax.text(lon, lat0, _fmt_lon(lon), transform=tr, ha="center", va="top",
                    size=size, clip_on=False)


def _inner_box(projection, extent, n=40):
    """
    Largest axis-aligned rectangle (in projected coordinates) that lies inside
    the lon/lat box `extent`, so a rectangular map has no empty corners.
    """
    import cartopy.crs as ccrs
    from shapely.geometry import LineString, Polygon

    verts = _domain_path(extent).vertices
    xy = projection.transform_points(ccrs.PlateCarree(), verts[:, 0], verts[:, 1])[:, :2]
    poly = Polygon(xy[np.all(np.isfinite(xy), axis=1)]).buffer(0)
    xmin, ymin, xmax, ymax = poly.bounds

    def x_interval(y):
        cut = poly.intersection(LineString([(xmin - 1, y), (xmax + 1, y)]))
        if cut.is_empty:
            return None
        parts = getattr(cut, "geoms", [cut])
        seg = max(parts, key=lambda g: g.length)
        return seg.bounds[0], seg.bounds[2]

    ys = np.linspace(ymin, ymax, n + 2)[1:-1]
    best, best_area = None, 0.0
    for i, y0 in enumerate(ys):
        for y1 in ys[i + 1:]:
            spans = [x_interval(y) for y in np.linspace(y0, y1, 12)]
            if any(sp is None for sp in spans):
                continue
            x0, x1 = max(sp[0] for sp in spans), min(sp[1] for sp in spans)
            area = (x1 - x0) * (y1 - y0)
            if x1 > x0 and area > best_area:
                best, best_area = (x0, x1, y0, y1), area
    x0, x1, y0, y1 = best
    dx, dy = 0.005 * (x1 - x0), 0.005 * (y1 - y0)   # tiny safety margin
    return [x0 + dx, x1 - dx, y0 + dy, y1 - dy]


def _style_geo_axes(ax, extent, coastlines, gridlines, resolution, curved=None, frame=None):
    import cartopy.crs as ccrs

    curved = False
    if extent is not None:
        if frame is None:
            frame = "box" if isinstance(ax.projection, ccrs.PlateCarree) else "curved"
        if frame == "fill":
            ax.set_extent(_inner_box(ax.projection, extent), crs=ax.projection)
        else:
            ax.set_extent(extent, crs=ccrs.PlateCarree())
        if frame == "curved":
            ax.set_boundary(_domain_path(extent), transform=ccrs.PlateCarree())
            curved = True
    else:
        ax.set_global()
    if coastlines:
        ax.coastlines(resolution=resolution, linewidth=1.0, color="k")
    if gridlines:
        ss = ax.get_subplotspec()
        label_left = ss is None or ss.is_first_col()
        label_bottom = ss is None or ss.is_last_row()
        line_style = dict(crs=ccrs.PlateCarree(), linestyle=plt.rcParams["grid.linestyle"],
                          linewidth=0.8, color="grey", alpha=0.7)
        if extent is not None and curved:
            # Curved frame: cartopy's label placement differs between versions
            # (labels on the wrong edges, or crashes), so label the edges ourselves.
            xlocs, ylocs = _nice_ticks(extent[0], extent[1]), _nice_ticks(extent[2], extent[3])
            ax.gridlines(draw_labels=False, xlocs=xlocs, ylocs=ylocs, **line_style)
            _edge_labels(ax, extent, xlocs, ylocs, label_left, label_bottom)
        else:
            gl = ax.gridlines(draw_labels=True, x_inline=False, y_inline=False,
                              rotate_labels=False, **line_style)
            gl.top_labels = gl.right_labels = False
            gl.bottom_labels = "x" if label_bottom else False   # longitudes at the bottom
            gl.left_labels = "y" if label_left else False       # latitudes on the left
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
            coastlines=True, gridlines=True, resolution="110m",
            vgap=0.6, hgap=0.35, frame=None, curved=None, **subplot_kw):
    """
    Create (or style) one or several cartopy map axes.

    projection : cartopy CRS for the *map*, default PlateCarree().
                 Your data's CRS goes in `transform=` when plotting.
    extent     : [lon_min, lon_max, lat_min, lat_max] in degrees, or None for global.
    nrows/ncols: grid of maps sharing the same projection and extent.
    width      : figure width in inches; the height is chosen to fit the maps
                 tightly. Don't call tight_layout() on the result.
    ax         : existing GeoAxes (or array of them) to style instead.
    frame      : how to frame a map with an `extent`:
                 "curved" - follow the lon/lat edges of the domain (fan shape for
                            conic/polar projections); default for non-PlateCarree.
                 "fill"   - largest rectangle that lies inside the domain: no
                            empty corners and no white space, but the outer
                            edges of the domain are cropped.
                 "box"    - rectangle around the whole domain (empty corners);
                            default for PlateCarree.
    curved     : old name; curved=True means frame="curved", False means "box".
    vgap/hgap  : space between maps in inches. Increase vgap (e.g. 1.4) when
                 every map gets its own colorbar underneath.
    Returns (fig, ax) like plt.subplots: a single axes, or an array for a grid.
    """
    import cartopy.crs as ccrs  # imported here so the module works without cartopy

    if frame is None and curved is not None:
        frame = "curved" if curved else "box"

    if ax is not None:
        for a in np.ravel(ax):
            _style_geo_axes(a, extent, coastlines, gridlines, resolution, frame=frame)
        return np.ravel(ax)[0].figure, ax

    projection = projection or ccrs.PlateCarree()
    fig, axes = plt.subplots(nrows, ncols, squeeze=False,
                             subplot_kw={"projection": projection, **subplot_kw})
    for a in axes.flat:
        _style_geo_axes(a, extent, coastlines, gridlines, resolution, frame=frame)

    # Fixed layout in inches: figure height follows the map aspect ratio,
    # so there is no dead space. Bottom margin leaves room for a colorbar.
    left, right, top, bottom = 0.7, 0.3, 0.55, 1.4
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