import matplotlib
from cartopy.mpl.gridliner import LongitudeFormatter, LatitudeFormatter
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from matplotlib.ticker import ScalarFormatter, AutoMinorLocator


def global_map(ax):
    ax.set_extent([-180, 180, -90, 90], crs=ccrs.PlateCarree())
    ax.set_xticks([-180, -120, -60, 0, 60, 120, 180], crs=ccrs.PlateCarree())
    ax.set_yticks([-90, -60, -30, 0, 30, 60, 90], crs=ccrs.PlateCarree())
    lon_formatter = LongitudeFormatter(zero_direction_label=True)
    lat_formatter = LatitudeFormatter()
    ax.xaxis.set_major_formatter(lon_formatter)
    ax.yaxis.set_major_formatter(lat_formatter)
    #     ax.yaxis.set_major_formatter(ScalarFormatter())
    ax.yaxis.set_minor_locator(AutoMinorLocator(3))
    ax.xaxis.set_minor_locator(AutoMinorLocator(2))
    ax.coastlines()
    ax.set_ylabel(" ")
    ax.set_xlabel(" ")


def create_facet_plot(
    nplots: int,
    figsize: tuple = None,
    subplot_kw: dict = {"projection": ccrs.PlateCarree()},
    last_axis_plain: bool = False,
    create_cax=True
):
    layout, default_figsize, cax_loc = _get_layouts(nplots)
    if not figsize:
        figsize = default_figsize
    fig, axes = plt.subplot_mosaic(layout, figsize=figsize, subplot_kw=subplot_kw)
    last = list(axes.keys())[-1]

    if last_axis_plain == True:
        pos = axes[last].get_position()
        axes[last] = fig.add_subplot(
            axes[last].get_gridspec().nrows, axes[last].get_gridspec().ncols, nplots
        )
    if create_cax:
        cax = fig.add_axes(cax_loc)
    else: 
        cax=None
    return fig, axes, cax


def _get_layouts(nplots: int):
    cax_loc = None
    figsize = None
    if nplots == 1:
        layout = """
            A
            """
        figsize = (6, 5)
    elif nplots == 2:
        layout = """
            AB
            """
        figsize = (7, 4)
    elif nplots == 3:
        layout = """
            AABB
            .CC.
        """
        figsize = (8, 6)
    elif nplots == 4:
        layout = """
            AB
            CD
        """
        figsize = (8, 6)
    elif nplots == 5:
        layout = """
            AABB
            CCDD
            .EE.
        """
        figsize = (14, 12)
        cax_loc = [0.94,0.2,0.02,0.62]
    elif nplots == 6:
        layout = """
            AB
            CD
            EF
        """
        figsize = (14, 12)
        cax_loc = [0.94,0.2,0.02,0.62]
    elif nplots == 7:
        layout = """
            ABC
            DE.
            FG.
        """
        figsize = (16, 10)
        cax_loc = [0.94,0.2,0.02,0.62]
    elif nplots == 8:
        layout = """
            ABC
            DEF
            GH.
        """
        figsize = (16, 10)
        cax_loc = [0.94,0.2,0.02,0.62]
    elif nplots == 9:
        layout = """
            ABC
            DEF
            GHI
        """
        figsize = (16, 10)
        cax_loc = [0.94,0.2,0.02,0.62]
    elif nplots == 10:
        layout == """
            ABC
            DEF
            GHI
            L
        """
    elif nplots == 11:
        layout = """
            ABC
            DEF
            GHI
            LLM
        
        """
    elif nplots == 12:
        layout = """
            ABCD
            EFGH
            ILMN
        """
    else:
        raise (NotImplementedError(f"Layout undefined for nplots = {nplots}"))
    if not cax_loc:
        cax_loc=[0.94,0.2,0.02,0.62]

    return layout, figsize, cax_loc
