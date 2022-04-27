import matplotlib
from cartopy.mpl.gridliner import LongitudeFormatter, LatitudeFormatter
import matplotlib.pyplot as plt
import cartopy.crs as ccrs 
from matplotlib.ticker import ScalarFormatter,AutoMinorLocator

def global_map(ax):
    ax.set_extent([-180,180,-90,90], crs=ccrs.PlateCarree())
    ax.set_xticks([-180,-120,-60,0,60,120,180], crs=ccrs.PlateCarree())
    ax.set_yticks([-90,-60,-30,0,30,60,90], crs=ccrs.PlateCarree())
    lon_formatter = LongitudeFormatter(zero_direction_label=True)
    lat_formatter = LatitudeFormatter()
    ax.xaxis.set_major_formatter(lon_formatter)
    ax.yaxis.set_major_formatter(lat_formatter)
#     ax.yaxis.set_major_formatter(ScalarFormatter())
    ax.yaxis.set_minor_locator(  AutoMinorLocator(3))
    ax.xaxis.set_minor_locator(  AutoMinorLocator(2))
    ax.coastlines()
    ax.set_ylabel(' ')
    ax.set_xlabel(' ')




def create_facet_plot(nplots: int, figsize: tuple, 
                        subplot_kw : dict ={'projection':ccrs.PlateCarree()},
                        last_axis_plain: bool = False):
    layout = _get_layouts(nplots)
    fig ,axes = plt.subplot_mosaic(layout,figsize=figsize,subplot_kw=subplot_kw)
    last = list(axes.keys())[-1]
    
    if last_axis_plain == True:
        pos = axes[last].get_position()
        axes[last] = fig.add_subplot(axes[last].get_gridspec().nrows, 
                                        axes[last].get_gridspec().ncols, nplots)
    return fig, axes
    



def _get_layouts(nplots: int):
    if nplots==1:
        layout="""
            A
            """
    elif nplots==2:
        layout="""
            AB
            """
    elif nplots==3:
        layout="""
            AABB
            .CC.
        """
    elif nplots==4:
        layout="""
            AB
            CD
        """
    elif nplots==5:
        layout="""
            AABB
            CCDD
            .EE.
        """
    elif nplots==6:
        layout="""
            AB
            CD
            EF
        """
    elif nplots==7:
        layout="""
            ABC
            DE.
            FG.
        """
    elif nplots==8:
        layout="""
            ABC
            DEF
            GH.
        """
    elif nplots==9:
        layout="""
            ABC
            DEF
            GHI
        """
    elif nplots==10:
        layout =="""
            ABC
            DEF
            GHI
            L
        """
    elif nplots==11:
        layout="""
            ABC
            DEF
            GHI
            LLM
        
        """
    elif nplots==12:
        layout="""
            ABCD
            EFGH
            ILMN
        """
    else:
        raise(NotImplementedError(f'Layout undefined for nplots = {nplots}'))

    return layout