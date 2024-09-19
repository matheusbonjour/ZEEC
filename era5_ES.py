import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import geopandas as gpd
import pandas as pd
import matplotlib
from pyproj import Proj, transform
from scipy.ndimage import gaussian_filter
import matplotlib.colors as mcolors


# PARTE DEDICADA A LEITURA DE DADOS # 

marem_path = "data/MAREM/"
clname = 'BR_MAREM_OP_pol_v2'
bs_gdf = gpd.read_file(marem_path + f'{clname}.shp')


ds = xr.open_dataset("data/era5_ES.nc")

#======================================================#

# CONFIGURACOES E ESTILIZACAO DO PLOT # 

plt.rcParams['font.size'] = 14
plt.rcParams['font.weight'] = 'bold'
plt.rcParams['axes.labelweight'] = 'bold'
plt.rcParams['axes.titleweight'] = 'bold'
plt.rcParams['axes.titlesize'] = 14
plt.rcParams['xtick.labelsize'] = 14
plt.rcParams['ytick.labelsize'] = 14
plt.rcParams['legend.fontsize'] = 14
plt.rcParams['legend.title_fontsize'] = 14

#================================================================#

# FUNCOES PARA VISUALIZACAO RAPIDA DE DADOS DO ERA5 # 

def plot_wind(ds, time_idx=0):

    fig, ax = plt.subplots(subplot_kw={'projection': ccrs.PlateCarree()})
    u = ds.u10.isel(time=time_idx)
    v = ds.v10.isel(time=time_idx)
    wind_speed = np.sqrt(u ** 2 + v ** 2)
    colors = ["#e7f2f4", "#ceeaee", "#b6e2e8", "#abdcff", "#a4d685", 
            "#abcf2a", "#ffba00", "#ffa200"]
    cmap2 = matplotlib.colors.ListedColormap(colors)
    cmap2.set_over('#ff8c00')
    cmap2.set_under('#fffafa')
    levels = np.arange(0, wind_speed.max()+1, 0.5)
    contour = ax.contourf(ds.longitude, ds.latitude, wind_speed, levels=levels, cmap=cmap2, transform=ccrs.PlateCarree(), extent='both')
    contorno = ax.contour(ds.longitude, ds.latitude, wind_speed,  colors='white', linewidths=0.5, levels=levels, zorder =50)
    cbar = plt.colorbar(contour, ax=ax, orientation='vertical')
    cbar.set_label('Magnitude do Vento (m/s)')
    ax.quiver(ds.longitude, ds.latitude, u, v, transform=ccrs.PlateCarree(), scale=100, regrid_shape=20)
    ax.coastlines()
    ax.add_feature(cfeature.BORDERS, linestyle=':')
    time_str = str(ds.time[time_idx].values)[:10]  # Formato de data: AAAA-MM-DD
    gl = ax.gridlines(draw_labels=True, linestyle='--', color='gray', alpha=0.7)
    gl.top_labels = False
    gl.right_labels = False
    gl.xlabel_style = {'size': 14, 'weight': 'bold'}
    gl.ylabel_style = {'size': 14, 'weight': 'bold'}
    plt.title(f'{time_str}')
    plt.show()

def plot_t2m(ds, time_idx=0):

    fig, ax = plt.subplots(subplot_kw={'projection': ccrs.PlateCarree()})
    t2m = ds.t2m.isel(time=time_idx)
    t2mC = t2m - 273.15
    levels = np.linspace(t2mC.min(), t2mC.max(), 50)
    contour = ax.contourf(ds.longitude, ds.latitude, t2mC, levels=levels, cmap='coolwarm', transform=ccrs.PlateCarree())
    cbar = plt.colorbar(contour, ax=ax, orientation='vertical')
    cbar.set_label('Temperatura a 2 metros [°C]', fontweight='bold')
    ax.coastlines()
    ax.add_feature(cfeature.BORDERS, linestyle=':')
    time_str = str(ds.time[time_idx].values)[:10]
    gl = ax.gridlines(draw_labels=True, linestyle='--', color='gray', alpha=0.7)
    gl.top_labels = False
    gl.right_labels = False
    gl.xlabel_style = {'size': 14, 'weight': 'bold'}
    gl.ylabel_style = {'size': 14, 'weight': 'bold'}
    plt.title(f'{time_str}', fontweight='bold')
    plt.show()

def plot_tp(ds, time_idx=0):

    fig, ax = plt.subplots(subplot_kw={'projection': ccrs.PlateCarree()})
    tp = ds.tp.isel(time=time_idx)
    tp_normalized = tp * 1000  
    levels = np.linspace(0, tp_normalized.max(), 40)
    contour = ax.contourf(ds.longitude, ds.latitude, tp_normalized, levels=levels, cmap='Blues', transform=ccrs.PlateCarree())
    cbar = plt.colorbar(contour, ax=ax, orientation='vertical')
    cbar.set_label('Precipitação total ($10^{-3}$ mm)', fontweight='bold')
    ax.coastlines()
    ax.add_feature(cfeature.BORDERS, linestyle=':')
    gl = ax.gridlines(draw_labels=True, linestyle='--', color='gray', alpha=0.7)
    gl.top_labels = False
    gl.right_labels = False
    gl.xlabel_style = {'size': 14, 'weight': 'bold'}
    gl.ylabel_style = {'size': 14, 'weight': 'bold'}
    time_str = str(ds.time[time_idx].values)[:10]
    plt.title(f'{time_str}', fontweight='bold')
    plt.show()



plot_wind(ds, time_idx=0)

plot_t2m(ds, time_idx=0)

plot_tp(ds, time_idx=0)