import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import numpy as np
import xarray as xr
import geopandas as gpd
import cartopy.feature as cfeature
from pyproj import Proj, transform
from scipy.ndimage import gaussian_filter
import matplotlib.colors as mcolors


# PALETA DE COR PERSONALIZADA # 

levels = [0, 30, 60, 90, 120, 150, 180, 210, 240, 270, 300]
colors = ["white", "#add8e6", "#0000ff", "#0000ff", "#0000cc", "#000099", "#660066", "#800080", "#800080", "#800080"] 
cmap, norm = mcolors.from_levels_and_colors(levels, colors, extend='both')


# ===================================================== #
# PARTE DEDICADA A LEITURA DE DADOS #

marem_path = "data/MAREM/"
clname = 'BR_MAREM_OP_pol_v2'
bs_gdf = gpd.read_file(marem_path + f'{clname}.shp')
ds = xr.open_dataset("data/river_2022.nc")

#======================================================#


dec_dis24 = ds['dis24'].sel(time='2022-12-01').mean(dim='number')

# VISUALIZACAO RAPIDA DA DESCARGA DOS RIOS #

fig, ax = plt.subplots(figsize=(10, 8), subplot_kw={'projection': ccrs.PlateCarree()})
ax.set_extent([ds.longitude.min(), ds.longitude.max(), ds.latitude.min(), ds.latitude.max()], crs=ccrs.PlateCarree())
ax.coastlines(color='lightgrey') 
discharge_plot = ax.contourf(ds['longitude'], ds['latitude'], dec_dis24, levels=levels, cmap=cmap, norm=norm, extend='both')
cbar = plt.colorbar(discharge_plot, ax=ax, orientation='vertical', extend='both')
cbar.set_label('Transporte Volumétrico')
ax.gridlines(draw_labels=True)
plt.title('Dezembro - 2022')
plt.show()