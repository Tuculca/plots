from netCDF4 import Dataset
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib.cm import get_cmap
import cartopy
import cartopy.crs as crs
from cartopy.feature import NaturalEarthFeature
from wrf import (to_np, getvar, smooth2d, get_cartopy, cartopy_xlim, cartopy_ylim, latlon_coords, ALL_TIMES)
import numpy as np
import os
import sys
import pandas
import cartopy.io.shapereader as shpreader
import cartopy.feature as cfeature


def plot(fileInput, out_folder):
	fileInput = Dataset(fileInput)
	var = "RAINNC"
	times = getvar(fileInput, "Times")
	get_var = getvar(fileInput, var)
	#get_var = get_var[0]
	#lats     = getvar(fileInput, "XLAT_U")
	#lons     = getvar(fileInput, "XLONG_U")
	Vmin = to_np(get_var).min()
	Vmax = to_np(get_var).max()
	lats, lons = latlon_coords(get_var)
	cart_proj = get_cartopy(get_var)
	
	
	colours=['#d2fffe', '#88fefd', '#00c6ff', '#1996ff', '#3c41ff', '#3cbc3d', '#a5d71f', '#ffe600', '#ffc300', '#ff7d00', '#ff0000', '#c80000', '#d464c3', '#b5199d', '#840094', '#dcdcdc', '#b4b4b4', '#8c8c8c', '#5a5a5a']
	cmap = (mpl.colors.ListedColormap(colours).with_extremes(over='black', under='white'))
	bounds = [0.2,1,3,5,7,10,15,20,25,30,40,50,60,70,80,100,125,150,175,200]
	norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
	
	fig = plt.figure(figsize=(16, 12), dpi=120)
	ax = fig.add_subplot(projection=cart_proj)
	get_var_np = to_np(get_var)
	#get_var_np[get_var_np<=3] = np.nan
	ax.set_xlim(cartopy_xlim(get_var))
	ax.set_ylim(cartopy_ylim(get_var))
	#clr = ax.pcolormesh(to_np(lons), to_np(lats), get_var_np, transform=crs.PlateCarree(), cmap='plasma', shading='gouraud', zorder=3, vmin=Vmin, vmax=Vmax)
	clr = ax.contourf(to_np(lons), to_np(lats), get_var_np, levels=bounds, cmap=cmap, norm=norm, zorder=3, transform=crs.PlateCarree(), alpha=0.85, extent=(lons.min(), lons.max(), lats.min(), lats.max()))
	
	#land_50m = cfeature.NaturalEarthFeature('physical', 'land', '50m', edgecolor='face', facecolor=cfeature.COLORS['land'])
	#ax.add_feature(land_50m)
	
	ocean_50m = cfeature.NaturalEarthFeature('physical', 'ocean', '50m', edgecolor='face')
	ax.add_feature(ocean_50m)
	
	ax.add_feature(cfeature.BORDERS)
	
	fname = 'C:/Users/Aless/Downloads/gadm41_ITA_shp/gadm41_ITA_1.shp'
	adm1_shapes = list(shpreader.Reader(fname).geometries())
	ax.add_geometries(adm1_shapes, crs.PlateCarree(), edgecolor='black', facecolor='gray', alpha=0.15)

	plt.title(get_var.description[0].upper() + get_var.description[1:].lower() + " " + \
	str(pandas.to_datetime(to_np(times))) + ", " + get_var.units)
	
	cbar = fig.colorbar(mpl.cm.ScalarMappable(cmap=cmap, norm=norm),
	extend='both',
	extendfrac='auto',
	ticks=bounds,
	spacing='uniform')
	#cbar=fig.colorbar(clr)
	cbar.set_label(get_var.units+'/h', rotation=0)
	ax.gridlines(color="black", linestyle="dotted")
	plt.savefig(out_folder + str(pandas.to_datetime(to_np(times))).replace(':00:00','') + '.jpg', bbox_inches='tight')

files = sys.argv[1]
output_folder = sys.argv[2] if len(sys.argv)> 2 else 'D:\Model/'
#files = sorted(os.listdir(input_folder))
#for file in files:
#	file_path = input_folder + file
#	ncfiles.append(Dataset(file_path))

	
plot(files, output_folder)
