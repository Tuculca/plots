from netCDF4 import Dataset
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib.cm import get_cmap
import cartopy
import cartopy.crs as crs
from cartopy.feature import NaturalEarthFeature
from wrf import (to_np, getvar, smooth2d, get_cartopy, cartopy_xlim, cartopy_ylim, latlon_coords, ALL_TIMES)
import numpy
import os
import sys
import pandas
import cartopy.io.shapereader as shpreader
import cartopy.feature as cfeature



def plot(fileInput, out_folder):
	fileInput = Dataset(fileInput)
	times = getvar(fileInput, "Times")
	uv10 = getvar(fileInput, "uvmet10", units = "km h-1")
	wspd = getvar(fileInput, "uvmet10_wspd_wdir")
	#get_var = get_var[0]
	#lats     = getvar(fileInput, "XLAT_U")
	#lons     = getvar(fileInput, "XLONG_U")
	lats, lons = latlon_coords(uv10)
	cart_proj = get_cartopy(uv10)
	
	cmap = plt.get_cmap('magma_r')
	cmap = colors.LinearSegmentedColormap.from_list('trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=0, b=0.85), cmap(numpy.linspace(0, 0.85, 8))).with_extremes(over='black')
	#cmap=plt.get_cmap('magma_r', 8).with_extremes(over='black')
	bounds = numpy.linspace(0, 40, 9)
	norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
	
	fig = plt.figure(figsize=(12, 9), dpi=100)
	ax = fig.add_subplot(projection=cart_proj)

	ax.set_xlim(cartopy_xlim(uv10))
	ax.set_ylim(cartopy_ylim(uv10))

	clr = ax.contourf(to_np(lons), to_np(lats), to_np(wspd[0]), cmap=cmap, norm=norm, zorder=1, transform=crs.PlateCarree())
	#ax.barbs(to_np(lons)[::20,::20], to_np(lats)[::20,::20], to_np(uv10[0])[::20,::20], to_np(uv10[1])[::20,::20], zorder=1, transform=crs.PlateCarree())
	ax.quiver(to_np(lons)[::20,::20], to_np(lats)[::20,::20], to_np(-numpy.sin(-numpy.radians(wspd[1])))[::20,::20], to_np(-numpy.cos(numpy.radians(wspd[1])))[::20,::20], zorder=1, transform=crs.PlateCarree())
	
	
	ax.add_feature(cfeature.BORDERS)
	
	fname = 'C:/Users/Aless/Downloads/gadm41_ITA_shp/gadm41_ITA_1.shp'
	adm1_shapes = list(shpreader.Reader(fname).geometries())
	ax.add_geometries(adm1_shapes, crs.PlateCarree(), edgecolor='black', facecolor='none', zorder=2) #faceocolor era "gray"

	plt.title("Vento a 10 metri " + str(pandas.to_datetime(to_np(times))))
	
	cbar = fig.colorbar(mpl.cm.ScalarMappable(cmap=cmap, norm=norm),
	extend='max',
	extendfrac='auto',
	ticks=bounds,
	spacing='uniform')
	cbar.set_label("km/h", rotation=0)
	#ax.gridlines(color="black", linestyle="dotted")
	plt.savefig(out_folder + str(pandas.to_datetime(to_np(times))).replace(':00:00','') + '.jpg', bbox_inches='tight')

files = sys.argv[1]
output_folder = sys.argv[2] if len(sys.argv)> 2 else 'D:\Model/'
#files = sorted(os.listdir(input_folder))
#for file in files:
#	file_path = input_folder + file
#	ncfiles.append(Dataset(file_path))

	
plot(files, output_folder)
