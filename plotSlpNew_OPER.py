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
	var = "slp"
	times = getvar(fileInput, "Times")
	get_var = getvar(fileInput, var)
	#get_var = get_var[0]
	#lats     = getvar(fileInput, "XLAT_U")
	#lons     = getvar(fileInput, "XLONG_U")
	Vmin = to_np(get_var).min()
	Vmax = to_np(get_var).max()
	lats, lons = latlon_coords(get_var)
	cart_proj = get_cartopy(get_var)
	
	
	#colours=['#d2fffe', '#88fefd', '#00c6ff', '#1996ff', '#3c41ff', '#3cbc3d', '#a5d71f', '#ffe600', '#ffc300', '#ff7d00', '#ff0000', '#c80000', '#d464c3', '#b5199d', '#840094', '#dcdcdc', '#b4b4b4', '#8c8c8c', '#5a5a5a']
	#cmap = (mpl.colors.ListedColormap(colours).with_extremes(over='black', under='white'))
	
	Vmin=numpy.floor(Vmin)
	Vmax=numpy.ceil(Vmax)
	cmap=plt.get_cmap('viridis', int(Vmax-Vmin))
	bounds = numpy.linspace(Vmin, Vmax, int(Vmax-Vmin+1))
	norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
	
	fig = plt.figure(figsize=(16, 12), dpi=120)
	ax = fig.add_subplot(projection=cart_proj)
	get_var_np = to_np(get_var)
	#get_var_np[get_var_np<=3] = numpy.nan
	ax.set_xlim(cartopy_xlim(get_var))
	ax.set_ylim(cartopy_ylim(get_var))
	#clr = ax.pcolormesh(to_np(lons), to_np(lats), get_var_np, transform=crs.PlateCarree(), cmap='plasma', shading='gouraud', zorder=3, vmin=Vmin, vmax=Vmax)
	ax.contourf(to_np(lons), to_np(lats), get_var_np, levels=bounds, cmap=cmap, norm=norm, zorder=1, transform=crs.PlateCarree())
	cs = ax.contour(to_np(lons), to_np(lats), get_var_np, levels=bounds, colors='none', zorder=4, transform=crs.PlateCarree())
	
	#land_50m = cfeature.NaturalEarthFeature('physical', 'land', '50m', edgecolor='face', facecolor=cfeature.COLORS['land'])
	#ax.add_feature(land_50m)
	
	#ocean_50m = cfeature.NaturalEarthFeature('physical', 'ocean', '50m', edgecolor='face')
	#ax.add_feature(ocean_50m)
	
	ax.add_feature(cfeature.BORDERS)
	ax.add_feature(cfeature.COASTLINE, zorder=3, alpha=0.35)
	
	fname = 'C:/Users/Aless/Downloads/gadm41_ITA_shp/gadm41_ITA_1.shp'
	adm1_shapes = list(shpreader.Reader(fname).geometries())
	ax.add_geometries(adm1_shapes, crs.PlateCarree(), edgecolor='black', facecolor='none', alpha=0.35, zorder=2) #faceocolor era "gray"

	plt.title("Pressione a livello del mare " + str(pandas.to_datetime(to_np(times))))
	
	cbar = fig.colorbar(mpl.cm.ScalarMappable(cmap=cmap, norm=norm),
	extend='neither',
	extendfrac='auto',
	ticks=bounds,
	spacing='uniform')
	#cbar=fig.colorbar(clr)
	cbar.set_label("hPa", rotation=0)
	ax.clabel(cs, inline=True, fontsize=10, colors='black')
	#ax.gridlines(color="black", linestyle="dotted")
	plt.savefig(out_folder + str(pandas.to_datetime(to_np(times))).replace(':00:00','') + '.jpg', bbox_inches='tight')

files = sys.argv[1]
output_folder = sys.argv[2] if len(sys.argv)> 2 else 'D:\Model/'
#files = sorted(os.listdir(input_folder))
#for file in files:
#	file_path = input_folder + file
#	ncfiles.append(Dataset(file_path))

	
plot(files, output_folder)
