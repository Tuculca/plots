from netCDF4 import Dataset
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib.cm import get_cmap
import cartopy
import cartopy.crs as crs
from cartopy.feature import NaturalEarthFeature
from wrf import (to_np, getvar, smooth2d, get_cartopy, cartopy_xlim, cartopy_ylim, latlon_coords, interplevel)
import numpy
import os
import sys
import pandas
import cartopy.io.shapereader as shpreader
import cartopy.feature as cfeature


def plot(fileInput, out_folder):
	fileInput = Dataset(fileInput)
	times = getvar(fileInput, "Times")
	slp = getvar(fileInput, "slp")
	slp = smooth2d(slp, 8)
	p = getvar(fileInput, "pressure")
	z = getvar(fileInput, 'z')
	lats, lons = latlon_coords(z)
	cart_proj = get_cartopy(z)
	
	
	#colours=['#d2fffe', '#88fefd', '#00c6ff', '#1996ff', '#3c41ff', '#3cbc3d', '#a5d71f', '#ffe600', '#ffc300', '#ff7d00', '#ff0000', '#c80000', '#d464c3', '#b5199d', '#840094', '#dcdcdc', '#b4b4b4', '#8c8c8c', '#5a5a5a']
	#cmap = (mpl.colors.ListedColormap(colours).with_extremes(over='black', under='white'))
	
	ht_500 = interplevel(z, p, 500)
	ht_500 = to_np(ht_500)/10
	Vmin = to_np(slp).min()
	Vmax = to_np(slp).max()
	Vmin=numpy.floor(Vmin)
	Vmax=numpy.ceil(Vmax)
	#cmap=plt.get_cmap('viridis', int(Vmax-Vmin))
	boundsSLP = numpy.linspace(Vmin, Vmax, int(Vmax-Vmin)+1)
	cmap=plt.get_cmap('jet', 20)
	bounds = numpy.linspace(480,600,21)
	norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
	
	fig = plt.figure(figsize=(16, 12), dpi=120)
	ax = fig.add_subplot(projection=cart_proj)
	get_var_np = to_np(ht_500)
	#get_var_np[get_var_np<=3] = numpy.nan
	ax.set_xlim(cartopy_xlim(z))
	ax.set_ylim(cartopy_ylim(z))
	#clr = ax.pcolormesh(to_np(lons), to_np(lats), get_var_np, transform=crs.PlateCarree(), cmap='plasma', shading='gouraud', zorder=3, vmin=Vmin, vmax=Vmax)
	ax.contourf(to_np(lons), to_np(lats), get_var_np, levels=bounds, cmap=cmap, norm=norm, zorder=1, transform=crs.PlateCarree())
	cs = ax.contour(to_np(lons), to_np(lats), to_np(slp), levels=boundsSLP, colors='white', zorder=4, linewidths=2, transform=crs.PlateCarree())
	
	#land_50m = cfeature.NaturalEarthFeature('physical', 'land', '50m', edgecolor='face', facecolor=cfeature.COLORS['land'])
	#ax.add_feature(land_50m)
	
	#ocean_50m = cfeature.NaturalEarthFeature('physical', 'ocean', '50m', edgecolor='face')
	#ax.add_feature(ocean_50m)
	
	ax.add_feature(cfeature.COASTLINE)
	#ax.add_feature(cfeature.BORDERS)
	
	fname = 'C:/Users/Aless/Downloads/gadm41_ITA_shp/gadm41_ITA_1.shp'
	adm1_shapes = list(shpreader.Reader(fname).geometries())
	ax.add_geometries(adm1_shapes, crs.PlateCarree(), edgecolor='black', facecolor='none', zorder=2) #faceocolor era "gray"

	plt.title("Altezza superficie 500hPa e Pressione al livello del mare" + str(pandas.to_datetime(to_np(times))))
	
	cbar = fig.colorbar(mpl.cm.ScalarMappable(cmap=cmap, norm=norm),
	extend='neither',
	extendfrac='auto',
	ticks=bounds,
	spacing='uniform')
	#cbar=fig.colorbar(clr)
	cbar.set_label("dam", rotation=0)
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
