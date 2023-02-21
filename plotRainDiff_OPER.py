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
import sys
import pandas
import cartopy.io.shapereader as shpreader
import cartopy.feature as cfeature


def plot(file1, file2, out_folder):
	file1 = Dataset(file1)
	file2 = Dataset(file2)
	var = "RAINNC"
	hgt = getvar(file1, "HGT")
	times = getvar(file2, "Times")
	get_var1 = getvar(file1, var)
	get_var2 = getvar(file2, var)
	get_var = get_var2 - get_var1
	get_var_np = to_np(get_var)
	#get_var_np[get_var_np<=2] = np.nan
	#get_var = get_var[0]
	#lats     = getvar(fileInput, "XLAT_U")
	#lons     = getvar(fileInput, "XLONG_U")
	Vmin = to_np(get_var).min()
	Vmax = to_np(get_var).max()
	lats, lons = latlon_coords(hgt)
	cart_proj = get_cartopy(hgt)
	
	
	colours=['#d2fffe', '#88fefd', '#00c6ff', '#1996ff', '#3c41ff', '#3cbc3d', '#a5d71f', '#ffe600', '#ffc300', '#ff7d00', '#ff0000', '#c80000', '#d464c3', '#b5199d', '#840094', '#dcdcdc', '#b4b4b4', '#8c8c8c', '#5a5a5a']
	cmap = (mpl.colors.ListedColormap(colours).with_extremes(over='black', under='white'))
	bounds = [0.2,1,3,5,7,10,15,20,25,30,40,50,60,70,80,100,125,150,175,200]
	norm = mpl.colors.BoundaryNorm(bounds, cmap.N, clip=True)
	
	fig = plt.figure(figsize=(16, 12), dpi=120)
	ax = fig.add_subplot(projection=cart_proj)

	#clr = ax.pcolormesh(to_np(lons), to_np(lats), get_var_np, transform=crs.PlateCarree(), cmap='plasma', shading='gouraud', zorder=3, vmin=Vmin, vmax=Vmax)
	clr = ax.contourf(to_np(lons), to_np(lats), get_var_np, levels=bounds, colors=colours, zorder=3, transform=crs.PlateCarree(), alpha=0.85)
	
	
	land_50m = cfeature.NaturalEarthFeature('physical', 'land', '50m', edgecolor='face', facecolor=cfeature.COLORS['land'])
	ax.add_feature(land_50m)
	
	ocean_50m = cfeature.NaturalEarthFeature('physical', 'ocean', '50m', edgecolor='face')
	ax.add_feature(ocean_50m)
	
	ax.add_feature(cfeature.BORDERS)
	
	fname = 'C:/Users/Aless/Downloads/gadm41_ITA_shp/gadm41_ITA_1.shp'
	adm1_shapes = list(shpreader.Reader(fname).geometries())
	ax.add_geometries(adm1_shapes, crs.PlateCarree(), edgecolor='black', facecolor='gray', alpha=0.15)

	#plt.title(get_var2.description[0].upper() + get_var2.description[1:].lower() + " " + \
	#str(pandas.to_datetime(to_np(times))) + ", " + get_var2.units)
	plt.title("Precipitazione in 3 ore " + \
	str(pandas.to_datetime(to_np(times))))
	
	cbar = fig.colorbar(mpl.cm.ScalarMappable(cmap=cmap, norm=norm),
	extend='both',
	extendfrac='auto',
	ticks=bounds,
	spacing='uniform')
	#cbar=fig.colorbar(clr)
	cbar.set_label(get_var2.units+'/3h', rotation=0)
	
	plt.savefig(out_folder + str(pandas.to_datetime(to_np(times))).replace(':00:00','') + '.jpg', bbox_inches='tight')


output_folder = sys.argv[3] if len(sys.argv)> 3 else 'D:\Model/'

fileStart = sys.argv[1]
fileEnd = sys.argv[2]
	
plot(fileStart, fileEnd, output_folder)
