#!/usr/local/bin/python

from netCDF4 import Dataset,date2num,num2date
from mpl_toolkits.basemap import Basemap, shiftgrid, cm
import numpy as np
import matplotlib.pyplot as plt
import matplotlib


infile = '/usr1/ymao/other/GRACE/plotting/trend_map/input/latlon.xy'
f=np.genfromtxt(infile,dtype=float)
lats=np.unique(f[:,0])
lats=-lats
lats.sort()
lats=-lats
lons=np.unique(f[:,1])
lons, lats = np.meshgrid(lons, lats)
infile2 = '/usr1/ymao/other/GRACE/USGS_well_data/output/trend_mm_yr.asc'
trends=np.genfromtxt(infile2,skip_header=6,dtype=float)
trends = np.ma.masked_where(trends==-9999, trends)
trends = np.ma.masked_where(trends==0, trends)
fig = plt.figure()
ax = fig.add_axes([0.1,0.1,0.8,0.8])
m = Basemap(llcrnrlon=-120.,llcrnrlat=20.,urcrnrlon=-60.,urcrnrlat=50.,            rsphere=(6378137.00,6356752.3142),            resolution='l',area_thresh=1000.,projection='lcc',            lat_1=50.,lon_0=-107.,ax=ax)
x,y=m(lons,lats)
p = plt.pcolormesh(x,y,trends,cmap=cm.GMT_no_green_r,vmin=-20,vmax=20)
m.drawcoastlines()
m.drawcountries()
m.drawstates()
parallels = np.arange(0.,80,10.)
m.drawparallels(parallels,labels=[1,0,0,1])
meridians = np.arange(10.,360.,10.)
m.drawmeridians(meridians,labels=[1,0,0,1])
#ticks = np.array([-1.2, -1.0, -0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2], np.float)
cbar = plt.colorbar(p, orientation='vertical', fraction=0.05, pad=0.02)
cbar.set_label('Trend (mm/year)', fontsize=16)
plt.title('Aug GW storage trend, CRN, \n2002-2013, sites deeper than 100 ft')  ############# change this ##############
plt.savefig('/usr1/ymao/other/GRACE/USGS_well_data/plots/trend_basin_Aug_daily_climNet_depthGT100.png', format='png') ############ change this ###########


infile = '/usr1/ymao/other/GRACE/plotting/trend_map/input/latlon.xy'
f=np.genfromtxt(infile,dtype=float)
lats=np.unique(f[:,0])
lats=-lats
lats.sort()
lats=-lats
lons=np.unique(f[:,1])
infile2 = '/usr1/ymao/other/GRACE/USGS_well_data/output/sd_resid_mm.asc'
data=np.genfromtxt(infile2,skip_header=6,dtype=float)
fig = plt.figure()
ax = fig.add_axes([0.1,0.1,0.8,0.8])
datam = np.ma.masked_where(data==-9999, data)
datam = np.ma.masked_where(datam==0, datam)
print infile+' min = '+str(datam.min())+' max = '+str(datam.max())
m = Basemap(llcrnrlon=-120.,llcrnrlat=20.,urcrnrlon=-60.,urcrnrlat=50.,            rsphere=(6378137.00,6356752.3142),            resolution='l',area_thresh=1000.,projection='lcc',            lat_1=50.,lon_0=-107.,ax=ax)
lons, lats = np.meshgrid(lons, lats)
x,y=m(lons,lats)
p = plt.pcolormesh(x,y,datam,cmap=cm.GMT_no_green_r,vmin=0,vmax=100)
m.drawcoastlines()
m.drawcountries()
m.drawstates()
parallels = np.arange(0.,80,10.)
m.drawparallels(parallels,labels=[1,0,0,1])
meridians = np.arange(10.,360.,10.)
m.drawmeridians(meridians,labels=[1,0,0,1])
cbar = plt.colorbar(p, orientation='vertical', fraction=0.05, pad=0.02)
cbar.set_label('Standard deviation (mm)', fontsize=16)
plt.title('Standard deviation in residual Aug GW storage, \nCRN, 2002-2013, sites deeper than 100 ft')
plt.savefig('/usr1/ymao/other/GRACE/USGS_well_data/plots/resid_basin_Aug_daily_CRN_depthGT100.png', format='png') 

