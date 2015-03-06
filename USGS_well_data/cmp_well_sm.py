#!/usr/local/bin/python

import numpy as np
import scipy.stats
import os
import datetime as dt
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import statsmodels.api as sm

def cmap_drought_listed(vmin=-1, vmax=1):
    '''listed drought colormap'''

    vmin = float(vmin)
    vmax = float(vmax)

    bounds = np.array([-1.0, -0.8, -0.6, -0.4, -0.2, 0.2, 0.4, 0.6, 0.8, 1.0], np.float)
    red = np.array([230, 230, 254, 254, 255, 170, 76, 56, 20],
                   np.float)
    green = np.array([0, 152, 211, 254, 255, 245, 230, 168, 90],
                     np.float)
    blue = np.array([0, 0, 127, 0, 255, 150, 0, 0, 0], np.float)

    ticks = bounds
    colors = np.transpose(np.array([red/255., green/255., blue/255.]))
#    bounds = (vmax - vmin) / (bounds.max()-bounds.min()) * bounds + vmin
    norm = mpl.colors.BoundaryNorm(bounds, colors.shape[0])
    cmap = mpl.colors.ListedColormap(colors, 'drought')
    cmap.set_bad('lightgray')

    return cmap, norm, ticks, bounds

#########################################################################
############################# user defined ##############################
#########################################################################
well_data_path = '/usr1/ymao/other/GRACE/USGS_well_data/output/ts_Aug_CRN_and_FLD' # [siteID] [lat] [lon] [well depth (ft)] [12 values; Aug GW storage (mm)]; This is one file for all sites
vic_Aug_output_dir = '/usr1/ymao/other/GRACE/USGS_well_data/soil_moist_data/vic_2002_2013_Aug_all_old_data'  # directory for vic output; only contain data for Aug 1 each year during 2002-2013
plot_output_dir = '/usr1/ymao/other/GRACE/USGS_well_data/plots'
missing_value = -9999
nyear = 12
frac = 0.65  # f factor in LOWESS smothing

#########################################################################
############################# load data #################################
#########################################################################
well_data = np.loadtxt(well_data_path)
nsite = np.shape(well_data)[0]

#########################################################################
############################# process well data #########################
#########################################################################
# determine the 0.5 deg grid cell that this site falls into
well_data_in_grid = np.empty([nsite, 15]) # [lat_grid; lon_grid; well_depth(ft); 12 GW storage values (mm)]
for i in range(nsite):
	lat = well_data[i,1]
	lon = well_data[i,2]
	if lat>=round(lat):
		lat_grid = round(lat) + 0.25
	else:
		lat_grid = round(lat) - 0.25
	if lon>=round(lon):
		lon_grid = round(lon) + 0.25
	else:
		lon_grid = round(lon) - 0.25
	well_data_in_grid[i,0] = lat_grid
	well_data_in_grid[i,1] = lon_grid
	well_data_in_grid[i,2:15] = well_data[i,3:16]

# find sites in the same grid cell and take medium of well depth and GW storage values
well_data_uniq_list = []
flag = np.zeros(nsite)
for i in range(nsite):
	if flag[i]==0:
		temp = []
		temp.append(well_data_in_grid[i,:])
		ind = []
		ind.append(i)
		for j in range(nsite):
			if i!=j and well_data_in_grid[i,0]==well_data_in_grid[j,0] and well_data_in_grid[i,1]==well_data_in_grid[j,1]:  # if another site is in the same grid cell
				temp.append(well_data_in_grid[j,:])
				ind.append(j)
		if len(temp)==1:  # if this site is the only site in the grid cell
			well_data_uniq_list.append(well_data_in_grid[i,:])
			flag[i] = 1
		else:  # if this site is in the same grid cell as one or more other sites
			temp = np.asarray(temp)
			well_depth = np.median(temp[:,2])
			ts = np.empty(nyear)
			for j in range(nyear):
				ts[j] = np.median(temp[:,j+3])
			new = np.empty(15)
			new[0] = temp[0,0]
			new[1] = temp[0,1]
			new[2] = well_depth
			new[3:15] = ts
			well_data_uniq_list.append(new)
			for j in range(len(ind)):
				flag[ind[j]] = 1

# delete grid cells that do not have soil moisture data
well_data_cell_list = []
for i in range(len(well_data_uniq_list)):
	filename = '%s/fluxes_%.4f_%.4f' %(vic_Aug_output_dir, well_data_uniq_list[i][0], well_data_uniq_list[i][1])
	if os.path.isfile(filename)==True:  # if this grid cell exist
		well_data_cell_list.append(well_data_uniq_list[i])

well_data_cell = np.asarray(well_data_cell_list)
ncell = np.shape(well_data_cell)[0]

#########################################################################
####################### load VIC soil moisture data #####################
#########################################################################
vic_sm = np.empty([ncell, nyear])
for i in range(ncell):
	temp = np.loadtxt('%s/fluxes_%.4f_%.4f' %(vic_Aug_output_dir, well_data_cell[i,0], well_data_cell[i,1]))
	vic_sm[i] = temp[:,8] + temp[:,9] + temp[:,10]  # sm of three layers

#########################################################################
################### calculate correlation at each cell ##################
#########################################################################
well_sm_tau = np.empty([ncell, 5]) # ncell; [lat_grid; lon_grid; well depth (ft); tau; if significant (1 for sig and 0 for not-sig)]
for i in range(ncell):
	# select valid data
	data_valid = []  # nyear_valid; [well data; soil moisture]
	for j in range(nyear):
		if abs(well_data_cell[i,j+3]-missing_value)>0.001:  # if well data is not missing
			data_valid.append([well_data_cell[i,j+3], vic_sm[i,j]])
	data_valid = np.asarray(data_valid)
	n = np.shape(data_valid)[0]  # number of all valid pairs
	# calculate Kendall's tau
	P = 0  # number of concordant pairs
	M = 0  # number of inconcordant pairs
	for j in range(n-1):
		for k in range(j+1, n):
			if (data_valid[j,0]>data_valid[k,0] and data_valid[j,1]>data_valid[k,1]) or (data_valid[j,0]<data_valid[k,0] and data_valid[j,1]<data_valid[k,1]):  # if concordant
				P = P + 1
			elif (data_valid[j,0]>data_valid[k,0] and data_valid[j,1]<data_valid[k,1]) or (data_valid[j,0]<data_valid[k,0] and data_valid[j,1]>data_valid[k,1]):  # if discordant
				M = M + 1
	S = P - M
	tau = S / (n*(n-1)/2.0)
	# test if correlation is significant
#	# alpha = 0.05, 2-sided
#	if n==12:
#		crit_tau = 0.45
#	elif n==11:
#		crit_tau = 0.48
#	elif n==10:
#		crit_tau = 0.50
#	elif n==9:
#		crit_tau = 0.55
	# alpha = 0.1, 1-sided
	if n==12:
		crit_tau = 0.30
	elif n==11:
		crit_tau = 0.32
	elif n==10:
		crit_tau = 0.345
	elif n==9:
		crit_tau = 0.375
	if tau>=crit_tau:
		sig_flag = 1
	else:
		sig_flag = 0
	# put results in one matrix
	well_sm_tau[i,0] = well_data_cell[i,0]
	well_sm_tau[i,1] = well_data_cell[i,1]
	well_sm_tau[i,2] = well_data_cell[i,2]
	well_sm_tau[i,3] = tau
	well_sm_tau[i,4] = sig_flag

#######################################################################
############################# plotting ################################
#######################################################################
# plot scatter plot
fig = plt.figure()
plt.plot(well_sm_tau[:,2], well_sm_tau[:,3], 'o', markersize=5)
xx = range(int(np.min(well_sm_tau[:,2]))-1, int(np.max(well_sm_tau[:,2]))+1, 10)
plt.plot(xx, crit_tau*np.ones(len(xx)), 'r--')
plt.plot(xx, -crit_tau*np.ones(len(xx)), 'r--')
plt.xlabel('Well depth (feet)', fontsize=16)
plt.ylabel('kendall\'s tau', fontsize=16)
plt.title('Kendall\'s tau between Aug well obs. and soil moist. from VIC\nalpha = 0.1, 1-sided', fontsize=16)
fig.savefig('%s/scatter_tau_wellDepth_Aug.png' %plot_output_dir, format='png')

# plot scatter plot with LOWESS smoothing
yest = sm.nonparametric.lowess(well_sm_tau[:,3], well_sm_tau[:,2], frac=frac)
fig = plt.figure()
plt.plot(well_sm_tau[:,2], well_sm_tau[:,3], 'o', markersize=5)
plt.plot(yest[:,0], yest[:,1], 'k-', linewidth=3, label='LOWESS smooth, f=%.2f' %frac)
xx = range(int(np.min(well_sm_tau[:,2]))-1, int(np.max(well_sm_tau[:,2]))+1, 10)
plt.plot(xx, crit_tau*np.ones(len(xx)), 'r--', label='Confidence boundary')
plt.plot(xx, -crit_tau*np.ones(len(xx)), 'r--')
plt.xlabel('Well depth (feet)', fontsize=16)
plt.ylabel('kendall\'s tau', fontsize=16)
plt.legend(prop={'size':16})
plt.title('Kendall\'s tau between Aug well obs. and soil moist. from VIC\nalpha = 0.1, 1-sided', fontsize=16)
fig.savefig('%s/scatter_LOWESS_tau_wellDepth_Aug.png' %plot_output_dir, format='png')

exit()

# plot tau map, all cells
fig = plt.figure(figsize=(10,8))
ax = plt.axes([0, 0.08, 1, 0.75])
m = Basemap(projection='mill', llcrnrlat=24, urcrnrlat=50,\
            llcrnrlon=-120, urcrnrlon=-66, resolution='l')
m.drawcoastlines()
m.drawparallels(np.arange(-90., 91., 10.), labels=[True,True,False,False])
m.drawmeridians(np.arange(-180., 181., 10.), labels=[False,False,True,True])
m.drawmapboundary(fill_color='0.85')
m.fillcontinents(zorder=0, color='0.75')
m.drawcountries()
m.drawstates()

x, y = m(well_sm_tau[:,1], well_sm_tau[:,0])

cmap, norm, ticks, bounds = cmap_drought_listed(-1, 1)
for i in range(ncell):
	tau = well_sm_tau[i,3]
	for j in range(len(bounds)-1):
		if tau>=bounds[j] and tau<bounds[j+1]:
			break
		if np.absolute(tau-bounds[len(bounds)-1])<0.0001:
			i = len(bounds)-2
			break
	cc = cmap(j)
	ax.plot(x[i], y[i], color=cc, marker='s', markersize=7)

colormap = mpl.cm.ScalarMappable(norm, cmap)
colormap.set_array(range(-1,1))
cbar = plt.colorbar(colormap, orientation='horizontal', fraction=0.045, pad=0.06, ticks=ticks)
cbar.set_label('Kendall\'s tau', fontsize=14)
plt.text(0.5, 1.1, 'Kendall\'s tau between Aug well obs. and soil moist. from VIC\nAll grid cells', horizontalalignment='center', \
         fontsize=16, transform = ax.transAxes)
fig.savefig('%s/map_tau_all_sites_Aug.png' %plot_output_dir, format='png')

# plot tau map, only significent-correlation cells
fig = plt.figure(figsize=(10,8))
ax = plt.axes([0, 0.08, 1, 0.75])
m = Basemap(projection='mill', llcrnrlat=24, urcrnrlat=50,\
            llcrnrlon=-120, urcrnrlon=-66, resolution='l')
m.drawcoastlines()
m.drawparallels(np.arange(-90., 91., 10.), labels=[True,True,False,False])
m.drawmeridians(np.arange(-180., 181., 10.), labels=[False,False,True,True])
m.drawmapboundary(fill_color='0.85')
m.fillcontinents(zorder=0, color='0.75')
m.drawcountries()
m.drawstates()

x, y = m(well_sm_tau[:,1], well_sm_tau[:,0])

cmap, norm, ticks, bounds = cmap_drought_listed(-1, 1)
for i in range(ncell):
	tau = well_sm_tau[i,3]
	if tau>=crit_tau:
		for j in range(len(bounds)-1):
			if tau>=bounds[j] and tau<bounds[j+1]:
				break
			if np.absolute(tau-bounds[len(bounds)-1])<0.0001:
				i = len(bounds)-2
				break
		cc = cmap(j)
		ax.plot(x[i], y[i], color=cc, marker='s', markersize=7)

colormap = mpl.cm.ScalarMappable(norm, cmap)
colormap.set_array(range(-1,1))
cbar = plt.colorbar(colormap, orientation='horizontal', fraction=0.045, pad=0.06, ticks=ticks)
cbar.set_label('Kendall\'s tau', fontsize=14)
plt.text(0.5, 1.1, 'Kendall\'s tau between Aug well obs. and soil moist. from VIC\nOnly grid cells with significant correlation', horizontalalignment='center', \
         fontsize=16, transform = ax.transAxes)
fig.savefig('%s/map_tau_sigCorr_Aug.png' %plot_output_dir, format='png')



