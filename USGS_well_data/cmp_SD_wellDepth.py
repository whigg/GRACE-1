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
data_path = '/usr1/ymao/other/GRACE/USGS_well_data/output/trend_sd_Aug_CRN_and_FLD' # [siteID] [lat] [lon] [storage trend (mm/yr)] [storage SD of residual (mm)] [well depth (ft)]; This is one file for all sites
plot_output_dir = '/usr1/ymao/other/GRACE/USGS_well_data/plots'
missing_value = -9999
nyear = 12
frac = 0.65  # f factor in LOWESS smothing

#########################################################################
############################# load data #################################
#########################################################################
data = np.loadtxt(data_path)
nsite = np.shape(data)[0]

#######################################################################
############################# plotting ################################
#######################################################################
# plot scatter plot
fig = plt.figure()
plt.plot(data[:,5], data[:,4], 'o', markersize=5)
plt.xlabel('Well depth (feet)', fontsize=16)
plt.ylabel('Standard deviation (mm)', fontsize=16)
plt.title('Aug GW standard deviation of trend residual v.s. well depth\n CRN and FLD', fontsize=16)
fig.savefig('%s/scatter_SD_wellDepth_Aug.png' %plot_output_dir, format='png')

# plot scatter plot with LOWESS smoothing
yest = sm.nonparametric.lowess(data[:,4], data[:,5], frac=frac)
fig = plt.figure()
plt.plot(data[:,5], data[:,4], 'o', markersize=5)
plt.plot(yest[:,0], yest[:,1], 'k-', linewidth=3, label='LOWESS smooth, f=%.2f' %frac)
plt.xlabel('Well depth (feet)', fontsize=16)
plt.ylabel('Standard deviation (mm)', fontsize=16)
plt.legend(prop={'size':16})
plt.title('Aug GW standard deviation of trend residual v.s. well depth\n CRN and FLD', fontsize=16)
fig.savefig('%s/scatter_LOWESS_SD_wellDepth_Aug.png' %plot_output_dir, format='png')

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



