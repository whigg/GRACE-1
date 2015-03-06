#!/usr/local/bin/python

import numpy as np
import scipy.stats
import os
import datetime as dt
import matplotlib as mpl
import matplotlib.pyplot as plt

############################# user defined ##############################
basin_list_path = '/usr1/ymao/other/GRACE/input/basin.list'
trend_resid_data_dir = '/usr1/ymao/other/GRACE/USGS_well_data/output/thiessenPoly_trend_sd_Aug_CRN' # [] [] [siteID] [lat] [lon] [trend (mm/yr)] [SD (mm)] [Thiessen polygon area (m^2)]
basin_value_ori_path = '/usr1/ymao/other/GRACE/input/basin.values.ori' # original basin values for plotting
output_dir = '/usr1/ymao/other/GRACE/USGS_well_data/output'

############################# load data #################################
basin_list = []
f = open(basin_list_path, 'r')
while 1:
	line = f.readline().rstrip('\n')
	if line=="":
		break
	basin_list.append(line)
f.close()
nbasin = len(basin_list)

trend_resid_data = []
for i in range(nbasin):
	trend_resid_data.append([])
	if os.path.isfile('%s/%s' %(trend_resid_data_dir,basin_list[i]))==True: # if file exist
		trend_resid_data[i] = np.loadtxt('%s/%s' %(trend_resid_data_dir,basin_list[i]), skiprows=1, usecols=(2,3,4,5,6,7))
		if np.shape(np.shape(trend_resid_data[i]))[0]==1:  # if the file only has single data line
			trend_resid_data[i] = trend_resid_data[i].reshape(1, np.shape(trend_resid_data[i])[0])

f = open(basin_value_ori_path, 'r')
basin_value_ori = []
while 1:
	line = f.readline().rstrip("\n")
	if line=="":
		break
	basin_value_ori.append(line)
f.close()

####################### Thiessen polygon average over each basin ###################
trend_ave = np.zeros(nbasin)
sd_resid_ave = np.zeros(nbasin)
for i in range(nbasin):
	if trend_resid_data[i]!=[]:
		print basin_list[i]
		nsite = np.shape(trend_resid_data[i])[0]
		print nsite
		tot_area = np.sum(trend_resid_data[i][:,5])
		for j in range(nsite):
			trend_ave[i] = trend_ave[i] + trend_resid_data[i][j,3]*trend_resid_data[i][j,5]
			sd_resid_ave[i] = sd_resid_ave[i] + trend_resid_data[i][j,4]*trend_resid_data[i][j,5]
		trend_ave[i] = trend_ave[i] / tot_area
		sd_resid_ave[i] = sd_resid_ave[i] / tot_area

# write basin storage trend and residual SD into files
f = open('%s/trend_mm_yr' %output_dir, 'w')
for i in range(nbasin):
	f.write('%s %.6f\n' %(basin_value_ori[i], trend_ave[i]))
f.close()

f = open('%s/sd_resid_mm' %output_dir, 'w')
for i in range(nbasin):
	f.write('%s %.6f\n' %(basin_value_ori[i], sd_resid_ave[i]))
f.close()


######################## simple average over each basin #######################
#trend_ave = np.zeros(nbasin)
#sd_resid_ave = np.zeros(nbasin)
#for i in range(nbasin):
#	if trend_resid_data[i]!=[]:
#		print basin_list[i]
#		nsite = np.shape(trend_resid_data[i])[0]
#		print nsite
#		for j in range(nsite):
#			trend_ave[i] = trend_ave[i] + trend_resid_data[i][j,3]
#			sd_resid_ave[i] = sd_resid_ave[i] + trend_resid_data[i][j,4]
#		trend_ave[i] = trend_ave[i] / nsite
#		sd_resid_ave[i] = sd_resid_ave[i] / nsite
#
## write basin storage trend and residual SD into files
#f = open('%s/trend_mm_yr' %output_dir, 'w')
#for i in range(nbasin):
#	f.write('%s %.6f\n' %(basin_value_ori[i], trend_ave[i]))
#f.close()
#
#f = open('%s/sd_resid_mm' %output_dir, 'w')
#for i in range(nbasin):
#	f.write('%s %.6f\n' %(basin_value_ori[i], sd_resid_ave[i]))
#f.close()








