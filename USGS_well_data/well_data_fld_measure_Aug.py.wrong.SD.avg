#!/usr/local/bin/python

import numpy as np
import datetime as dt
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap, cm
from collections import Counter
from scipy.stats import norm

######################## user defined #######################
well_data_dir = '/usr1/ymao/other/GRACE/USGS_well_data/ori_data_from_Liz'  # SiteID,Date level measured,Time level measured,Time Datum,Water level value(ft) below land surface,Water level value (feet) above specific vertical datum,Referenced vertical datum, The status of the site at the time the water level was measured
site_info_dir = '/usr1/ymao/other/GRACE/USGS_well_data/site_info/site_info_fld_have_data_v2'  # Site identification number,Decimal latitude,Decimal longitude,Altitude of Gage/land surface,Local aquifer type code,National aquifer code,Well depth,Field water-level measurements count
basin_list_path = '/usr1/ymao/other/GRACE/input/basin.list'  # list of 11 basin names
plots_output_dir = '/usr1/ymao/other/GRACE/USGS_well_data/plots'
output_dir = '/usr1/ymao/other/GRACE/USGS_well_data/output'
siteID_list_climNet_path = '/usr1/ymao/other/GRACE/USGS_well_data/site_info/site_num_climNet_allSites'
Sy_list_path = '/usr1/ymao/other/GRACE/USGS_well_data/site_info/Sy_national_aquifer_code' # [three digits (in national aquifer code)] [Sy]

sig_level = 0.05  # significance level for trend analysis; 2-sided
missing_value = -9999

freq = 1  # final data frequency; unit: month
start_date = dt.datetime(year=2002, month=1, day=1)
end_date = dt.datetime(year=2013, month=12, day=31)
uni_window = 2 # check if there's a 'uni_window'-year window where there's no data; if there is such window, then not uniformly sampled; unit: year

start_year = start_date.year
end_year = end_date.year
nyear = end_year - start_year + 1

basin_list = []
f = open(basin_list_path, 'r')
while 1:
	line = f.readline().rstrip('\n')
	if line=="":
		break
	basin_list.append(line)
f.close()
nbasin = len(basin_list)

siteID_list_climNet = []
f = open(siteID_list_climNet_path, 'r')
while 1:
	line = f.readline().rstrip('\n')
	if line=="":
		break
	siteID_list_climNet.append(line)
f.close()

Sy_list = []
f = open(Sy_list_path, 'r')
while 1:
	line = f.readline().rstrip('\n')
	if line=="":
		break
	Sy_list.append(line.split())
f.close()

######################## load data #########################
print 'Loading data...'

# load well observation data
print 'Loading well observation data...'
well_data_ori = []  # nbasin; ndata; data
for i in range(nbasin):
#	print 'Loading basin %d...' %(i+1)
	well_data_ori.append([])
	f = open('%s/%s.csv' %(well_data_dir, basin_list[i]), 'r')
	count = 0
	while 1:
		line = f.readline().rstrip('\n').rstrip('\r')
		if line=="":
			break
		if count!=0:
			well_data_ori[i].append(line.split(','))
		count = count + 1
	f.close()

# load well site info
print 'Loading site info...'
site_info = []  # nbasin; site_info
for i in range(nbasin):
	site_info.append([])
	f = open('%s/%s.csv' %(site_info_dir, basin_list[i]), 'r')
	while 1:
		line = f.readline().rstrip('\n').rstrip('\r')
		if line=="":
			break
		site_info[i].append(line.split(','))
	f.close()

################## process data - check if data type is consistent for each site ################
print 'Processing data...'
well_data = []   # basin; site; [siteID; date; water level above a datum (ft); lat; lon; siteType; status_flag; national aquifer code; well depth (ft)]
for i in range(nbasin): 
	# initialize
	well_data.append([])
	siteID = well_data_ori[i][0][0]
	well_data[i].append([])
	count_site = 0
	# check consistency of data
	if well_data_ori[i][0][4]!='':
		flag = 0  # 'water level below surface' exists
	elif well_data_ori[i][0][5]!='':
		flag = 1 # 'water level above datum' exists
	# find surface height, lat, lon, siteType, aquifer code and well depth of the site
	for k in range(len(site_info[i])): 
		if site_info[i][k][0]==siteID:
			if site_info[i][k][3]!='':
				surface_height = float(site_info[i][k][3])
			else:
				surface_height = ''
			if site_info[i][k][1]!='':
				lat = float(site_info[i][k][1])
			else:
				lat = ''
			if site_info[i][k][2]!='':
				lon = float(site_info[i][k][2])
			else:
				lon = ''
			site_type = site_info[i][k][4]
			national_aquifer_code = site_info[i][k][5]
			if  site_info[i][k][6]!='':
				well_depth = float(site_info[i][k][6])
			else:
				well_depth = ''
			break

	for j in range(len(well_data_ori[i])):
		if well_data_ori[i][j][0]!=siteID:  # if change to the next site
			siteID = well_data_ori[i][j][0]
			well_data[i].append([])
			count_site = count_site + 1
			if well_data_ori[i][j][4]!='':
				flag = 0  # 'water level below surface' exists
			elif well_data_ori[i][j][5]!='':
				flag = 1 # 'water level above datum' exists
			for k in range(len(site_info[i])):  # find surface height, lat, lon, siteType, national aquifer code and well depth of the site
				if site_info[i][k][0]==siteID:
					if site_info[i][k][3]!='':
						surface_height = float(site_info[i][k][3])
					else:
						surface_height = ''
					if site_info[i][k][1]!='':
						lat = float(site_info[i][k][1])
					else:
						lat = ''
					if site_info[i][k][2]!='':
						lon = float(site_info[i][k][2])
					else:
						lon = ''
					site_type = site_info[i][k][4]
					national_aquifer_code = site_info[i][k][5]
					if  site_info[i][k][6]!='':
						well_depth = float(site_info[i][k][6])
					else:
						well_depth = ''
					break

		if well_data_ori[i][j][7]=='':  # if there is no status code
			if well_data_ori[i][j][4]!='':  # if 'water level below surface' is available
				if flag!=0:
					print 'Warning1: not consistant at %s, %s' %(basin_list[i], siteID)
				well_data[i][count_site].append([well_data_ori[i][j][0],well_data_ori[i][j][1],surface_height-float(well_data_ori[i][j][4]), lat, lon, site_type, 0, national_aquifer_code, well_depth]) # convert 'below' to 'datum'
			elif well_data_ori[i][j][5]!='':  # if 'water level above datum' is available
				if flag!=1:
					print 'Warning2: not consistent at %s, %s' %(basin_list[i], site_ID)
				well_data[i][count_site].append([well_data_ori[i][j][0],well_data_ori[i][j][1],float(well_data_ori[i][j][5]), lat, lon, site_type, 0, national_aquifer_code, well_depth])  # above datum
		else:  # if there is status code
			well_data[i][count_site].append([well_data_ori[i][j][0],well_data_ori[i][j][1],well_data_ori[i][j][4], lat, lon, site_type, 1, national_aquifer_code, well_depth])

	# delete invalid sites
	for jj in range(len(well_data[i])-1,-1,-1):
		# delete no-valid-data site
		if len(well_data[i][jj])==0: 
			del well_data[i][jj]
		# delete no-lat-lon site
		elif well_data[i][jj][0][3]=='' or well_data[i][jj][0][4]=='': 
			del well_data[i][jj]
		# delete sites with status code in any observation
		else: 
			for k in range(len(well_data[i][jj])):
				if well_data[i][jj][k][6]==1: # if there is a status code at this observation
					del well_data[i][jj]
					break

################ select data (select the first August data in each year) ####################
print '\nSelecting data...'
well_data_sel = []  # basin; site; [siteID; date; water level above a datum (ft); lat; lon; siteType;national aquifer code; well depth (ft)]
for i in range(nbasin):
	well_data_sel.append([])
	for j in range(len(well_data[i])):
		well_data_sel[i].append([])
		start_search_ind = 0
		for t in range(nyear):  # find the first August data in each year
			for k in range(start_search_ind, len(well_data[i][j])):
				date_str = well_data[i][j][k][1].split('/')
				if len(date_str)==3:
					date = dt.datetime(year=int(date_str[2]), month=int(date_str[0]), day=int(date_str[1]))
					if date.year==start_year+t and date.month==8: # if this date is in August
						well_data_sel[i][j].append([well_data[i][j][k][0],well_data[i][j][k][1],well_data[i][j][k][2],well_data[i][j][k][3],well_data[i][j][k][4],well_data[i][j][k][5],well_data[i][j][k][7], well_data[i][j][k][8]])
						start_search_ind = k + 1
						break
					elif (date-dt.datetime(year=start_year+t,month=8,day=31)).days>0: # if this date is after this August
						break

########################## check uniformity ##########################
print '\nChecking uniformity...'
well_data_uni = []
for i in range(nbasin):
	well_data_uni.append([])
	count_site = 0
	for j in range(len(well_data_sel[i])):
		if len(well_data_sel[i][j])>0:
			flag = 0

			date_str_first_data = well_data_sel[i][j][0][1].split('/')
			if len(date_str_first_data)==3:
				date_first_data = dt.datetime(year=int(date_str_first_data[2]), month=int(date_str_first_data[0]), day=int(date_str_first_data[1]))
				if (date_first_data-start_date).days>uni_window*365:  # if the first day is too late
					flag = 1

			date_str_last_data = well_data_sel[i][j][len(well_data_sel[i][j])-1][1].split('/')
			if len(date_str_last_data)==3:
				date_last_data = dt.datetime(year=int(date_str_last_data[2]), month=int(date_str_last_data[0]), day=int(date_str_last_data[1]))
				if (end_date-date_last_data).days>uni_window*365:  # if the last day is too early
					flag = 1

			for k in range(len(well_data_sel[i][j])-1):
				date_str = well_data_sel[i][j][k][1].split('/')
				date_str_next = well_data_sel[i][j][k+1][1].split('/')
				if len(date_str)==3 and len(date_str_next)==3:
					date = dt.datetime(year=int(date_str[2]), month=int(date_str[0]), day=int(date_str[1]))
					date_next = dt.datetime(year=int(date_str_next[2]), month=int(date_str_next[0]), day=int(date_str_next[1]))
					if (date_next-date).days>uni_window*365:  # if more then the window-size data is missing, throw out this site
						flag = 1
						break
			if flag==0:
				well_data_uni[i].append([])
				for k in range(len(well_data_sel[i][j])):
					well_data_uni[i][count_site].append([well_data_sel[i][j][k][0],well_data_sel[i][j][k][1],well_data_sel[i][j][k][2],well_data_sel[i][j][k][3],well_data_sel[i][j][k][4],well_data_sel[i][j][k][5],well_data_sel[i][j][k][6],well_data_sel[i][j][k][7]])
				count_site = count_site + 1

for i in range(nbasin):
	print basin_list[i], len(well_data[i]), len(well_data_sel[i]), len(well_data_uni[i])

################################## plot well depth ################################
#fig = plt.figure(figsize=(10,10))
#ax = fig.add_axes([0.1,0.1,0.8,0.8])
#m = Basemap(llcrnrlon=-120., llcrnrlat=20., urcrnrlon=-60., urcrnrlat=50., rsphere=(6378137.00,6356752.3142), resolution='l', area_thresh=1000.,projection='lcc', lat_1=50.,lon_0=-107.,ax=ax)
#m.drawcoastlines()
#m.drawparallels(np.arange(-90., 91., 5.), labels=[1,0,0,1])
#m.drawmeridians(np.arange(-180., 181., 5.), labels=[1,0,0,1])
#m.drawmapboundary(fill_color='0.85')
#m.fillcontinents(zorder=0, color='0.75')
#m.drawcountries()
#m.drawstates()
#
#for i in range(nbasin):
#	if len(well_data_uni[i])>0:
#		for j in range(len(well_data_uni[i])):
#			x, y = m(well_data_uni[i][j][0][4], well_data_uni[i][j][0][3])
#			cs = plt.scatter(x, y, s=20, c=well_data_uni[i][j][0][7], cmap='jet', vmax=200, vmin=0, linewidths=0)
#cbar = plt.colorbar(cs, fraction=0.045)
#cbar.set_label('Well depth (feet)', fontsize=16)
#plt.title('Well depth, field measurement sites', fontsize=16)
#
#fig.savefig('%s/well_depth_Aug_freq%dmon_window%dyear.png' %(plots_output_dir,freq,uni_window), format='png')


######################### calculate trend - each site ###########################
print '\nCalculating and plotting trend at each site...\n'
well_trend = []  # basin; site; [siteID, lat, lon, trend(mm/yr), siteType, p-value, Sy, SD_residual(mm), well_depth(ft)]
for i in range(nbasin):
	well_trend.append([])
	if len(well_data_uni[i])>0:
		for j in range(len(well_data_uni[i])):
			well_trend[i].append([])
			# calculate trend
			data_temp = np.empty([len(well_data_uni[i][j]), 2])
			for k in range(len(well_data_uni[i][j])):
				date_str = well_data_uni[i][j][k][1].split('/')
				date = dt.datetime(year=int(date_str[2]), month=int(date_str[0]), day=int(date_str[1]))
				date_ind = (date-start_date).days 
				data_temp[k,0] = date_ind # unit: day
				data_temp[k,1] = well_data_uni[i][j][k][2]  # unit: ft
			xi = data_temp[:,0]
			A = np.array([xi, np.ones(np.shape(xi)[0])])
			y = data_temp[:,1] 
			m, c = np.linalg.lstsq(A.T,y)[0]  # y = mx + c; trend: m (ft/day); intercept: c (ft)
			trend = m * 12 * 25.4 * 365.25 # convert unit to: mm/yr
			# calculate standard deviation of residuals
			residual = np.empty(len(well_data_uni[i][j]))
			for k in range(len(well_data_uni[i][j])):
				date_str = well_data_uni[i][j][k][1].split('/')
				date = dt.datetime(year=int(date_str[2]), month=int(date_str[0]), day=int(date_str[1]))
				date_ind = (date-start_date).days
				wl_est = m * date_ind + c  # ft
				residual[k] = well_data_uni[i][j][k][2] - wl_est  # ft
			residual = residual * 12 * 25.4  # convert to mm
			std_resid = np.std(residual)  # mm
			# calculate p-value using Mann-Kendall test
			n = np.shape(y)[0]
			s = 0
			for ii in range(n-1):
				for jj in range(ii+1,n):
					if y[jj]>y[ii]:
						s = s + 1
					elif y[jj]<y[ii]:
						s = s - 1
			lst = Counter(y).most_common()
			if lst[0][1]==1:  # if no repeating values
				var_s = (n*(n-1)*(2*n+5)) / 18.0
			else:  # if there are repeating values
				tt = []
				for ii in range(len(lst)):
					tt.append(lst[ii][1])
				tie = Counter(tt).most_common()  # first column: tie index (i); second column: number of this tie (ti)
				var_s = n*(n-1)*(2*n+5)
				for ii in range(len(tie)):
					ti = tie[ii][1]
					tie_i = tie[ii][0]
					var_s = var_s - ti * tie_i * (tie_i-1) * (2*tie_i+5)
				var_s = var_s / 18.0
			sigma_s = np.sqrt(var_s)
			if s>0:
				z = float((s-1)) / sigma_s
			elif s<0:
				z = float((s+1)) / sigma_s
			else:
				z = 0
			p_value = norm.cdf(z)
			# get Sy
			national_aquifer_code = well_data_uni[i][j][0][6]
			flag = 0
			for kk in range(len(Sy_list)):
				if national_aquifer_code[1:4]==Sy_list[kk][0]:
					Sy = float(Sy_list[kk][1])
					flag = 1
					break
			if flag==0:  # if no valid aquifer type, estimate Sy as 0.1
				Sy = 0.1
			# put all data in well_trend
			well_trend[i][j].append([well_data_uni[i][j][0][0],well_data_uni[i][j][0][3],well_data_uni[i][j][0][4],trend, well_data_uni[i][j][0][5], p_value, Sy, std_resid, well_data_uni[i][j][0][7]])

# write trend and residual results into files; one file for each basin, for Thiessen polygon basin average
for i in range(nbasin):
	f = open('%s/%s.txt' %(output_dir,basin_list[i]), 'w')
	f.write('siteID,lat,lon,trend_mm_yr,sd_mm\n')
	if len(well_trend[i])>0:
		for j in range(len(well_trend[i])):
			Sy = well_trend[i][j][0][6]
			f.write('%s,%.6f,%.6f,%.4f,%.4f\n' %(well_trend[i][j][0][0], well_trend[i][j][0][1], well_trend[i][j][0][2], well_trend[i][j][0][3]*Sy, well_trend[i][j][0][7]*Sy))  # siteID, lat, lon, storage trend (mm/yr), storage SD of residual (mm)
	f.close()

exit()

# writ trend, SD of residual and well depth into one file
f = open('%s/trend_sd_Aug_FLD' %(output_dir), 'w')
for i in range(nbasin):
	if len(well_trend[i])>0:
		for j in range(len(well_trend[i])):
			Sy = well_trend[i][j][0][6]
			if well_trend[i][j][0][8]!='':
				f.write('%s %.6f %.6f %.4f %.4f %.1f\n' %(well_trend[i][j][0][0], well_trend[i][j][0][1], well_trend[i][j][0][2], well_trend[i][j][0][3]*Sy, well_trend[i][j][0][7]*Sy, well_trend[i][j][0][8]))  # siteID, lat, lon, storage trend (mm/yr), storage SD of residual (mm), well depth (ft)
f.close()

# write time series into one file
f = open('%s/ts_Aug_FLD' %output_dir, 'w')
for i in range(nbasin):
	if len(well_data_uni[i])>0:
		for j in range(len(well_data_uni[i])):
			flag = 0
			for kk in range(len(Sy_list)):
				if well_data_uni[i][j][0][6][1:4]==Sy_list[kk][0]:
					Sy = float(Sy_list[kk][1])
					flag = 1
					break
			if flag==0:  # if no valid aquifer type, estimate Sy as 0.1
				Sy = 0.1
			if well_data_uni[i][j][0][7]!='':
				f.write('%s %.6f %.6f %.1f ' %(well_data_uni[i][j][0][0], well_data_uni[i][j][0][3], well_data_uni[i][j][0][4], well_data_uni[i][j][0][7])) # siteID; lat; lon; well depth (ft)
				ts = np.ones(12) * missing_value
				for k in range(len(well_data_uni[i][j])):
					date_str = well_data_uni[i][j][k][1].split('/')
					date = dt.datetime(year=int(date_str[2]), month=int(date_str[0]), day=int(date_str[1]))
					year_ind = date.year - start_year
					ts[year_ind] = well_data_uni[i][j][k][2] *12 * 25.4 * Sy  # storage, mm
				for y in range(nyear):
					if ts[y]!=missing_value:
						f.write('%.1f ' %ts[y])
					else:
						f.write('%d ' %ts[y])
				f.write('\n')
f.close()

################################# plot storage trend ################################
# plot all sites
fig = plt.figure(figsize=(10,10))
ax = fig.add_axes([0.1,0.1,0.8,0.8])
m = Basemap(llcrnrlon=-120., llcrnrlat=20., urcrnrlon=-60., urcrnrlat=50., rsphere=(6378137.00,6356752.3142), resolution='l', area_thresh=1000.,projection='lcc', lat_1=50.,lon_0=-107.,ax=ax)
m.drawcoastlines()
m.drawparallels(np.arange(-90., 91., 5.), labels=[1,0,0,1])
m.drawmeridians(np.arange(-180., 181., 5.), labels=[1,0,0,1])
m.drawmapboundary(fill_color='0.85')
m.fillcontinents(zorder=0, color='0.75')
m.drawcountries()
m.drawstates()

for i in range(nbasin):
	for j in range(len(well_trend[i])):
		x, y = m(well_trend[i][j][0][2], well_trend[i][j][0][1])
		cs = plt.scatter(x, y, s=20, c=well_trend[i][j][0][3]*well_trend[i][j][0][6], cmap=cm.GMT_no_green_r, vmax=20, vmin=-20, linewidths=0)
cbar = plt.colorbar(cs, fraction=0.045)
cbar.set_label('Trend (mm/year)', fontsize=16)
plt.title('Aug GW Storage trend, field measurement, \n2002-2013, all sites', fontsize=16)

fig.savefig('%s/storage_Aug_freq%dmon_window%dyear.png' %(plots_output_dir,freq,uni_window), format='png')

## plot sites <= 100 feet
#fig = plt.figure(figsize=(10,10))
#ax = fig.add_axes([0.1,0.1,0.8,0.8])
#m = Basemap(llcrnrlon=-120., llcrnrlat=20., urcrnrlon=-60., urcrnrlat=50., rsphere=(6378137.00,6356752.3142), resolution='l', area_thresh=1000.,projection='lcc', lat_1=50.,lon_0=-107.,ax=ax)
#m.drawcoastlines()
#m.drawparallels(np.arange(-90., 91., 5.), labels=[1,0,0,1])
#m.drawmeridians(np.arange(-180., 181., 5.), labels=[1,0,0,1])
#m.drawmapboundary(fill_color='0.85')
#m.fillcontinents(zorder=0, color='0.75')
#m.drawcountries()
#m.drawstates()
#
#for i in range(nbasin):
#	for j in range(len(well_trend[i])):
#		x, y = m(well_trend[i][j][0][2], well_trend[i][j][0][1])
#		if well_trend[i][j][0][8]<=100:
#			cs = plt.scatter(x, y, s=20, c=well_trend[i][j][0][3]*well_trend[i][j][0][6], cmap=cm.GMT_no_green_r, vmax=20, vmin=-20, linewidths=0)
#cbar = plt.colorbar(cs, fraction=0.045)
#cbar.set_label('Trend (mm/year)', fontsize=16)
#plt.title('Aug GW Storage trend, field measurement, \n2002-2013, well depth shallower than 100 ft', fontsize=16)
#
#fig.savefig('%s/storage_Aug_depthLT100_freq%dmon_window%dyear.png' %(plots_output_dir,freq,uni_window), format='png')

## plot sites > 100 feet
#fig = plt.figure(figsize=(10,10))
#ax = fig.add_axes([0.1,0.1,0.8,0.8])
#m = Basemap(llcrnrlon=-120., llcrnrlat=20., urcrnrlon=-60., urcrnrlat=50., rsphere=(6378137.00,6356752.3142), resolution='l', area_thresh=1000.,projection='lcc', lat_1=50.,lon_0=-107.,ax=ax)
#m.drawcoastlines()
#m.drawparallels(np.arange(-90., 91., 5.), labels=[1,0,0,1])
#m.drawmeridians(np.arange(-180., 181., 5.), labels=[1,0,0,1])
#m.drawmapboundary(fill_color='0.85')
#m.fillcontinents(zorder=0, color='0.75')
#m.drawcountries()
#m.drawstates()
#
#for i in range(nbasin):
#	for j in range(len(well_trend[i])):
#		x, y = m(well_trend[i][j][0][2], well_trend[i][j][0][1])
#		if well_trend[i][j][0][8]>100:
#			cs = plt.scatter(x, y, s=20, c=well_trend[i][j][0][3]*well_trend[i][j][0][6], cmap=cm.GMT_no_green_r, vmax=20, vmin=-20, linewidths=0)
#cbar = plt.colorbar(cs, fraction=0.045)
#cbar.set_label('Trend (mm/year)', fontsize=16)
#plt.title('Aug GW Storage trend, field measurement, \n2002-2013, well depth deeper than 100 ft', fontsize=16)
#
#fig.savefig('%s/storage_Aug_depthGT100_freq%dmon_window%dyear.png' %(plots_output_dir,freq,uni_window), format='png')


## plot confined sites
#fig = plt.figure(figsize=(10,10))
#ax = fig.add_axes([0.1,0.1,0.8,0.8])
#m = Basemap(llcrnrlon=-120., llcrnrlat=20., urcrnrlon=-60., urcrnrlat=50., rsphere=(6378137.00,6356752.3142), resolution='l', area_thresh=1000.,projection='lcc', lat_1=50.,lon_0=-107.,ax=ax)
#m.drawcoastlines()
#m.drawparallels(np.arange(-90., 91., 5.), labels=[1,0,0,1])
#m.drawmeridians(np.arange(-180., 181., 5.), labels=[1,0,0,1])
#m.drawmapboundary(fill_color='0.85')
#m.fillcontinents(zorder=0, color='0.75')
#m.drawcountries()
#m.drawstates()
#
#for i in range(nbasin):
#	for j in range(len(well_trend[i])):
#		x, y = m(well_trend[i][j][0][2], well_trend[i][j][0][1])
#		if well_trend[i][j][0][4]=='C':
#			cs = plt.scatter(x, y, s=20, c=well_trend[i][j][0][3]*well_trend[i][j][0][6], cmap=cm.GMT_no_green_r, vmax=20, vmin=-20, linewidths=0)
#cbar = plt.colorbar(cs, fraction=0.045)
#cbar.set_label('Trend (mm/year)', fontsize=16)
#plt.title('Aug GW storage trend, field measurement, \n2002-2013, confined sites', fontsize=16)
#
#fig.savefig('%s/storage_Aug_confined_freq%dmon_window%dyear.png' %(plots_output_dir,freq,uni_window), format='png')
#
## plot unconfined sites
#fig = plt.figure(figsize=(10,10))
#ax = fig.add_axes([0.1,0.1,0.8,0.8])
#m = Basemap(llcrnrlon=-120., llcrnrlat=20., urcrnrlon=-60., urcrnrlat=50., rsphere=(6378137.00,6356752.3142), resolution='l', area_thresh=1000.,projection='lcc', lat_1=50.,lon_0=-107.,ax=ax)
#m.drawcoastlines()
#m.drawparallels(np.arange(-90., 91., 5.), labels=[1,0,0,1])
#m.drawmeridians(np.arange(-180., 181., 5.), labels=[1,0,0,1])
#m.drawmapboundary(fill_color='0.85')
#m.fillcontinents(zorder=0, color='0.75')
#m.drawcountries()
#m.drawstates()
#
#for i in range(nbasin):
#	for j in range(len(well_trend[i])):
#		x, y = m(well_trend[i][j][0][2], well_trend[i][j][0][1])
#		if well_trend[i][j][0][4]=='U':
#			cs = plt.scatter(x, y, s=20, c=well_trend[i][j][0][3]*well_trend[i][j][0][6], cmap=cm.GMT_no_green_r, vmax=20, vmin=-20, linewidths=0)
#cbar = plt.colorbar(cs, fraction=0.045)
#cbar.set_label('Trend (mm/year)', fontsize=16)
#plt.title('Aug GW storage trend, field measurement, \n2002-2013, unconfined sites', fontsize=16)
#
#fig.savefig('%s/storage_Aug_unconfined_freq%dmon_window%dyear.png' %(plots_output_dir,freq,uni_window), format='png')
#
## plot unknown-type sites
#fig = plt.figure(figsize=(10,10))
#ax = fig.add_axes([0.1,0.1,0.8,0.8])
#m = Basemap(llcrnrlon=-120., llcrnrlat=20., urcrnrlon=-60., urcrnrlat=50., rsphere=(6378137.00,6356752.3142), resolution='l', area_thresh=1000.,projection='lcc', lat_1=50.,lon_0=-107.,ax=ax)
#m.drawcoastlines()
#m.drawparallels(np.arange(-90., 91., 5.), labels=[1,0,0,1])
#m.drawmeridians(np.arange(-180., 181., 5.), labels=[1,0,0,1])
#m.drawmapboundary(fill_color='0.85')
#m.fillcontinents(zorder=0, color='0.75')
#m.drawcountries()
#m.drawstates()
#
#for i in range(nbasin):
#	for j in range(len(well_trend[i])):
#		x, y = m(well_trend[i][j][0][2], well_trend[i][j][0][1])
#		if well_trend[i][j][0][4]!='C' and well_trend[i][j][0][4]!='U':
#			cs = plt.scatter(x, y, s=20, c=well_trend[i][j][0][3]*well_trend[i][j][0][6], cmap=cm.GMT_no_green_r, vmax=20, vmin=-20, linewidths=0)
#cbar = plt.colorbar(cs, fraction=0.045)
#cbar.set_label('Trend (mm/year)', fontsize=16)
#plt.title('Aug GW storage trend, 2002-2013, \nunknown-type sites', fontsize=16)
#
#fig.savefig('%s/storage_Aug_unknown_freq%dmon_window%dyear.png' %(plots_output_dir,freq,uni_window), format='png')
#
#
## plot only significant trends, all sites
#fig = plt.figure(figsize=(10,10))
#ax = fig.add_axes([0.1,0.1,0.8,0.8])
#m = Basemap(llcrnrlon=-120., llcrnrlat=20., urcrnrlon=-60., urcrnrlat=50., rsphere=(6378137.00,6356752.3142), resolution='l', area_thresh=1000.,projection='lcc', lat_1=50.,lon_0=-107.,ax=ax)
#m.drawcoastlines()
#m.drawparallels(np.arange(-90., 91., 5.), labels=[1,0,0,1])
#m.drawmeridians(np.arange(-180., 181., 5.), labels=[1,0,0,1])
#m.drawmapboundary(fill_color='0.85')
#m.fillcontinents(zorder=0, color='0.75')
#m.drawcountries()
#m.drawstates()
#
#for i in range(nbasin):
#	for j in range(len(well_trend[i])):
#		x, y = m(well_trend[i][j][0][2], well_trend[i][j][0][1])
#		if well_trend[i][j][0][5]>(1.0-sig_level/2.0) or well_trend[i][j][0][5]<sig_level/2.0:
#			cs = plt.scatter(x, y, s=20, c=well_trend[i][j][0][3]*well_trend[i][j][0][6], cmap=cm.GMT_no_green_r, vmax=20, vmin=-20, linewidths=0)
#cbar = plt.colorbar(cs, fraction=0.045)
#cbar.set_label('Trend (mm/year)', fontsize=16)
#plt.title('Aug GW storage trend, 2002-2013, \nall sites with significant trend', fontsize=16)
#
#fig.savefig('%s/storage_Aug_sig%.2f_freq%dmon_window%dyear.png' %(plots_output_dir,sig_level,freq,uni_window), format='png')
#
## plot only sites in Climate Response Network
#fig = plt.figure(figsize=(10,10))
#ax = fig.add_axes([0.1,0.1,0.8,0.8])
#m = Basemap(llcrnrlon=-120., llcrnrlat=20., urcrnrlon=-60., urcrnrlat=50., rsphere=(6378137.00,6356752.3142), resolution='l', area_thresh=1000.,projection='lcc', lat_1=50.,lon_0=-107.,ax=ax)
#m.drawcoastlines()
#m.drawparallels(np.arange(-90., 91., 5.), labels=[1,0,0,1])
#m.drawmeridians(np.arange(-180., 181., 5.), labels=[1,0,0,1])
#m.drawmapboundary(fill_color='0.85')
#m.fillcontinents(zorder=0, color='0.75')
#m.drawcountries()
#m.drawstates()
#
#for i in range(nbasin):
#	for j in range(len(well_trend[i])):
#		x, y = m(well_trend[i][j][0][2], well_trend[i][j][0][1])
#		for k in range(len(siteID_list_climNet)):
#			if int(well_trend[i][j][0][0])==int(siteID_list_climNet[k]):
#				cs = plt.scatter(x, y, s=20, c=well_trend[i][j][0][3]*well_trend[i][j][0][6], cmap=cm.GMT_no_green_r, vmax=20, vmin=-20, linewidths=0)
#				break
#cbar = plt.colorbar(cs, fraction=0.045)
#cbar.set_label('Trend (mm/year)', fontsize=16)
#plt.title('Aug GW storage trend, field measurement, \n2002-2013, Climite Response Network sites', fontsize=16)
#
#fig.savefig('%s/storage_Aug_ClimNet_fldMeas_freq%dmon_window%dyear.png' %(plots_output_dir,freq,uni_window), format='png')


########################### plot SD of residulas ################################
# plot all sites
fig = plt.figure(figsize=(10,10))
ax = fig.add_axes([0.1,0.1,0.8,0.8])
m = Basemap(llcrnrlon=-120., llcrnrlat=20., urcrnrlon=-60., urcrnrlat=50., rsphere=(6378137.00,6356752.3142), resolution='l', area_thresh=1000.,projection='lcc', lat_1=50.,lon_0=-107.,ax=ax)
m.drawcoastlines()
m.drawparallels(np.arange(-90., 91., 5.), labels=[1,0,0,1])
m.drawmeridians(np.arange(-180., 181., 5.), labels=[1,0,0,1])
m.drawmapboundary(fill_color='0.85')
m.fillcontinents(zorder=0, color='0.75')
m.drawcountries()
m.drawstates()

for i in range(nbasin):
	for j in range(len(well_trend[i])):
		x, y = m(well_trend[i][j][0][2], well_trend[i][j][0][1])
		cs = plt.scatter(x, y, s=20, c=well_trend[i][j][0][7]*well_trend[i][j][0][6], cmap=cm.GMT_no_green_r, vmin=0, vmax=100, linewidths=0)
cbar = plt.colorbar(cs, fraction=0.045)
cbar.set_label('Standard deviation (mm)', fontsize=16)
plt.title('Standard deviation in residual Aug GW, \nfield measurement, 2002-2013, all sites', fontsize=16)

fig.savefig('%s/sd_resid_Aug_freq%dmon_window%dyear.png' %(plots_output_dir,freq,uni_window), format='png')


## plot sites <=100 ft
#fig = plt.figure(figsize=(10,10))
#ax = fig.add_axes([0.1,0.1,0.8,0.8])
#m = Basemap(llcrnrlon=-120., llcrnrlat=20., urcrnrlon=-60., urcrnrlat=50., rsphere=(6378137.00,6356752.3142), resolution='l', area_thresh=1000.,projection='lcc', lat_1=50.,lon_0=-107.,ax=ax)
#m.drawcoastlines()
#m.drawparallels(np.arange(-90., 91., 5.), labels=[1,0,0,1])
#m.drawmeridians(np.arange(-180., 181., 5.), labels=[1,0,0,1])
#m.drawmapboundary(fill_color='0.85')
#m.fillcontinents(zorder=0, color='0.75')
#m.drawcountries()
#m.drawstates()
#
#for i in range(nbasin):
#	for j in range(len(well_trend[i])):
#		x, y = m(well_trend[i][j][0][2], well_trend[i][j][0][1])
#		if well_trend[i][j][0][8]<=100:
#			cs = plt.scatter(x, y, s=20, c=well_trend[i][j][0][7]*well_trend[i][j][0][6], cmap='jet_r', vmin=0, vmax=100, linewidths=0)
#cbar = plt.colorbar(cs, fraction=0.045)
#cbar.set_label('Standard deviation (mm)', fontsize=16)
#plt.title('Standard deviation in residual Aug GW, \nfield measurement, 2002-2013, sites shallower than 100 ft', fontsize=16)
#
#fig.savefig('%s/sd_resid_Aug_depthLT100_freq%dmon_window%dyear.png' %(plots_output_dir,freq,uni_window), format='png')
#
## plot sites >100 ft
#fig = plt.figure(figsize=(10,10))
#ax = fig.add_axes([0.1,0.1,0.8,0.8])
#m = Basemap(llcrnrlon=-120., llcrnrlat=20., urcrnrlon=-60., urcrnrlat=50., rsphere=(6378137.00,6356752.3142), resolution='l', area_thresh=1000.,projection='lcc', lat_1=50.,lon_0=-107.,ax=ax)
#m.drawcoastlines()
#m.drawparallels(np.arange(-90., 91., 5.), labels=[1,0,0,1])
#m.drawmeridians(np.arange(-180., 181., 5.), labels=[1,0,0,1])
#m.drawmapboundary(fill_color='0.85')
#m.fillcontinents(zorder=0, color='0.75')
#m.drawcountries()
#m.drawstates()
#
#for i in range(nbasin):
#	for j in range(len(well_trend[i])):
#		x, y = m(well_trend[i][j][0][2], well_trend[i][j][0][1])
#		if well_trend[i][j][0][8]>100:
#			cs = plt.scatter(x, y, s=20, c=well_trend[i][j][0][7]*well_trend[i][j][0][6], cmap='jet_r', vmin=0, vmax=100, linewidths=0)
#cbar = plt.colorbar(cs, fraction=0.045)
#cbar.set_label('Standard deviation (mm)', fontsize=16)
#plt.title('Standard deviation in residual Aug GW, \nfield measurement, 2002-2013, sites deeper than 100 ft', fontsize=16)
#
#fig.savefig('%s/sd_resid_Aug_depthGT100_freq%dmon_window%dyear.png' %(plots_output_dir,freq,uni_window), format='png')

############################## plot time series in some regions ########################
#print 'Calculating anamolies...\n'
## calculate anomalies
#well_data_anom = []
#for i in range(nbasin):
#	well_data_anom.append([])
#	for j in range(len(well_data_uni[i])):
#		well_data_anom[i].append([])
#		# calculate mean water level for this site
#		ave_water_level = 0
#		for k in range(len(well_data_uni[i][j])):
#			ave_water_level = ave_water_level + well_data_uni[i][j][k][2]
#		ave_water_level = ave_water_level / len(well_data_uni[i][j])
#		# calculate anamaly for this site
#		for k in range(len(well_data_uni[i][j])):
#			anom = well_data_uni[i][j][k][2] - ave_water_level # water level
#			if well_trend[i][j][0][6]!='':
#				well_data_anom[i][j].append([well_data_uni[i][j][k][0],well_data_uni[i][j][k][1],anom*well_trend[i][j][0][6],well_data_uni[i][j][k][3],well_data_uni[i][j][k][4],well_data_uni[i][j][k][5]])  # anamoly of water storage
#
#
###### plot southern part of Florida ######
#ave_anom_tm = np.zeros(nyear)
#count = np.zeros(nyear)
#for i in range(nbasin):
#	for j in range(len(well_data_anom[i])):
#		if len(well_data_anom[i][j])>0:
#			lat = well_data_anom[i][j][0][3]
#			lon = well_data_anom[i][j][0][4]
#			if lat>25 and lat<30 and lon>-85 and lon<-80:  # if in the region
#				for k in range(len(well_data_anom[i][j])):
#					date_str = well_data_uni[i][j][k][1].split('/')
#					date = dt.datetime(year=int(date_str[2]), month=int(date_str[0]), day=int(date_str[1]))
#					for t in range(nyear):
#						if date.year==start_year+t and date.month==8: # if in this time segment
#							ave_anom_tm[t] = ave_anom_tm[t] + well_data_anom[i][j][k][2]
#							count[t] = count[t] + 1
#							break
#ave_anom_tm = ave_anom_tm / count  # unit: ft
#ave_anom_tm = ave_anom_tm * 12 * 25.4  # unit: mm
#
#fig = plt.figure()
#year = []
#for y in range(start_year, end_year+1):
#	year.append(dt.datetime(year=y, month=8, day=1))
#plt.plot_date(year, ave_anom_tm)
#plt.plot_date(year, ave_anom_tm, 'b-')
#
#xi = []
#for t in range(nyear):
#	xi.append((dt.datetime(year=start_year+t,month=8,day=1)-start_date).days)
#A = np.array([xi, np.ones(np.shape(xi)[0])])
#y = ave_anom_tm
#w = np.linalg.lstsq(A.T,y)[0] # y = w[0]* x + w[1]; x: days
#x_plot = []
#y_plot = []
#for t in range(nyear):
#	x_plot.append(start_date+dt.timedelta(days=xi[t]))
#	y_plot.append(w[0]*xi[t]+w[1])
#plt.plot_date(x_plot, y_plot, 'k--')
#plt.xlabel('Year', fontsize=16)
#plt.ylabel('Water level anomaly (mm)', fontsize=16)
#plt.title('Average Aug storage anomaly, southern part of Florida', fontsize=16)
#fig.savefig('%s/ts_Aug_anom_fld_allSites_freq%dmon_southern_florida.png' %(plots_output_dir, freq), format='png')
#
###### plot North Dakota ######
#ave_anom_tm = np.zeros(nyear)
#count = np.zeros(nyear)
#for i in range(nbasin):
#	for j in range(len(well_data_anom[i])):
#		if len(well_data_anom[i][j])>0:
#			lat = well_data_anom[i][j][0][3]
#			lon = well_data_anom[i][j][0][4]
#			if lat>45 and lat<50 and lon>-105 and lon<-95:  # if in the region
#				for k in range(len(well_data_anom[i][j])):
#					date_str = well_data_uni[i][j][k][1].split('/')
#					date = dt.datetime(year=int(date_str[2]), month=int(date_str[0]), day=int(date_str[1]))
#					for t in range(nyear):
#						if date.year==start_year+t and date.month==8: # if in this time segment
#							ave_anom_tm[t] = ave_anom_tm[t] + well_data_anom[i][j][k][2]
#							count[t] = count[t] + 1
#							break
#ave_anom_tm = ave_anom_tm / count  # unit: ft
#ave_anom_tm = ave_anom_tm * 12 * 25.4  # unit: mm
#
#fig = plt.figure()
#year = []
#for y in range(start_year, end_year+1):
#	year.append(dt.datetime(year=y, month=8, day=1))
#plt.plot_date(year, ave_anom_tm)
#plt.plot_date(year, ave_anom_tm, 'b-')
#
#xi = []
#for t in range(nyear):
#	xi.append((dt.datetime(year=start_year+t,month=8,day=1)-start_date).days)
#A = np.array([xi, np.ones(np.shape(xi)[0])])
#y = ave_anom_tm
#w = np.linalg.lstsq(A.T,y)[0] # y = w[0]* x + w[1]; x: days
#x_plot = []
#y_plot = []
#for t in range(nyear):
#	x_plot.append(start_date+dt.timedelta(days=xi[t]))
#	y_plot.append(w[0]*xi[t]+w[1])
#plt.plot_date(x_plot, y_plot, 'k--')
#plt.xlabel('Year', fontsize=16)
#plt.ylabel('Water level anomaly (mm)', fontsize=16)
#plt.title('Average Aug storage anomaly, North Dakota', fontsize=16)
#fig.savefig('%s/ts_Aug_anom_fld_allSites_freq%dmon_north_dakota.png' %(plots_output_dir, freq), format='png')
#
###### plot Virginia ######
#ave_anom_tm = np.zeros(nyear)
#count = np.zeros(nyear)
#for i in range(nbasin):
#	for j in range(len(well_data_anom[i])):
#		if len(well_data_anom[i][j])>0:
#			lat = well_data_anom[i][j][0][3]
#			lon = well_data_anom[i][j][0][4]
#			if lat>35 and lat<39 and lon>-80 and lon<-75:  # if in the region
#				for k in range(len(well_data_anom[i][j])):
#					date_str = well_data_uni[i][j][k][1].split('/')
#					date = dt.datetime(year=int(date_str[2]), month=int(date_str[0]), day=int(date_str[1]))
#					for t in range(nyear):
#						if date.year==start_year+t and date.month==8: # if in this time segment
#							ave_anom_tm[t] = ave_anom_tm[t] + well_data_anom[i][j][k][2]
#							count[t] = count[t] + 1
#							break
#ave_anom_tm = ave_anom_tm / count  # unit: ft
#ave_anom_tm = ave_anom_tm * 12 * 25.4  # unit: mm
#
#fig = plt.figure()
#year = []
#for y in range(start_year, end_year+1):
#	year.append(dt.datetime(year=y, month=8, day=1))
#plt.plot_date(year, ave_anom_tm)
#plt.plot_date(year, ave_anom_tm, 'b-')
#
#xi = []
#for t in range(nyear):
#	xi.append((dt.datetime(year=start_year+t,month=8,day=1)-start_date).days)
#A = np.array([xi, np.ones(np.shape(xi)[0])])
#y = ave_anom_tm
#w = np.linalg.lstsq(A.T,y)[0] # y = w[0]* x + w[1]; x: days
#x_plot = []
#y_plot = []
#for t in range(nyear):
#	x_plot.append(start_date+dt.timedelta(days=xi[t]))
#	y_plot.append(w[0]*xi[t]+w[1])
#plt.plot_date(x_plot, y_plot, 'k--')
#plt.xlabel('Year', fontsize=16)
#plt.ylabel('Water level anomaly (mm)', fontsize=16)
#plt.title('Average Aug storage anomaly, Virginia', fontsize=16)
#fig.savefig('%s/ts_Aug_anom_fld_allSites_freq%dmon_virginia.png' %(plots_output_dir, freq), format='png')





