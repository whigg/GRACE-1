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
site_info_dir = '/usr1/ymao/other/GRACE/USGS_well_data/site_info/site_info_have_data'  # Site identification number,Decimal latitude,Decimal longitude,Altitude of Gage/land surface [ft],Local aquifer type code,Field water-level measurements count
basin_list_path = '/usr1/ymao/other/GRACE/input/basin.list'  # list of 11 basin names
plots_output_dir = '/usr1/ymao/other/GRACE/USGS_well_data/plots'
siteID_list_climNet_path = '/usr1/ymao/other/GRACE/USGS_well_data/site_info/site_num_climNet_allSites'

sig_level = 0.05  # significance level for trend analysis; 2-sided

freq = 1  # final data frequency; unit: month
start_date = dt.datetime(year=2002, month=1, day=1)
end_date = dt.datetime(year=2013, month=12, day=31)
uni_window = 1 # check if there's a 'uni_window'-year window where there's no data; if there is such window, then not uniformly sampled; unit: year

nmonth = (end_date.year-start_date.year+1) * 12
nseg = nmonth / freq
# create a list of start day and end day of each segment
first_day_year = np.empty(nseg)
first_day_month = np.empty(nseg)
first_day_day = np.empty(nseg)

first_day = []
last_day = []

first_day_year[0] = start_date.year
first_day_month[0] = start_date.month
first_day_day[0] = start_date.day
first_day.append(dt.datetime(year=int(first_day_year[0]),month=int(first_day_month[0]),day=int(first_day_day[0])))
for i in range(1,nseg):
	if first_day_month[i-1]+freq<=12:
		first_day_year[i] = first_day_year[i-1]
		first_day_month[i] = first_day_month[i-1] + freq
		first_day_day[i] = 1
	else:
		first_day_year[i] = first_day_year[i-1] + 1
		first_day_month[i] = first_day_month[i-1] + freq - 12
		first_day_day[i] = 1
	first_day.append(dt.datetime(year=int(first_day_year[i]),month=int(first_day_month[i]),day=int(first_day_day[i])))
for i in range(nseg-1):
	last_day.append(first_day[i+1]-dt.timedelta(days=1))
last_day.append(end_date)

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
well_data = []   # basin; site; [siteID; date; water level below surface (ft); lat; lon; siteType; status_flag]
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
	# find surface height, lat, lon and siteType of the site
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
			for k in range(len(site_info[i])):  # find surface height, lat, lon and siteType of the site
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
					break

		if well_data_ori[i][j][7]=='':  # if there is no status code
			if well_data_ori[i][j][4]!='':  # if 'water level below surface' is available
				if flag!=0:
					print 'Warning1: not consistant at %s, %s' %(basin_list[i], siteID)
					well_data[i][count_site].append([well_data_ori[i][j][0],well_data_ori[i][j][1],surface_height-float(well_data_ori[i][j][4]), lat, lon, site_type, 0]) # convert 'below' to 'datum'
				else:
					well_data[i][count_site].append([well_data_ori[i][j][0],well_data_ori[i][j][1],float(well_data_ori[i][j][4]), lat, lon, site_type, 0])
			elif well_data_ori[i][j][5]!='':  # if 'water level above datum' is available, convert to 'water level below surface'
				if flag!=1:
					print 'Warning2: not consistent at %s, %s' %(basin_list[i], site_ID)
					well_data[i][count_site].append([well_data_ori[i][j][0],well_data_ori[i][j][1],surface_height-float(well_data_ori[i][j][5]), lat, lon, site_type, 0])  # convert 'datum' to 'below'
				else:
					well_data[i][count_site].append([well_data_ori[i][j][0],well_data_ori[i][j][1],float(well_data_ori[i][j][5]), lat, lon, site_type, 0])
		else:  # if there is status code
			well_data[i][count_site].append([well_data_ori[i][j][0],well_data_ori[i][j][1],well_data_ori[i][j][4], lat, lon, site_type, 1])

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

#################### plot histogram of the counts of data at all sites #################
#print '\nPlotting histogram...'
#count = []
#for i in range(nbasin): # for each basin
#	count.append([])
#	for j in range(len(well_data[i])):  # for each site
#		count[i].append([])
#		count[i][j] = 0
#		for k in range(len(well_data[i][j])):
#			date_str = well_data[i][j][k][1].split('/')
#			if len(date_str)==3:
#				date = dt.datetime(year=int(date_str[2]), month=int(date_str[0]), day=int(date_str[1]))
#				if date.year>=2002 and date.year<=2013:
#					count[i][j] = count[i][j] + 1
#
#count_tot = []
#for i in range(nbasin):
#	for j in range(len(well_data[i])):
#		count_tot.append(count[i][j])
#
#fig = plt.figure()
#plt.hist(count_tot, bins=50)
#plt.xlabel('Number of data at each site', fontsize=16)
#plt.ylabel('Frequency', fontsize=16)
#plt.title('Number of data at each site (Eastern U.s., 2002-2013)', fontsize=16)
#fig.savefig('%s/hist_data_count_tot.png' %plots_output_dir, format='png')
#
#fig = plt.figure()
#plt.hist(count_tot, bins=50, range=(0,200))
#plt.xlabel('Number of data at each site', fontsize=16)
#plt.ylabel('Frequency', fontsize=16)
#plt.title('Number of data at each site (Eastern U.S., 2002-2013)', fontsize=16)
#fig.savefig('%s/hist_data_count_tot_less200.png' %plots_output_dir, format='png')
#
#for i in range(nbasin):
#	fig = plt.figure()
#	plt.hist(count[i], bins=50)
#	plt.xlabel('Number of data at each site', fontsize=16)
#	plt.ylabel('Frequency', fontsize=16)
#	plt.title('Number of data at each site (%s, 2002-2013)' %basin_list[i], fontsize=16)
#	fig.savefig('%s/%s.hist_data_count.png' %(plots_output_dir,basin_list[i]), format='png')

################ select data (select the first data in every time segment) ####################
print '\nSelecting data...'
well_data_sel = []  # basin; site; [siteID; date; water level below surface (ft); lat; lon; siteType]
for i in range(nbasin):
	well_data_sel.append([])
	for j in range(len(well_data[i])):
		well_data_sel[i].append([])
		start_search_ind = 0
		for t in range(nseg):  # find the first data in this time segment
			for k in range(start_search_ind, len(well_data[i][j])):
				date_str = well_data[i][j][k][1].split('/')
				if len(date_str)==3:
					date = dt.datetime(year=int(date_str[2]), month=int(date_str[0]), day=int(date_str[1]))
					if (date-first_day[t]).days>=0 and (last_day[t]-date).days>=0: # if this date is within this time segment
						well_data_sel[i][j].append([well_data[i][j][k][0],well_data[i][j][k][1],well_data[i][j][k][2],well_data[i][j][k][3],well_data[i][j][k][4],well_data[i][j][k][5]])
						start_search_ind = k + 1
						break
					elif (date-last_day[t]).days>0: # if this date is after this time segment
						break

########################## check uniformity ##########################
print '\nChecking uniformity...'
well_data_uni = []
for i in range(nbasin):
	well_data_uni.append([])
	count_site = 0
	for j in range(len(well_data_sel[i])):
		flag = 0
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
				well_data_uni[i][count_site].append([well_data_sel[i][j][k][0],well_data_sel[i][j][k][1],well_data_sel[i][j][k][2],well_data_sel[i][j][k][3],well_data_sel[i][j][k][4],well_data_sel[i][j][k][5]])
			count_site = count_site + 1

for i in range(nbasin):
	print basin_list[i], len(well_data[i]), len(well_data_sel[i]), len(well_data_uni[i])


######################### calculate and plot trend - each site ###########################
print '\nCalculating and plotting trend at each site...'
well_trend = []  # basin; site; [siteID, lat, lon, trend (mm/yr), siteType, p-value]
for i in range(nbasin):
	well_trend.append([])
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
		w = np.linalg.lstsq(A.T,y)[0]  # trend: w[0]; unit: ft/day
		trend = w[0] * 12 * 25.4 * 365.25 # convert unit to: mm/yr
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
		# put all data in well_trend
		well_trend[i][j].append([well_data_uni[i][j][0][0],well_data_uni[i][j][0][3],well_data_uni[i][j][0][4],trend, well_data_uni[i][j][0][5], p_value])



## plot all sites
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
##		print basin_list[i], well_trend[i][j][0][0], well_trend[i][j][0][3]
#		x, y = m(well_trend[i][j][0][2], well_trend[i][j][0][1])
#		cs = plt.scatter(x, y, s=20, c=well_trend[i][j][0][3], cmap=cm.GMT_no_green_r, vmax=100, vmin=-100, linewidths=0)
#cbar = plt.colorbar(cs, fraction=0.045)
#cbar.set_label('Trend (mm/year)', fontsize=16)
#plt.title('Well data trend, 2002-2013, all sites', fontsize=16)
#
#fig.savefig('%s/well_trend_map_freq%dmon_window%dyear.png' %(plots_output_dir,freq,uni_window), format='png')
#
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
##		print basin_list[i], well_trend[i][j][0][0], well_trend[i][j][0][3]
#		x, y = m(well_trend[i][j][0][2], well_trend[i][j][0][1])
#		if well_trend[i][j][0][4]=='C':
#			cs = plt.scatter(x, y, s=20, c=well_trend[i][j][0][3], cmap=cm.GMT_no_green_r, vmax=100, vmin=-100, linewidths=0)
#cbar = plt.colorbar(cs, fraction=0.045)
#cbar.set_label('Trend (mm/year)', fontsize=16)
#plt.title('Well data trend, 2002-2013, confined sites', fontsize=16)
#
#fig.savefig('%s/well_trend_map_confined_freq%dmon_window%dyear.png' %(plots_output_dir,freq,uni_window), format='png')
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
##		print basin_list[i], well_trend[i][j][0][0], well_trend[i][j][0][3]
#		x, y = m(well_trend[i][j][0][2], well_trend[i][j][0][1])
#		if well_trend[i][j][0][4]=='U':
#			cs = plt.scatter(x, y, s=20, c=well_trend[i][j][0][3], cmap=cm.GMT_no_green_r, vmax=100, vmin=-100, linewidths=0)
#cbar = plt.colorbar(cs, fraction=0.045)
#cbar.set_label('Trend (mm/year)', fontsize=16)
#plt.title('Well data trend, 2002-2013, unconfined sites', fontsize=16)
#
#fig.savefig('%s/well_trend_map_unconfined_freq%dmon_window%dyear.png' %(plots_output_dir,freq,uni_window), format='png')
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
##		print basin_list[i], well_trend[i][j][0][0], well_trend[i][j][0][3]
#		x, y = m(well_trend[i][j][0][2], well_trend[i][j][0][1])
#		if well_trend[i][j][0][4]!='C' and well_trend[i][j][0][4]!='U':
#			cs = plt.scatter(x, y, s=20, c=well_trend[i][j][0][3], cmap=cm.GMT_no_green_r, vmax=100, vmin=-100, linewidths=0)
#cbar = plt.colorbar(cs, fraction=0.045)
#cbar.set_label('Trend (mm/year)', fontsize=16)
#plt.title('Well data trend, 2002-2013, unknown-type sites', fontsize=16)
#
#fig.savefig('%s/well_trend_map_unknown_freq%dmon_window%dyear.png' %(plots_output_dir,freq,uni_window), format='png')
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
##		print basin_list[i], well_trend[i][j][0][0], well_trend[i][j][0][3]
#		x, y = m(well_trend[i][j][0][2], well_trend[i][j][0][1])
#		if well_trend[i][j][0][5]>(1.0-sig_level/2.0) or well_trend[i][j][0][5]<sig_level/2.0:
#			cs = plt.scatter(x, y, s=20, c=well_trend[i][j][0][3], cmap=cm.GMT_no_green_r, vmax=100, vmin=-100, linewidths=0)
#cbar = plt.colorbar(cs, fraction=0.045)
#cbar.set_label('Trend (mm/year)', fontsize=16)
#plt.title('Well data trend, 2002-2013, all sites with significant trend', fontsize=16)
#
#fig.savefig('%s/well_trend_map_sig%.2f_freq%dmon_window%dyear.png' %(plots_output_dir,sig_level,freq,uni_window), format='png')


# plot only sites in Climate Response Network
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
		for k in range(len(siteID_list_climNet)):
			if int(well_trend[i][j][0][0])==int(siteID_list_climNet[k]):
				cs = plt.scatter(x, y, s=20, c=well_trend[i][j][0][3], cmap=cm.GMT_no_green_r, vmax=100, vmin=-100, linewidths=0)
				break
cbar = plt.colorbar(cs, fraction=0.045)
cbar.set_label('Trend (mm/year)', fontsize=16)
plt.title('Well data trend, field measurement, 2002-2013\nClimite Response Network sites', fontsize=16)

fig.savefig('%s/well_trend_map_ClimNet_fldMeas_freq%dmon_window%dyear.png' %(plots_output_dir,freq,uni_window), format='png')



