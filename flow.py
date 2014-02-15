#!/usr/local/bin/python

###########################################
# HEADER
###########################################
import numpy as np
import scipy.stats
from os import listdir
import datetime as dt
import matplotlib as mpl
import matplotlib.pyplot as plt
import pdb

##########################################
# user defined
##########################################
data_path = "/usr1/ymao/GRACE/data_streamflow/streamflow_1988_2012"  # data unit: cfs
start_year = 2003 
end_year = 2012
K = 45   # characteristic time scale; unit: days
area_index = 1  # 0 for using basin area by Brutsaert; 1 for using drainage area by USGS

##########################################
# load data
##########################################
station = {}
for filename in listdir(data_path):
	station[filename] = np.loadtxt("%s/%s" %(data_path, filename))

area_Bru = np.loadtxt("/usr1/ymao/GRACE/data_streamflow/station_area_Bru")
area_USGS = np.loadtxt("/usr1/ymao/GRACE/data_streamflow/station_area_Bru")

station_order = np.loadtxt("/usr1/ymao/GRACE/station.order")

#########################################
# calculate 7-day low flow
# yL7[std][year][index] contains 7-day low flow and corresponding month
# for each year at each station
# if data is missing, then flow and month is -1 
#########################################
yL7 = {}
for std in station:
	yL7[std] = np.empty((end_year-start_year+1, 2))

for std in (station):
	for yr in range(start_year, end_year+1):
		first_day = (dt.date(yr,4,15)-dt.date(1988,1,1)).days
		last_day = (dt.date(yr,10,31)-dt.date(1988,1,1)).days
		# check the number of missing data; if >=50, omit the year
		count = 0
		for i in range(first_day, last_day+1):
			if station[std][i][3]<0:
				count = count + 1
		if count>=50:
			yL7[std][yr-start_year][0] = -1
			yL7[std][yr-start_year][1] = -1
			continue
		# if missing data <50, calculate 7-day low flow
		min = 10000000
		for i in range(first_day, last_day-6):
			window = np.empty(7)
			count = 0
			for j in range(7):
				window[j] = station[std][i+j][3]
				if window[j]<0: 
					count = count + 1
			if count==0:  # if no missing data in the window
				ave = np.mean(window)
			elif count==1: # if only one missing data in the window
				sum = 0
				for j in range(7):
					if window[j]>=0:
						sum = sum + window[j]
				ave = sum/6
			elif count>1:  # if more than one data is missing in the window, omit
				continue
			if ave<min:
				min = ave 
				low_month = station[std][i+3][1]  # keep track of the 7-day low flow month
		yL7[std][yr-start_year][0] = min
		yL7[std][yr-start_year][1] = low_month

#######################################
# calculate temporal trends by simple linear regression
#######################################
trend_flow = {}
flow = {}
c = {}
m = {}
for std in station:
	# copy the non-missing station data and corresponding years into a new array
	count = 0
	for i in range(np.shape(yL7[std])[0]):
		if np.absolute(yL7[std][i][0]+1)>0.0001:  # if not missing
			count = count + 1
	flow[std] = np.empty((count,2))
	count = 0
	for i in range(np.shape(yL7[std])[0]):
		if np.absolute(yL7[std][i][0]+1)>0.0001:  # if not missing
			flow[std][count][0] = start_year + i
			flow[std][count][1] = yL7[std][i][0]
			count = count + 1
	# calculate trend by simple linear regression
	A = np.array([flow[std][:,0], np.ones(np.shape(flow[std])[0])])
	m[std], c[std] = np.linalg.lstsq(A.T, flow[std][:,1])[0]    # y = mx + c 
	trend_flow[std] = m[std]   # flow trend in cfs/yr

trend_storage = {}
if area_index==0:   # using basin area reported by Brutsaert
	area = area_Bru
else:   # using drainage reported by USGS
	area = area_USGS

for i in range(np.shape(area)[0]):
	std = '0' + str(int(area[i][0]))
	trend_storage[std] = trend_flow[std] * K / area[i][1]  # storage trend in (cfs*day)/(yr*km2)
	trend_storage[std] = trend_storage[std] * np.power(30.48,3) * 86400 / np.power(10,9) # storage trend in mm/yr

# calculate p value (two-tailed)
p = {}
for std in station:
	SSE = np.sum(np.power(flow[std][:,1]-(c[std]+m[std]*flow[std][:,0]), 2))
	s_2 = SSE / (np.shape(flow[std])[0] - 2)
	sB1 = np.sqrt(s_2 / np.sum(np.power(flow[std][:,0]-np.mean(flow[std][:,0]), 2)))
	t = trend_flow[std] / sB1
	tcdf = scipy.stats.t.cdf(t, np.shape(flow[std])[0])
	if t<0:  # if trend<0
		p[std] = tcdf * 2
	else:   # if trend.=0
		p[std] = (1-tcdf) * 2

# write storage trend and corresponding p value into file
f = open('/usr1/ymao/GRACE/output/trend_%s_%s_K%d' %(start_year,end_year,K), 'w')
for i in range(np.shape(station_order)[0]):
	std = '0' + str(int(station_order[i]))
	if p[std]>=0.005:
		f.write("%s %.5f %.2f\n" %(std, trend_storage[std], p[std]))
	else:
		f.write("%s %.5f %.1g\n" %(std, trend_storage[std], p[std]))
f.close()

## write yL7  (41 columns, each column is the yL7 for all years at one basin)
# # convert cfs to mm/d
#yL7_mm_d = {}
#for std in station:
#	yL7_mm_d[std] = np.empty(end_year-start_year+1)
#
#for i in range(np.shape(area)[0]):
#	std = '0' + str(int(area[i][0]))
#	for y in range(end_year-start_year+1):
#		if np.absolute(yL7[std][y][0]+1)<0.0001:   # if missing year
#			yL7_mm_d[std][y] = -1
#		else:  # if not missing
#			yL7_mm_d[std][y] = yL7[std][y][0] / area[i][1]   # cfs/km2
#			yL7_mm_d[std][y] = yL7_mm_d[std][y] * np.power(30.48,3) * 86400 / np.power(10,9) # mm/d 
#
# # write to file
#f = open('/usr1/ymao/GRACE/output/yL7_%s_%s' %(start_year,end_year), 'w')
#for y in range(end_year-start_year+1):
#	for i in range(np.shape(station_order)[0]):
#		std = '0' + str(int(station_order[i]))
#		if np.absolute(yL7[std][y][0]+1)<0.0001:  # if missing year
#			f.write("%d " %yL7[std][y][0])
#		else:  # if not missing
#			f.write("%f " %yL7_mm_d[std][y]) 
#	f.write("\n")
#f.close()
#
#
#
