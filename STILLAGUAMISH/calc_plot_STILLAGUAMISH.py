#!/usr/local/bin/python

#===========================================================================
# Calculate baseflow and ground storage for STILLAGUAMISH site
# Yixin Mao
# 11/19/2014
#===========================================================================

import numpy as np
import datetime as dt
import my_functions
import matplotlib.pyplot as plt
from matplotlib.dates import DateFormatter

#########################################################
# User defined
#########################################################
flow_data_path = '/usr1/ymao/other/GRACE/STILLAGUAMISH/data_streamflow/12167000' # year, month, day, flow (cfs)
output_dir = '/usr1/ymao/other/GRACE/STILLAGUAMISH/output'
start_year = 1929 # calculation start year
end_year = 2014 # calculation end year

#########################################################
# import data
#########################################################
flow_data = np.loadtxt(flow_data_path)

#########################################################
# Calculate 7-day low flow
#########################################################
low_flow = my_functions.calc_7day_lowflow(flow_data, start_year, end_year) # [year] [month] [day] [7-day low flow (cfs)]

#########################################################
# plot 7-day low flow
#########################################################
fig = plt.figure(figsize=(12,5))
ax = plt.axes()
plt.plot(low_flow[:,0], low_flow[:,3], 'k-o')
plt.ylabel('7-day low flow (cfs)', fontsize=16)
plt.title('Stillaguamish River (USGS 12167000)', fontsize=16)
plt.xlim([start_year, end_year])
my_functions.set_tick_fontsize(ax, fontsize=16)
fig.savefig('%s/ts_y7L_12167000.png' %output_dir, format='png')

#########################################################
# plot low flow starting date
#########################################################
start_date = []  # the start date of low flow windows in each year; since we only need the month and day, we use a fake year (2011)
for i in range(len(low_flow)):
	start_date.append(dt.date(2011, int(low_flow[i,1]), int(low_flow[i,2])))

fig = plt.figure(figsize=(12,5))
ax = plt.axes()
plt.plot_date(low_flow[:,0], start_date, 'k-o', xdate=False, ydate=True)
plt.ylabel('Starting date of 7-day low flow', fontsize=16)
plt.title('Stillaguamish River (USGS 12167000)', fontsize=16)
plt.xlim([start_year, end_year])
my_functions.set_tick_fontsize(ax, fontsize=16)
ax.yaxis.set_major_formatter(DateFormatter('%b %d'))
fig.savefig('%s/ts_y7L_date_12167000.png' %output_dir, format='png')

#########################################################
# plot histogram of low flow month
#########################################################
fig = plt.figure()
ax = plt.axes()
plt.hist(low_flow[:,1])
plt.xlabel('Low flow month', fontsize=16)
plt.ylabel('Frequncy', fontsize=16)
my_functions.set_tick_fontsize(ax, fontsize=16)
fig.savefig('%s/hist_y7L_month_12167000.png' %output_dir, format='png')









