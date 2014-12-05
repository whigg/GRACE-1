#!/usr/local/bin/python

#===============================================================================
#===============================================================================

def calc_7day_lowflow(daily_flow, start_year, end_year):
	'''Calculate annual 7-day low flow

	Input arguments:
		daily_flow: np array of daily flow 
			- [year] [month] [day] [flow (unit)] (flow<0 for missing or invalid data; date must be continuous)
		start_year: calculation start year
		end_year: calculation end year

	Return:
		an np array: [year] [month] [day] [7-day low flow (same unit as input)]

	Required packages:
		import numpy as np
		import datetime as dt

	Note:
		start_year and end_year must be within daily_flow data range;
		Calculate 7-day low flow between Apr 15 and Oct 31
	'''

	import numpy as np
	import datetime as dt

	############ check validity of input arguments ##########
	data_start_date = dt.date(int(daily_flow[0,0]), int(daily_flow[0,1]), int(daily_flow[0,2]))
	data_end_date = dt.date(int(daily_flow[-1,0]), int(daily_flow[-1,1]), int(daily_flow[-1,2]))
	# check if the data dates are continuous
	if ((data_end_date-data_start_date).days+1) != len(daily_flow):
		print 'Data dates are not continuous!'
		return None 
	# check if the calculation period is within the available data period
	if (dt.date(start_year, 4, 15) - data_start_date).days < 0:
		print 'Calculation start year earlier than available data!'
		return None
	elif (dt.date(end_year, 10, 31) - data_end_date).days > 0:
		print 'Calculation end year later than available data!'
		return None

	############ param set up #############
	nyear = end_year - start_year + 1
	yL7 = np.empty([nyear, 4]) # [year] [month] [day] [7-day low flow] (date is the first day)

	############ calculate 7-day low flow and corresponding month ##############
	for yr in range(start_year, end_year+1):  # loop through each year to be calculated
		first_day = (dt.date(yr,4,15)-data_start_date).days
		last_day = (dt.date(yr,10,31)-data_start_date).days
		yL7[yr-start_year][0] = yr
		# check the number of missing data; if >=50, omit the year
		count = 0
		for i in range(first_day, last_day+1):
			if daily_flow[i, 3]<0:
				count = count + 1
		if count>=50:
			yL7[yr-start_year][1] = -1
			yL7[yr-start_year][2] = -1
			yL7[yr-start_year][3] = -1
			continue
		# if missing data <50, calculate 7-day low flow
		min = 10000000
		for i in range(first_day, last_day-6):
			window = np.empty(7)
			count = 0
			for j in range(7):
				window[j] = daily_flow[i+j][3]
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
				low_month = daily_flow[i][1]  # keep track of the 7-day low flow date (the first day of the window)
				low_day = daily_flow[i][2]  # keep track of the 7-day low flow date (the first day of the window)
		yL7[yr-start_year][3] = min
		yL7[yr-start_year][1] = low_month
		yL7[yr-start_year][2] = low_day
	return yL7

#===============================================================================
#===============================================================================

def set_tick_fontsize(ax, fontsize=16, set_xaxis=True, set_yaxis=True):
	''' Set plot tick fontsize'''
	if set_xaxis==True:
		for tick in ax.xaxis.get_major_ticks():
			tick.label.set_fontsize(fontsize)
	if set_yaxis==True:
		for tick in ax.yaxis.get_major_ticks():
			tick.label.set_fontsize(fontsize)
	return ax



