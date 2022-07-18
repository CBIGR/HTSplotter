import numpy as np
import os
from grupping import Groupping
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm
from matplotlib import ticker
from matplotlib.backends.backend_pdf import PdfPages
import math
import scipy.interpolate as inter
from mpl_toolkits.mplot3d import Axes3D
from scipy.optimize import curve_fit
from doseresponsecurvefunctionscommon import DoseResponseFunctions


class ResponseCurve:

	def __init__(self, time_point, popt, pcov, compound, zip=False):

		self.popt = popt
		self.pcov = pcov

		self.m = popt[0]
		self.lambd = popt[1]

		if not zip:
			if popt[0] < 0:
				self.Emin = 0
			else:
				self.Emin = popt[2]
			self.m = popt[0]
			self.lambd = popt[1]
			self.Emax = popt[3]

		self.time_point = time_point
		self.compoundname = compound.grup.name
		self.concentrationrange = compound.concentration_range


	def four_parameter_logistic(self, x):

		# for LOEWE
		# b = Emax
		# l = Emin
		# x0 = m
		# k = lambda

		func = self.Emax + ((self.Emin - self.Emax) / (1 + ((x / self.m) ** self.lambd)))

		return func

class DoseResponse:

	# def __init__(self, data, std, grup, row, concentration, confluency, timeposition, timepoint, txt_path, pdf_pages,
	# 			 readout, readout_units):
	def __init__(self, combination):
		# self.data_range = None
		# self.std_range = None
		# self.inhib_data = None
		# self.inhib_std = None
		self.colun_num = 3

		self.grup = None
		self.concentration_range = None  # combination.concentration_range
		self.response_range = None  # combination.confluency_range
		self.response_rangestd = None
		self.confluency_range = None
		self.std_range = None

		# From combination script, compoundscreen or genetic-chemical
		self.row = combination.row_num
		self.time_position = combination.time_position
		self.time_selected = combination.time_selected
		self.time_selectedgr = combination.time_selectedgr
		self.time_positiongr = combination.position_gr
		self.alltimepoints = combination.elapse
		self.txt_path = combination.txt_path
		self.pdf_pages = combination.pdf_pages
		self.pdf_path = combination.pdf_path
		self.readout = combination.readout
		self.readout_units = combination.readout_units


		# get from doseresponse script
		self.curves = []
		self.curveszip = []

		# self.curve_specifictimepoints = []
		self.compoundname = []
		self.compoundconcentration = []
		self.compounds_grups = []

		# Atributes for on and off situation
		self.compounds_list_on = None
		self.compounds_list_off = None
		self.compounds_name_on_off = None
		self.dose_curve_on_off = None
		self.data_curve_on_off = None
		self.std_curve_on_off = None
		self.concentration_range_on_off = None
		self.ec_concentration_on_off = None
		self.ic_concentration_on_off = None
		self.concentration_value_on_off = None

	def get_curve_alltimepoint(self, combination, data, grup):
		self.grup = grup
		self.response_range = np.asarray(data)
		self.concentration_range = np.asarray(combination.concentration_range)
		# self.concentrationtesting = np.asarray(combination.concentrationtotest)
		self.compoundname.append(self.grup.name)
		self.compoundconcentration.append(self.grup.concentration[:, 0])
		self.curves.append([])
		# self.alltimepoints = combination.time_selected
		for tim in range(len(self.alltimepoints)):
			print(self.alltimepoints[tim])
			# get the data based on index position from concentration
			self.curves[-1].append(self.get_curve_fit(tim, self.alltimepoints[tim]))
		# self.curves[-1].append(
			# 	self.get_curve_fit(self.concentration_range, self.response_range[tim], self.alltimepoints[tim]))  # four-parameter logistic
		print('Normalcurvefit')

	def get_curveforzip_alltimepoint(self, combination, data, grup):
		self.grup = grup
		self.response_range = np.asarray(data[:, 1:])
		self.concentration_range = np.asarray(combination.concentration_range[1:])
		self.curveszip = []

		for tim in range(len(self.alltimepoints)):
			print(self.alltimepoints[tim])
			# get the data based on index position from concentration
			self.curveszip.append(self.get_curve_fit(tim, self.alltimepoints[tim], True))

		print('getcurveZIP')

	def get_curve_specifictimepoints(self, combination, data, grup):
		# this function gets the curves for specific time points, normally requested for
		# compound combination scripts and genetic-chemical perturbagem
		self.grup = grup
		self.response_range = data
		self.concentration_range = combination.concentration_range
		# self.compoundname.append(self.grup[0].name)
		# self.compoundconcentration.append(self.grup[0].concentration[:, 0])
		self.curves = []

		for tim in range(len(self.time_selected)):
			print(self.time_position[tim])
			# get the data based on index position from concentration
			self.curves.append(
				self.get_curve_fit(self.time_position[tim], self.time_selected[tim]))  # four-parameter logistic
		print('jweon')

	def get_curve_fit(self, tim, tp, zip=False):

		# if zip:
		# 	p1 = [np.median(self.concentration_range), 0.05]
		# else:
		if '100' in str(self.response_range[tim, :][0]):
			p1 = [np.median(self.concentration_range), 0.05, max(self.response_range[tim, 1:]),
				  min(self.response_range[tim, 1:])]
		else:
			p1 = [np.median(self.concentration_range), 0.05, min(self.response_range[tim, 1:]),
				  max(self.response_range[tim, 1:])]

		try:
			# if zip:
			# 	popt, pcov = curve_fit(self.four_parameter_logistic, self.concentration_range,
			# 						   self.response_range[tim, :], p1, method='lm') # ,
			# 						   # bounds=((np.min(self.concentration_range), 0),
			# 						   # 	   (0.90*np.max(self.concentration_range), np.inf)))

			popt, pcov = curve_fit(self.four_parameter_logistic, self.concentration_range[1:],
								   self.response_range[tim, 1:], p1, method='lm')

		except RuntimeError as error:

			# if zip:
			# 	popt, pcov = curve_fit(self.four_parameter_logistic, self.concentration_range,
			# 						   self.response_range[tim, :], p1, method='lm')  # ,
			# 						   # bounds=((np.min(self.concentration_range)), (np.max(self.concentration_range))),
			# 						   # maxfev=int(1e7))
			# else:
			popt, pcov = curve_fit(self.four_parameter_logistic, self.concentration_range[1:],
								   self.response_range[tim, 1:], p1, method='lm', maxfev=int(1e7))

		return ResponseCurve(tp, popt, pcov, self, zip)

	def get_curve_specifictimepointsGR(self, combination, data, grup):
		# this function gets the curves for specific time points, normally requested for
		# compound combination scripts and genetic-chemical perturbagem
		self.grup = grup
		self.response_range = data
		# self.response_range = data
		self.concentration_range = np.asarray(combination.concentration_range)
		# self.compoundname.append(self.grup[0].name)
		# self.compoundconcentration.append(self.grup[0].concentration[:, 0])
		self.curves = []

		for tim in range(len(self.time_selectedgr)):
			print(self.time_position[tim])
			# get the data based on index position from concentration
			self.curves.append(
				self.get_curve_fitGR(tim, self.time_selected[tim]))  # four-parameter logistic
		print('Gr dosecurve')

	def get_curve_fitGR(self, tim, tp):

		p1 = [np.median(self.concentration_range[1:]), 0.05, max(self.response_range[tim, 1:]),
			  min(self.response_range[tim, 1:])]

		try:
			# if zip:
			# 	popt, pcov = curve_fit(self.four_parameter_logistic, self.concentration_range,
			# 						   self.response_range[tim, :], p1, method='lm') # ,
			# 						   # bounds=((np.min(self.concentration_range), 0),
			# 						   # 	   (0.90*np.max(self.concentration_range), np.inf)))

			popt, pcov = curve_fit(self.four_parameter_logistic, self.concentration_range[1:],
								   self.response_range[tim, 1:], p1, method='lm')

		except RuntimeError as error:

			# if zip:
			# 	popt, pcov = curve_fit(self.four_parameter_logistic, self.concentration_range,
			# 						   self.response_range[tim, :], p1, method='lm')  # ,
			# 						   # bounds=((np.min(self.concentration_range)), (np.max(self.concentration_range))),
			# 						   # maxfev=int(1e7))
			# else:
			popt, pcov = curve_fit(self.four_parameter_logistic, self.concentration_range[1:],
								   self.response_range[tim, 1:], p1, method='lm', maxfev=int(1e7))

		return ResponseCurve(tp, popt, pcov, self)

	def get_confluency_range(self, compound, grup, grothrate=0, synergy=False):

		data_range = []  # this value is for DMSO control
		std_range = []

		std = compound.stdnormalized_perc
		time = list(range(len(self.alltimepoints)))
		controlinformation = compound.control_info[-1]
		if synergy:
			data = compound.inhib_data
		else:
			data = compound.normalized_perc

			# time = list(range(len(self.alltimepoints)))

		# time = compound.time_position
		controlinfo = compound.control_info
		# concentration_range = compound.concentration_range
		grup = grup
		# for k in range(len(self.time_position)):

		for k in range(len(time)):

			tmp1 = [data[controlinformation][time[k]]]  # this value is for DMSO control
			tmp2 = [std[controlinformation][time[k]]]
			for j in range(len(grup.concentration)):
				# if grup.concentration[j][0] == concentration_range[j + 1]:
				pos = int(grup.concentration[j][-1])
				# tmp1.append(data[time[k]][j])
				tmp1.append(data[pos][time[k]])
				tmp2.append(std[pos][time[k]])
			data_range.append(tmp1)
			std_range.append(tmp2)

		self.confluency_range = np.asarray(data_range)
		# self.concentration_range = np.asarray(concentration_range)
		self.std_range = np.asarray(std_range)

		return self.confluency_range, self.std_range

	def get_confluency_rangeGR(self, compound, grup, posi, posj, geneticchemical=0):

		data_range = []  # this value is for DMSO control
		std_range = []

		std = compound.stdinh
		controlinformation = compound.control_info[-1]
		data = compound.grresults
		controlinfo = compound.control_info
		grup = grup
		for k in range(len(self.time_selectedgr)):
			posgr = self.time_positiongr[k]
			tmp1 = [data[posi][posj][0][0][posgr]]  # this value is for DMSO control
			for j in range(len(grup.concentration)):
				tmp1.append(data[posi][posj][j+1][0][posgr])
			data_range.append(tmp1)
		for k in range(len(self.time_selectedgr)):
			stdpo = self.time_position[k]
			tmp2 = [std[controlinformation][stdpo]]
			for j in range(len(grup.concentration)):
				pos = int(grup.concentration[j][-1])
				tmp2.append(std[pos][stdpo])
			std_range.append(tmp2)
		self.confluency_rangeGR = np.asarray(data_range)
		self.std_rangeGR = np.asarray(std_range)

		return self.confluency_rangeGR, self.std_rangeGR

	def get_confluency_rangegeneticchemicalGR(self, compound, grup, posi, posj, posk):

		data_range = []  # this value is for DMSO control
		std_range = []

		std = compound.stdinh
		controlinformation = compound.control_info[-1]
		data = compound.grresultscombination

		controlinfo = compound.control_info

		grup = grup
		# for k in range(len(self.time_position)):

		for k in range(len(self.time_selectedgr)):
			posgr = self.time_positiongr[k]
			tmp1 = [data[posi][posj][posk][0][0][posgr]]  # this value is for DMSO control

			for j in range(len(grup.concentration)):
				tmp1.append(data[posi][posj][posk][j+1][0][posgr])
			data_range.append(tmp1)
		for k in range(len(self.time_selectedgr)):
			stdpo = self.time_position[k]
			tmp2 = [std[controlinformation][stdpo]]
			for j in range(len(grup.concentration)):
				pos = int(grup.concentration[j][-1])
				tmp2.append(std[pos][stdpo])
			std_range.append(tmp2)
		self.confluency_rangeGR = np.asarray(data_range)
		self.std_rangeGR = np.asarray(std_range)

		return self.confluency_rangeGR, self.std_rangeGR

	def get_confluency_rangecombination(self, compound, grup):
		concentration_range_zip = compound.concentration_rangeforzip
		data_range = []  # this value is for DMSO control
		std_range = []
		std = compound.stdnormalized_perc
		time = list(range(len(self.alltimepoints)))
		controlinformation = compound.control_info[-1]
		data = compound.inhib_data
		controlinfo = compound.control_info
		grup = grup
		for k in range(len(time)):
			tmp1 = [data[controlinformation][time[k]]]  # this value is for DMSO control
			tmp2 = [std[controlinformation][time[k]]]
			for j in range(len(grup.concentration)):
				pos = int(grup.concentration[j][-1])
				tmp1.append(data[pos][time[k]])
				tmp2.append(std[pos][time[k]])
			data_range.append(tmp1)
			std_range.append(tmp2)

		self.confluency_range = np.asarray(data_range)
		self.std_range = np.asarray(std_range)

		return self.confluency_range, self.std_range

	def get_concentration_range(self, grup):
		concentration_range = [0.001]  # this value is for DMSO control

		for j in range(len(grup.concentration)):
			concentration_range.append(grup.concentration[j][0])

		self.concentration_range = concentration_range

		return self.concentration_range

	@staticmethod
	def four_parameter_logistic(x, x0, k, l=0, b=1):

		# for LOEWE
		# b = Emax
		# l = Emin
		# x0 = m
		# k = lambda

		func = b + ((l - b) / (1 + ((x / x0) ** k)))

		return func

	def set_y_axis_name(self, inf):
		if inf == 0:
			y_axis_title = self.readout
		if inf == 1:
			y_axis_title = "relative " + self.readout
		if inf == 2:
			y_axis_title = "relative " + self.readout + " inhbition"

		return y_axis_title

	@staticmethod
	def set_y_axis_limt(inf):
		if inf == 0:
			y_axis_min = 0
			y_axis_max = 120

		if inf == 1:
			y_axis_min = -1
			y_axis_max = 1.5

		return y_axis_min, y_axis_max



