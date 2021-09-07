import numpy as np
import os
from HTSplotter.grupping import Groupping
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm
from matplotlib import ticker
from matplotlib.backends.backend_pdf import PdfPages
import math
import scipy.interpolate as inter
from mpl_toolkits.mplot3d import Axes3D
from scipy.optimize import curve_fit



class DoseResponse:

	def __init__(self, data, std, grup, row, concentration, confluency, timeposition, timepoint, txt_path, pdf_pages,
				 readout, readout_units):
		self.data_range = data
		self.std_range = std
		self.grup = grup
		self.row = row
		self.colun_num = 3
		self.concentration_range = concentration
		self.confluency_range = confluency
		self.time_position = timeposition
		self.time_selected = timepoint
		self.txt_path = txt_path
		self.pdf_pages = pdf_pages
		self.readout = readout
		self.readout_units = readout_units

		# Atributes for on and off situation

		self.compounds_list_on = None
		self.compounds_list_off = None
		self.compounds_name_on_off = None
		self.compounds_grups = []
		self.dose_curve_on_off = None
		self.data_curve_on_off = None
		self.std_curve_on_off = None
		self.concentration_range_on_off = None
		self.ec_concentration_on_off = None
		self.ic_concentration_on_off = None
		self.concentration_value_on_off = None

	# Methods for plotting requested on combination and singlecompound

	def plotting_1row_fourparameter(self, grup, info):
		cont = 0
		fitcurve_for_plot_all = []
		x_labe_fit = []
		compound_info = "compound: " + grup.name + "\n" + " cell line: " + grup.cell_name
		if len(self.time_position) == 1:
			rs, cols = self.row, 2
			fig, axs = plt.subplots(int(rs), int(2), figsize=(15, 5), sharey='row')
			fig.suptitle(compound_info, fontsize=10, fontweight='bold', x=0.3)
		else:
			rs, cols = self.row, 3
			position_titel = "center"
			fig, axs = plt.subplots(int(rs), int(cols), figsize=(15, 5))  # sharey='row',
			fig.suptitle(compound_info, fontsize=10, fontweight='bold',
						 ha='center')

		# this situation should be applied
		# if the number of plot require more than 1row have 1 row

		color = iter(cm.tab20b(np.linspace(0, 1, len(grup.concentration) + 1)))
		time_selected = ''
		data_for_plot = ''
		x_label_log = ''
		std_for_plot = ''
		for col in range(cols):
			if cont >= len(self.time_position):
				axs[col].axis('off')
				pass
			else:
				for k in range(len(self.confluency_range)):
					data_for_plot = self.confluency_range[cont]
					std_for_plot = self.std_range[cont]
					x_label_log = np.log10(self.concentration_range)
					time_selected = str(int(self.time_selected[cont]))

				title_for_each_plot = "Time: " + time_selected + " h"

				popt, pcov = self.get_curve_fit(self.concentration_range, data_for_plot)  # four-parameter logistic

				fit_curv_info, fit_curv_file = self.get_statistics(self.four_parameter_logistic,
																   self.concentration_range, data_for_plot,
																   std_for_plot, popt, pcov)

				area_under_curve = round(np.trapz(data_for_plot, self.concentration_range), 4)

				x_label_morepoints, tranf_x_label_morepoint = self.get_more_x_point(self.concentration_range)

				ic_selected, ic_concentration, ec_selected = self.get_ic(popt)

				title_x_label = "concentration range log$_{10}$ from: [" + str(self.concentration_range[1]) + \
								" ; " + str(self.concentration_range[-1]) + "] " + grup.unidade
				axs[col].set_title(title_for_each_plot, fontsize=8)
				axs[col].set_xlabel(title_x_label, fontsize=9)
				axs[col].set_ylabel(self.readout_units + ' '+self.readout)
				axs[col].set_ylim(0, 120)
				axs[col].yaxis.set_major_locator(ticker.FixedLocator(np.arange(0, 120, 20)))
				axs[col].xaxis.set_major_locator(ticker.FixedLocator(x_label_log))
				axs[col].xaxis.set_major_formatter(ticker.FixedFormatter(title_x_label))
				axs[col].set_xlim(self.concentration_range[1] * 0.6, self.concentration_range[-1] * 1.6)
				c = next(color)
				axs[col].errorbar(self.concentration_range[1:], data_for_plot[1:], yerr=std_for_plot[1:],
								  fmt='p', color=c)

				ic15_legeng, ic50_legeng, ec50_legend, ecmax_legend = self.get_ic_legend(grup, ic_concentration,
																						 ic_selected,
																						 ic_concentration[2],
																						 ic_concentration[9], popt[1],
																						 min(data_for_plot))

				self.get_icplotrow(axs, col, ic_concentration[2], ec_selected[2], 'indianred')
				if popt[1] > np.min(grup.concentration[:, 0]):
					if popt[1] < np.max(grup.concentration[:, 0]):
						axs[col].axvline(x=popt[1], linestyle="--", color="silver")
				self.get_icplotrow(axs, col, ic_concentration[9], ic_selected[9], 'black')
				axs[col].axhline(y=min(data_for_plot), linestyle="--", color="white")
				axs[col].semilogx(x_label_morepoints, self.four_parameter_logistic(x_label_morepoints,
																				   *popt), color=c)  # fit curve

				x_labe_fit.append(x_label_log)
				axs[col].legend((ic15_legeng, ec50_legend, ic50_legeng, ecmax_legend,
								 ("LL4=> " + fit_curv_info), "experiment_data"),
								loc="lower left", fontsize=7, frameon=False)

				fig.tight_layout()
				fig.subplots_adjust(top=0.88)
				self.save_ic_info(grup, cont, data_for_plot, ec50_legend, ic_selected, ic_concentration,
								  x_label_log, std_for_plot, area_under_curve,
								  fit_curv_info, popt)
				if info == 1:
					fitcurve_for_plot_all = self.four_parameter_logistic(x_label_morepoints, *popt)
					self.get_curve_information(grup, x_labe_fit, fitcurve_for_plot_all, x_label_morepoints,
											   data_for_plot, ic_selected, ec_selected, ic_concentration,
											   std_for_plot)

				fig.tight_layout()
				fig.subplots_adjust(top=0.88)
				cont += 1

		if len(self.time_position) % self.row:
			axs[-1, -1].axis('off')
		self.pdf_pages.savefig(fig)
		plt.close()
		plt.close("all")

	def plotting_rows_fourparameter(self, grup, info):
		cont = 0
		x_labe_fit = []
		rs, cols = self.row, 3
		compound_info = "compound: " + grup.name + "\n" + " cell line: " + grup.cell_name
		# this situation should be applied
		# if the number of plot require more than 1row have 1 row

		fig, axs = plt.subplots(int(rs), int(cols), figsize=(35, 15))  # sharey='row',
		fig.suptitle(compound_info, fontsize=18, fontweight='bold',
					 ha="center")
		color = iter(cm.tab20b(np.linspace(0, 1, len(grup.concentration)*10)))
		time_selected = ''
		data_for_plot = ''
		x_label_log = ''
		std_for_plot = ''
		for r in range(rs):

			for col in range(cols):
				if cont >= len(self.time_position):
					axs[r, col].axis('off')
					pass
				else:
					for k in range(len(self.confluency_range)):
						data_for_plot = self.confluency_range[cont]
						std_for_plot = self.std_range[cont]
						x_label_log = np.log10(self.concentration_range)
						time_selected = str(int(self.time_selected[cont]))

					title_for_each_plot = "Time: " + time_selected + " h"

					popt, pcov = self.get_curve_fit(self.concentration_range, data_for_plot)  # four-parameter logistic

					fit_curv_info, fit_curv_file = self.get_statistics(self.four_parameter_logistic,
																	   self.concentration_range, data_for_plot,
																	   std_for_plot, popt, pcov)

					area_under_curve = round(np.trapz(data_for_plot, self.concentration_range), 4)

					x_label_morepoints, tranf_x_label_morepoint = self.get_more_x_point(self.concentration_range)

					ic_selected, ic_concentration, ec_selected = self.get_ic(popt)

					title_x_label = "concentration range log$_{10}$ from: [" + str(self.concentration_range[1]) + \
									" ; " + str(self.concentration_range[-1]) + "] " + grup.unidade
					axs[r, col].set_title(title_for_each_plot, fontsize=16)
					axs[r, col].set_xlabel(title_x_label, fontsize=16)
					axs[r, col].set_ylabel(self.readout_units + ' ' + self.readout)
					axs[r, col].set_ylim(0, 120)
					axs[r, col].yaxis.set_major_locator(ticker.FixedLocator(np.arange(0, 120, 20)))
					# axs[col].yaxis.set_major_formatter(ticker.FixedFormatter("% confluence"))
					axs[r, col].xaxis.set_major_locator(ticker.FixedLocator(x_label_log))
					axs[r, col].xaxis.set_major_formatter(ticker.FixedFormatter(title_x_label))
					# axs[r, col].set_xticklabels(title_x_label, fontsize=16)
					axs[r, col].set_yticklabels(np.arange(0, 120, 20), fontsize=16)
					axs[r, col].set_xlim(self.concentration_range[1] * 0.6, self.concentration_range[-1] * 1.6)
					c = next(color)
					axs[r, col].errorbar(self.concentration_range[1:], data_for_plot[1:], yerr=std_for_plot[1:],
										 fmt='p', color=c)

					ic15_legeng, ic50_legeng, ec50_legend, ecmax_legend = self.get_ic_legend(grup, ic_concentration,
																							 ic_selected,
																							 ic_concentration[2],
																							 ic_concentration[9],
																							 popt[1],
																							 min(data_for_plot))

					self.get_icplotrows(axs, r, col, ic_concentration[2], ec_selected[2], 'indianred')
					if popt[1] > np.min(grup.concentration[:, 0]):
						if popt[1] < np.max(grup.concentration[:, 0]):
							axs[r, col].axvline(x=popt[1], linestyle="--", color="silver")
					self.get_icplotrows(axs, r, col, ic_concentration[9], ic_selected[9], 'black')
					axs[r, col].axhline(y=min(data_for_plot), linestyle="--", color="white")
					axs[r, col].semilogx(x_label_morepoints, self.four_parameter_logistic(x_label_morepoints,
																						  *popt), color=c)  # fit curve

					x_labe_fit.append(x_label_log)
					axs[r, col].legend((ic15_legeng, ec50_legend, ic50_legeng, ecmax_legend,
										("LL4=> " + fit_curv_info), "experiment_data"),
									   loc="lower left", fontsize=14, frameon=False)

					fig.tight_layout()
					fig.subplots_adjust(top=0.88)
					self.save_ic_info(grup, cont, data_for_plot, ec50_legend, ic_selected, ic_concentration,
									  x_label_log, std_for_plot, area_under_curve,
									  fit_curv_info, popt)
					if info == 1:
						fitcurve_for_plot_all = self.four_parameter_logistic(x_label_morepoints, *popt)
						self.get_curve_information(grup, x_labe_fit, fitcurve_for_plot_all,x_label_morepoints,
												   data_for_plot, ic_selected, ec_selected, ic_concentration,
												   std_for_plot)

					cont += 1

		if len(self.time_position) % self.row:
			axs[-1, -1].axis('off')
		self.pdf_pages.savefig(fig)
		plt.close()
		plt.close("all")


	@staticmethod
	def get_statistics(logfunc, x_label_log, data_for_plot, std_for_plot, popt, pcov):
		# gets standard deviations of the parameters (square root of the diagonal of the covariance
		stdves = np.sqrt(np.diag(pcov))
		# determine Chi-square

		chisq = sum(((data_for_plot - logfunc(x_label_log, *popt)) / std_for_plot) ** 2)

		dof = len(x_label_log) - len(popt)
		reduc_chisquared = chisq / dof

		# Residual sum of squares
		residuals = data_for_plot - logfunc(x_label_log, *popt)
		ss_res = np.sum(residuals ** 2)
		# total sum of squares
		ss_tot = sum((data_for_plot - np.mean(data_for_plot)) ** 2)
		ss_tot_sqrt = np.sqrt(ss_tot)

		# stadr error
		standar_err_mode = np.sqrt(
			sum(((logfunc(x_label_log, *popt) - data_for_plot) ** 2)) / \
			len(data_for_plot))

		# r_squared - value
		r_squared = 1 - (ss_res / ss_tot)
		perr = np.sqrt(np.diag(pcov))

		fit_curv_info = "fit_curve, $R^2$ = " + str("{:.4f}".format(r_squared)) + \
						" " + "$X^2$ = " + str("{:.1f}".format(chisq)) + \
						" Sum of squares = " + str("{:.2f}".format(ss_res))

		fit_curv_file = []

		return fit_curv_info, fit_curv_file

	def get_curve_fit(self, x_label_log, data_for_plot):
		uma = max(data_for_plot) - min(data_for_plot)

		p1 = [max(data_for_plot), np.median(x_label_log), 0.05, min(data_for_plot)]
		try:
			popt, pcov = curve_fit(self.four_parameter_logistic, x_label_log, data_for_plot, p1, method='lm')

		except RuntimeError as error:
			popt, pcov = curve_fit(self.four_parameter_logistic, x_label_log, data_for_plot, p1,
								   method='lm', maxfev=10000000)
		return popt, pcov

	def get_ic(self, popt):

		ic_selected = np.asarray([5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95])
		effect_inhibition = 100 - ic_selected
		ic_concentration = []
		ec_selected = 100 - ic_selected  # transform IC to EC--> effective concentration

		for k in range(len(effect_inhibition)):
			ic_concentration.append(popt[1] *
									((((popt[0] - popt[-1]) / (effect_inhibition[k] - popt[-1])) - 1)
									 ** (1 / popt[2])))

		return ic_selected, ic_concentration, ec_selected

	def save_ic_info(self, grup, cont, data_for_plot, ec_selected, ic_selected, ic_concentration,
					 x_label_log, std_for_plot, area_under_curve, fit_curv_info,
					 popt):

		f = open(self.txt_path, 'a')
		# ##get function

		popinfo = "pop2 = > input max(data_for_plot), np.median(10 ** x_label_log), 0.05, min(data_for_plot)" + \
				  '\t'
		fit_curv_info1 = fit_curv_info.split(" ")
		statistctitle = 'statistcs- four parameters' + fit_curv_info1[0]+ '\t'
		fit_curv_info1.pop(0)
		fit_curv_info1[6] = " ".join(fit_curv_info1[6:9])
		fit_curv_info1.pop(7)
		fit_curv_info1.pop(7)
		com_max = self.set_headers_writing(f, grup, cont, statistctitle, popinfo, data_for_plot,
										   ic_selected, fit_curv_info1)

		self.set_writing_file(f, com_max, x_label_log, data_for_plot, std_for_plot, ec_selected,
							  ic_selected, ic_concentration, area_under_curve, fit_curv_info1, popt)

		f.close()

	def set_headers_writing(self, f, grup, cont, statistctitle, popinfo, data_for_plot,
							ic_selected,fit_curv_info1):
		f.write('\n')
		f.write("Cell line: " + grup.cell_name + '\n')
		f.write("Compound name: " + grup.name + '\n')
		f.write('time_point' + '\t' + str(self.time_selected[cont]) + '\n')
		f.write('Concentration' + '\t')
		f.write('Concentration log scale' + '\t')
		f.write(self.readout + ' ' + self.readout_units + '\t')
		f.write('Std' + '\t')
		f.write('IC' + '\t')
		f.write('IC_concentration' + '\t')
		f.write(statistctitle )
		f.write('EC 50' + '\t')
		f.write('AUC' + '\t')
		f.write(popinfo)

		f.write('\n')
		com_max = max(len(self.concentration_range), len(data_for_plot), len(ic_selected))
		return com_max

	def set_writing_file(self, f, com_max, x_label_log, data_for_plot, std_for_plot, ec50_value,
						 ic_selected, ic_concentration, area_under_curve, fit_curv_info, popt):
		cont = 0
		cont1 = 0
		c = 0
		c2 = 3
		for i in range(com_max):
			if cont < len(self.concentration_range):
				f.write(str(self.concentration_range[cont]) + '\t')
				f.write(str(x_label_log[cont]) + '\t')
				f.write(str(data_for_plot[cont]) + '\t')
				f.write(str(std_for_plot[cont]) + '\t')
			if cont >= len(self.concentration_range):
				f.write('\t')
				f.write('\t')
				f.write('\t')
				f.write('\t')
			if cont1 < len(ic_selected):
				f.write(str(ic_selected[cont1]) + '\t')
				f.write(str(ic_concentration[cont1])+ '\t')
			if cont1 <= 3:
				f.write(" ".join(fit_curv_info[c:c2]) + '\t')
				c = c2
				c2= c2+3
			if i == 0:
				f.write(str(ec50_value) + '\t')
				f.write(str(area_under_curve) + '\t')
				f.write(" ".join(popt.astype(str)) + '\t')

			f.write('\n')
			cont += 1
			cont1 += 1

		return

	def get_ic_legend(self, grup, ic_concentration, ic_selected, concentration_value15,
					  concentration_value, ec50_value, ecmax):
		if math.isnan(ic_concentration[2]) or ic_concentration[2] < np.min(grup.concentration[:, 0]) \
				or ic_concentration[2] > np.max(grup.concentration[:, 0]):

			ic15_legeng = "IC_" + str(ic_selected[2]) + " Value outside the tested range!"
		else:
			ic15_legeng = "IC_" + str(ic_selected[2]) + " = " + str("{:.2f}".format(concentration_value15)) \
						  + ' ' + grup.unidade

		if math.isnan(ic_concentration[9]) or ic_concentration[9] < np.min(grup.concentration[:, 0]) \
				or ic_concentration[9] > np.max(grup.concentration[:, 0]): # self.concentration_range[-1] + 2:

			ic50_legeng = "IC_" + str(ic_selected[9]) + " Value outside the tested range!"
		else:
			ic50_legeng = "IC_" + str(ic_selected[9]) + " = " + str("{:.2f}".format(concentration_value)) \
						  + ' ' + grup.unidade
		if ec50_value < np.min(grup.concentration[:, 0])\
				or ec50_value > np.max(grup.concentration[:, 0]):
			ec50_legend = "EC_50" + "Not possible to calculate"
		else:
			ec50_legend = "EC_50 = " + str("{:.2f}".format(ec50_value)) + ' ' + grup.unidade
		ecmax_legend = "Emax = " + str("{:.2f}".format(ecmax)) + " % confluency"

		return ic15_legeng, ic50_legeng, ec50_legend, ecmax_legend

	def get_ic_on_off_legend(self, ic_on, ic_off, grup, compound_on, compound_off,
							 concentration_on, concentration_off):

		if math.isnan(concentration_on) or concentration_on < np.min(grup.concentration[:, 0]) \
				or concentration_on > np.max(grup.concentration[:, 0]):
			ic_on_legeng = "IC_" + str(ic_on) + " Value outside the tested range!"
		else:
			ic_on_legeng = "IC_" + str(ic_on) + " = " + str("{:.2f}".format(float(concentration_on))) + \
						   ' ' + grup.unidade
		if math.isnan(concentration_off) or concentration_off < np.min(grup.concentration[:, 0]) \
				or concentration_off > np.max(grup.concentration[:, 0]):

			ic_off_legeng = "IC_" + str(ic_off) + " Value outside the tested range!"
		else:
			ic_off_legeng = "IC_" + str(ic_off) + " = " + str("{:.2f}".format(float(concentration_off))) + \
							' ' + grup.unidade

		compound_on_legend_fit = "fit_curve_" + compound_on
		compound_off_legend_fit = "fit_curve_" + compound_off
		compound_on_legend_exper = "Expreiment_data_" + compound_on
		compound_off_legend_exper = "Expreiment_data_" + compound_off

		return ic_on_legeng, ic_off_legeng, compound_on_legend_fit,\
			compound_off_legend_fit, compound_on_legend_exper, compound_off_legend_exper

	def get_ic15plotrows(self, axs, r, col, ic_concentration, ic_selected):

		if np.any(ic_concentration[2] <= self.concentration_range[1:]):
			axs[r, col].plot(ic_concentration[2], ic_selected[2], "o", color="indianred")
		else:  # case IC15 is was not measured
			axs[r, col].plot(2 * ic_concentration[-1], ic_selected[2], "o", color="indianred")

	def get_ic15plotrow(self, axs, col, ic_concentration, ic_selected):

		if np.any(ic_concentration[2] <= self.concentration_range[1:]):
			axs[col].plot(ic_concentration[2], ic_selected[2], "o", color="dimgray")
		else:  # case IC15 is was not measured
			axs[col].plot(2 * ic_concentration[-1], ic_selected[2], "o", color="dimgray")

	def get_icplotrows(self, axs, r, col, ic_concentration, ic_selected, cor):
		if np.any(ic_concentration <= self.concentration_range[1:]):
			axs[r, col].plot(ic_concentration, ic_selected, "o", color=cor)
		else:  # case IC50 is was not measured
			axs[r, col].plot(2 * ic_concentration, ic_selected, "o", color=cor)

	def get_icplotrow(self, axs, col, ic_concentration, ic_selected, cor):
		if np.any(ic_concentration <= self.concentration_range[1:]):
			axs[col].plot(ic_concentration, ic_selected, "o", color=cor)
		else:  # case IC50 is was not measured
			axs[col].plot(2 * ic_concentration, ic_selected, "o", color=cor)

	def get_curve_information(self, grup, x_labe_fit, fitcurve_for_plot_all, x_label_morepoints,
							  data_for_plot, ic_selected, ec_selected, ic_concentration, std_for_plot):
		self.compounds_name_on_off.append(grup.name)
		self.concentration_range_on_off = x_label_morepoints
		self.dose_curve_on_off.append(fitcurve_for_plot_all)
		self.data_curve_on_off.append(data_for_plot)
		self.std_curve_on_off.append(std_for_plot)
		self.ec_concentration_on_off = ec_selected
		self.ic_concentration_on_off.append(ic_selected)
		self.concentration_value_on_off.append(ic_concentration)

	@staticmethod
	def get_more_x_point(x_label_log):
		x_label_morepoints = (np.arange(np.min(x_label_log), np.max(x_label_log) * 1.2,
										(np.max(x_label_log) * 1.2 - np.min(x_label_log)) / 10000))
		tranf_x_label_morepoint = 10 ** x_label_morepoints

		return x_label_morepoints, tranf_x_label_morepoint

	@staticmethod
	def four_parameter_logistic(x, l, x0, k, b):

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


class DoseResponseSingle(DoseResponse):
	def __init__(self, data, std, grup, row, concentration, confluency, timeposition, timepoint, txt_path, pdf_pages,
				 readout, readout_units):
		super().__init__(data, std, grup, row, concentration, confluency, timeposition, timepoint, txt_path, pdf_pages,
						 readout, readout_units)

		if self.row == 1:
			self.plotting_1row_fourparameter(grup, 0)
		else:
			self.plotting_rows_fourparameter(grup, 0)


class DoseResponseCondition(DoseResponse):
	def __init__(self, data, std, grup, row, concentration, confluency, timeposition, timepoint, txt_path, pdf_pages,
				 readout, readout_units):
		super().__init__(data, std, grup, row, concentration, confluency, timeposition, timepoint, txt_path, pdf_pages,
						 readout, readout_units)

	def singlecompoundplot(self, grup, compounds_name_on_off,  dose_curve_on_off,
						   data_curve_on_off, std_curve_on_off, ic_concentration_on_off,
						   concentration_value_on_off):
		self.compounds_name_on_off = compounds_name_on_off
		self.dose_curve_on_off = dose_curve_on_off
		self.data_curve_on_off = data_curve_on_off
		self.std_curve_on_off = std_curve_on_off
		self.ic_concentration_on_off = ic_concentration_on_off
		self.concentration_value_on_off = concentration_value_on_off

		if self.row == 1:
			self.plotting_1row_fourparameter(grup, 1)
		else:
			self.plotting_rows_fourparameter(grup, 1)

# method only for on and off situation

	def sigmoid_on_off_summary(self, grup, compounds_name_on_off,  dose_curve_on_off,
							   data_curve_on_off, std_curve_on_off, ic_concentration_on_off,
							   concentration_value_on_off, concentration_range):
		if self.row == 1:
			self.plotting_1row_onoff_parametric(grup)
		else:
			self.plotting_rows_onoff_parametric(grup)

	def plotting_1row_onoff_parametric(self, grup):
		cont = 0
		step = int(len(self.compounds_name_on_off) / 2)
		rs, cols = self.row, 3
		compound_info = "compound: " + self.compounds_name_on_off[0] + " and " + self.compounds_name_on_off[
			step] + "\n" + " cell line: " + grup.cell_name
		# this situation should be applied
		# if the number of plot require more than 1row have 1 row
		# axs = fig.add_subplot(row_numeber_per_figure, cols, cont+1)
		if len(self.time_position) == 1:
			rs, cols = self.row, 2
			fig, axs = plt.subplots(int(rs), int(2), figsize=(15, 5), sharey='row')
			fig.suptitle(compound_info, fontsize=10, fontweight='bold', x=0.3)
		else:
			rs, cols = self.row, 3
			fig, axs = plt.subplots(int(rs), int(cols), figsize=(15, 5))  # sharey='row',
			fig.suptitle(compound_info, fontsize=10, fontweight='bold',
						 ha='center')

		color = iter(cm.tab20b(np.linspace(0, 1, len(self.compounds_name_on_off) + 1)))
		ste = step
		con = 0
		for col in range(cols):
			if cont >= len(self.time_position):
				axs[col].axis('off')
				pass
			else:
				compound_on = self.compounds_name_on_off[con]
				compound_off = self.compounds_name_on_off[ste]
				data_for_plot_on = self.data_curve_on_off[con]
				data_for_plot_off = self.data_curve_on_off[ste]
				std_for_plot_on = self.std_curve_on_off[con]
				std_for_plot_off = self.std_curve_on_off[ste]
				curve_fit_on = self.dose_curve_on_off[con]
				curve_fit_off = self.dose_curve_on_off[ste]
				ic_on = self.ic_concentration_on_off[con][9]
				ic_off = self.ic_concentration_on_off[ste][9]

				concentration_on = self.concentration_value_on_off[con][9]
				concentration_off = self.concentration_value_on_off[ste][9]

				time_selected = str(int(self.time_selected[cont]))
				ste += 1
				con += 1
				title_for_each_plot = "Time: " + time_selected + " h"

				title_x_label = "concentration range log$_{10}$ from: [" + str(self.concentration_range[1]) + \
								" ; " + str(self.concentration_range[-1]) + "] " + grup.unidade
				axs[col].set_title(title_for_each_plot, fontsize=8)
				axs[col].set_xlabel(title_x_label, fontsize=9)
				axs[col].set_ylabel(self.readout_units + ' ' + self.readout)
				axs[col].set_ylim(0, 120)

				if (self.concentration_range[-1] * 1.6) / (self.concentration_range[1] * 0.6) >= 10.0:
					axs[col].set_xlim(self.concentration_range[1] * 0.6, self.concentration_range[-1] * 1.6)
				elif (self.concentration_range[-1] * 1.6) / (self.concentration_range[1] * 0.3) >= 10.0:

					axs[col].set_xlim(self.concentration_range[1] * 0.3, self.concentration_range[-1] * 1.6)
				else:
					axs[col].set_xlim(self.concentration_range[1] * 0.1, self.concentration_range[-1] * 1.6)

				c = next(color)

				axs[col].errorbar(self.concentration_range[1:], data_for_plot_on[1:], yerr=std_for_plot_on[1:], fmt='p',
								  color=c)

				axs[col].semilogx(self.concentration_range_on_off, curve_fit_on,
								  color=c)  # fit curve
				axs[col].plot(10 ** np.asarray(ic_on), self.ec_concentration_on_off[9], "o", color="dimgray")
				c = next(color)

				axs[col].errorbar(self.concentration_range[1:], data_for_plot_off[1:], yerr=std_for_plot_off[1:],
								  fmt='p',
								  color=c)
				axs[col].semilogx(self.concentration_range_on_off, curve_fit_off,
								  color=c)

				axs[col].plot(10 ** np.asarray(ic_off), self.ec_concentration_on_off[9], "o", color="black")

				ic_on_legend, ic_off_legend, compound_on_legend_fit, \
					compound_off_legend_fit, \
					compound_on_legend_exper, compound_off_legend_exper = self.get_ic_on_off_legend(ic_on, ic_off, grup,
																									compound_on,
																									compound_off,
																									concentration_on,
																									concentration_off)
				axs[col].legend((compound_on_legend_fit, ic_on_legend, compound_off_legend_fit, ic_off_legend,
								 compound_on_legend_exper, compound_off_legend_exper), loc="lower left", fontsize=7,
								frameon=False)
				fig.tight_layout()
				fig.subplots_adjust(top=0.88)

				cont += 1
		if len(self.time_position) % self.row:
			axs[-1, -1].axis('off')

		self.pdf_pages.savefig(fig)
		plt.close()
		plt.close("all")

	def plotting_rows_onoff_parametric(self, grup):
		cont = 0
		step = int(len(self.compounds_name_on_off) / 2)
		rs, cols = self.row, 3
		compound_info = "compound: " + self.compounds_name_on_off[0] + " and " + self.compounds_name_on_off[
			step] + "\n" + " cell line: " + grup.cell_name
		# this situation should be applied
		# if the number of plot require more than 1row have 1 row
		# axs = fig.add_subplot(row_numeber_per_figure, cols, cont+1)
		fig, axs = plt.subplots(rs, cols, figsize=(30, 15))  # sharey='row',
		fig.suptitle(compound_info, fontsize=18, fontweight='bold',
					 ha="center")
		color = iter(cm.tab20b(np.linspace(0, 1, len(self.compounds_name_on_off) + 1)))
		ste = step
		con = 0
		for r in range(rs):
			for col in range(cols):
				if cont >= len(self.time_position):
					axs[r, col].axis('off')
					pass
				else:
					compound_on = self.compounds_name_on_off[con]
					compound_off = self.compounds_name_on_off[ste]
					data_for_plot_on = self.data_curve_on_off[con]
					data_for_plot_off = self.data_curve_on_off[ste]
					std_for_plot_on = self.std_curve_on_off[con]
					std_for_plot_off = self.std_curve_on_off[ste]
					curve_fit_on = self.dose_curve_on_off[con]
					curve_fit_off = self.dose_curve_on_off[ste]
					ic_on = self.ic_concentration_on_off[con][9]
					ic_off = self.ic_concentration_on_off[ste][9]

					concentration_on = self.concentration_value_on_off[con]
					concentration_off = self.concentration_value_on_off[ste]

					time_selected = str(int(self.time_selected[cont]))
					ste += 1
					con += 1
					title_for_each_plot = "Time: " + time_selected + " h"

					title_x_label = "concentration range log$_{10}$ from: [" + str(self.concentration_range[1]) + \
									" ; " + str(self.concentration_range[-1]) + "] " + grup.unidade
					axs[r, col].set_title(title_for_each_plot, fontsize=16)
					axs[r, col].set_xlabel(title_x_label, fontsize=16)
					axs[r, col].set_ylabel(self.readout_units + ' ' + self.readout)
					axs[r, col].set_ylim(0, 120)

					if (self.concentration_range[-1] * 1.6) / (self.concentration_range[1] * 0.6) >= 10.0:
						axs[r, col].set_xlim(self.concentration_range[1] * 0.6, self.concentration_range[-1] * 1.6)
					elif (self.concentration_range[-1] * 1.6) / (self.concentration_range[1] * 0.3) >= 10.0:

						axs[r, col].set_xlim(self.concentration_range[1] * 0.3, self.concentration_range[-1] * 1.6)
					else:
						axs[r, col].set_xlim(self.concentration_range[1] * 0.1, self.concentration_range[-1] * 1.6)

					c = next(color)

					axs[r, col].errorbar(self.concentration_range[1:], data_for_plot_on[1:], yerr=std_for_plot_on[1:],
										 fmt='p', color=c)

					axs[r, col].semilogx(self.concentration_range_on_off, curve_fit_on, color=c)  # fit curve
					if concentration_on[9] != "nan":
						axs[r, col].plot(concentration_on[9], ic_on,"o", color="dimgray")

					c = next(color)

					axs[r, col].errorbar(self.concentration_range[1:], data_for_plot_off[1:],
										 yerr=std_for_plot_off[1:], fmt='p', color=c)
					axs[r, col].semilogx(self.concentration_range_on_off, curve_fit_off, color=c)
					if concentration_off[9] != "nan":
						axs[r, col].plot(concentration_off[9], ic_off, "o", color="black")

					ic_on_legend, ic_off_legend, compound_on_legend_fit, \
					compound_off_legend_fit, \
					compound_on_legend_exper, compound_off_legend_exper = self.get_ic_on_off_legend(ic_on, ic_off, grup,
																									compound_on,
																									compound_off,
																									concentration_on[9],
																									concentration_off[9])

					axs[r, col].legend((compound_on_legend_fit, ic_on_legend, compound_off_legend_fit, ic_off_legend,
										compound_on_legend_exper, compound_off_legend_exper), loc="lower left",
									   fontsize=14, frameon=False)
					fig.tight_layout()
					fig.subplots_adjust(top=0.88)

					cont += 1

		if len(self.time_position) % self.row:
			axs[-1, -1].axis('off')

		self.pdf_pages.savefig(fig)
		plt.close()
		plt.close("all")


