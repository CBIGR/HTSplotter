import numpy as np

import os
from grupping import Groupping
from grupping import Groupping, Data_group, Grouppingmedium
import matplotlib.pyplot as plt
from doseresponse import DoseResponse
from plotting import Overtime
from timepointselection import Timepointselection
from commonfunctions import CommonFunctions, OpenPdf, SaveTXTinfo, GrwothRateSaveTXT
from doseresponseplots import DoseResponsePlot, DoseResponseSinglePlotting, DoseResponseGrowthRate
from growthrate import GrowthRate, GrowthRateCompoundscreenSeveral
from matplotlib.pyplot import cm
from matplotlib.backends.backend_pdf import PdfPages
import math


class ExperimentSingleCompound:
    def __init__(self, branch, file_info, information_readout, readout_units, biologicalreplicate,
                 hdfcompound, file_names, i):

        # ##Usefule information for each experiment, information gave by the calls==>new_save_data_base:

        self.branch = branch
        self.readout = information_readout
        self.readout_units = readout_units
        self.elapse = file_info.elapsed
        self.date_info = file_info.date_info
        self.biologicalreplicateinfo = biologicalreplicate

        #  from hdfcompound object in main
        self.seeding = hdfcompound.seeding
        self.celine = hdfcompound.celline
        self.condition = hdfcompound.condition
        self.control = hdfcompound.control
        self.header = hdfcompound.fields
        self.stdinfo = hdfcompound.std_info

        self.originaldata = hdfcompound.data
        self.inhib_data = hdfcompound.inhibited
        self.stdinh = hdfcompound.std_inh

        self.mediumfiels = hdfcompound.fieldsmedium
        self.mediuminhfiels = hdfcompound.fieldsmediuminhibited
        self.mediumdata = hdfcompound.datamedium
        self.mediumstd = hdfcompound.stdmedium
        self.inhibitemediumdata = hdfcompound.inhibitedmedium
        self.translationmediumdata = hdfcompound.normalizedtranslationmedium
        self.inhibitemediumstd = hdfcompound.std_inhmedium

        self.compound_alone = hdfcompound.compoundalone

        # for Dose-response curve
        self.normalized_perc = hdfcompound.normalized_perc
        self.inhibitedperc = hdfcompound.inhibitedperc
        self.stdnormalized_perc = hdfcompound.std
        self.normtozero = hdfcompound.normtozero

        #  from hdfcompound object in main
        self.fileioriginaldatapath = file_names.fileioriginaldatapath
        self.fileinhibiteddatapath = file_names.fileinhibiteddatapath
        self.filenormalizedatapath = file_names.filenormalizedatapath
        self.filegrowthrateresults = file_names.filegrowthrateresults
        self.pdf_path = file_names.filepdfresultspath
        self.txt_path = file_names.fileictxtresultspath

        # #####
        self.experiment_name = i

        self.pdf_pages = None
        self.plot_titel = None
        self.time_position = None
        self.time_selected = None
        self.comalone_list_group = None
        self.data_plot = None
        self.std_plot = None
        self.plotting = None
        self.concentration_range = None
        self.doseresponsecurve = None
        self.startpdf = None
        self.control_info = None

        # growth rate
        self.grresults = None
        self.time_selectedgr = None
        self.grheader = None
        self.position_gr = []
        self.position_elapsegr = []
        self.grselectedelap = []
        self.grselectedtimepoint = []
        self.icgrowthratetxt_path = file_names.fileicgrtxtresultspath
        if len(self.elapse) >1:
            self.extremepoints = [self.elapse[1], self.elapse[-2]]
            self.extremepointspositions = []

        # get main time point
        timeselect = Timepointselection(self.elapse)
        self.time_selected = timeselect.time_selected
        self.time_position = timeselect.time_position

        # get from doseresponse script
        self.curves = None
        self.curvesgr = None

        self.make_empty_list_single()
        self.get_headers_single()

        self.startcommonfunctions = CommonFunctions(self)

        self.row_num = self.startcommonfunctions.row_num
        self.colun_num = self.startcommonfunctions.colun_num


    def get_control(self):
        for j in self.header:
            for i in range(len(self.control)):
                for k in range(len(self.control[i])):
                    # per cell line
                    if self.control[i][k][0] == j[0] and self.control[i][k][1] == j[1] and \
                            self.control[i][k][2] == j[2] and self.control[i][k][3] == j[3]:
                        self.control[i][k] = j

    def make_empty_list_single(self):
        list1 = []
        for j in range(len(self.compound_alone)):
            header_count = len(self.compound_alone[j])
            headers = [[] for i in range(0, header_count)]  # list_possible_combination
            list1.append(headers)
        self.comalone_list_group = list1

    def get_headers_single(self):

        for i in range(len(self.compound_alone)):

            for k in range(len(self.compound_alone[i])):
                list1 = []
                for u in self.header:
                    if u[0] == self.celine[i] and u[3] == self.compound_alone[i][k]:
                        list1.append(u)
                self.comalone_list_group[i][k].append(list1)

    def close_pdf(self):
        self.startpdf.close_pdf()


class SingleCompound(ExperimentSingleCompound):
    def __init__(self,  branch, elapsed, information_readout, readout_units, biologicalreplicate,
                 hdfcompound, file_names, i):
        super().__init__(branch, elapsed, information_readout, readout_units, biologicalreplicate,
                         hdfcompound, file_names, i)

        self.get_control()
        if len(self.elapse) > 1:
            self.position_gr = []
            gr = GrowthRateCompoundscreenSeveral(self)
            self.grresults = gr.grresults
            for eac in self.time_selected:
                self.position_gr.append(np.searchsorted(gr.xs, eac))
            for eac in self.elapse[1:-1]:
                self.position_elapsegr.append(np.searchsorted(gr.xs, eac))
            self.position_elapsegr.append(np.searchsorted(gr.xs, self.elapse[-2]))
            for eac in self.extremepoints:
                self.extremepointspositions.append(np.searchsorted(gr.xs, eac))
            if self.position_gr[-1] == self.position_gr[-1]:
                self.position_gr[-1] = np.searchsorted(gr.xs, gr.xs[-2])
            self.time_selectedgr = self.time_selected
            self.grheader = gr.grheader
            GrwothRateSaveTXT(self)
        SaveTXTinfo(self, self.fileioriginaldatapath, self.originaldata)
        SaveTXTinfo(self, self.filenormalizedatapath, self.normalized_perc)
        SaveTXTinfo(self, self.fileinhibiteddatapath, self.inhib_data)
        # this subclass is for single compound screen with several controls
        self.startpdf = OpenPdf(self)
        self.pdf_path = self.startpdf.pdf_path
        self.pdf_pages = self.startpdf.pdf_pages
        self.plotting = Overtime(self)

    def doseresponse(self, growthrate=0):
        self.data_plot = self.normalized_perc
        self.std_plot = self.stdnormalized_perc
        self.doseresponsecurve = DoseResponse(self)
        for i in range(len(self.comalone_list_group)):
            for j in range(len(self.comalone_list_group[i])):
                self.control_info = self.control[i][j]
                grup = [Groupping(c, self.compound_alone[i][j], self.celine[i], self.seeding,
                                  self.condition) for c in self.comalone_list_group[i][j]]
                if len(grup[0].concentration) > 2:
                    self.concentration_range = self.doseresponsecurve.get_concentration_range(grup[0])
                    datanormperce, std = self.doseresponsecurve.get_confluency_range(self, grup[0])
                    self.doseresponsecurve.get_curve_specifictimepoints(self, datanormperce, grup[0])
                    self.curves = self.doseresponsecurve.curves
                    DoseResponseSinglePlotting(self, grup[0], datanormperce, std)

                    # Growth rate
                    if growthrate== 0:
                        datagr, stdgr = self.doseresponsecurve.get_confluency_rangeGR(self, grup[0], i, j)
                        self.doseresponsecurve.get_curve_specifictimepointsGR(self, datagr, grup[0])
                        self.curvesgr = self.doseresponsecurve.curves
                        DoseResponseGrowthRate(self, grup[0], i, j, datagr, stdgr)

    def inhibitionovertime(self, data, std, info1=0, info2=0):
        self.data_plot = data
        self.std_plot = std

        for i in range(len(self.comalone_list_group)):
            # cell line
            # # plot all compounds alone
            for j in range(len(self.comalone_list_group[i])):
                self.control_info = self.control[i][j]
                grup = [Groupping(c, self.compound_alone[i][j], self.celine[i], self.seeding[i],
                                  self.condition[i]) for c in self.comalone_list_group[i][j]]
            #
                self.plotting.plot_single_compond(grup[0], self, info1, info2)
                if info1 == 0:
                    self.plotting.plot_grwothrate(self, grup[0], i, j, 3, 2)
            if info1 != 0:
                try:
                    for j in range(len(self.mediuminhfiels[i])):
                        self.control_info = self.control[i][j]
                        grup = [Grouppingmedium(self.mediuminhfiels[i][j],  self.mediuminhfiels[i][j][3], self.celine[i],
                                                self.seeding[i], self.condition[i])]

                        self.plotting.plot_control_mediuminhibited(grup[0], self, i, j, info1, info2)

                except IndexError:
                    pass
        # self.data_plot = self.normtozero
        # if info1 !=0:
        #     for i in range(len(self.comalone_list_group)):
        #         # cell line
        #         # # plot all compounds alone
        #         for j in range(len(self.comalone_list_group[i])):
        #             self.control_info = self.control[i][j]
        #             grup = [Groupping(c, self.compound_alone[i][j], self.celine[i], self.seeding[i],
        #                               self.condition[i]) for c in self.comalone_list_group[i][j]]
        #             #
        #             self.plotting.plot_single_compond(grup[0], self, info1, info2)

class SingleCompoundonecontrol(ExperimentSingleCompound):
    def __init__(self,  branch, elapsed, information_readout, readout_units, biologicalreplicate,
                 hdfcompound, file_names, i):
        super().__init__(branch, elapsed, information_readout, readout_units, biologicalreplicate,
                         hdfcompound, file_names, i)

        if len(self.elapse) > 1:
            self.position_gr = []
            gr = GrowthRateCompoundscreenSeveral(self, 1)
            self.grresults = gr.grresults
            for eac in self.time_selected:
                self.position_gr.append(np.searchsorted(gr.xs, eac))
            for eac in self.elapse[1:-1]:
                self.position_elapsegr.append(np.searchsorted(gr.xs, eac))
            self.position_elapsegr.append(np.searchsorted(gr.xs, self.elapse[-2]))
            for eac in self.extremepoints:
                self.extremepointspositions.append(np.searchsorted(gr.xs, eac))
            if self.position_gr[-1] == self.position_gr[-1]:
                self.position_gr[-1] = np.searchsorted(gr.xs, gr.xs[-2])
            self.time_selectedgr = self.time_selected
            self.grheader = gr.grheader
            GrwothRateSaveTXT(self)
        SaveTXTinfo(self, self.fileioriginaldatapath, self.originaldata)
        SaveTXTinfo(self, self.filenormalizedatapath, self.normalized_perc)
        SaveTXTinfo(self, self.fileinhibiteddatapath, self.inhib_data)

        self.startpdf = OpenPdf(self)
        self.pdf_path = self.startpdf.pdf_path
        self.pdf_pages = self.startpdf.pdf_pages
        self.plotting = Overtime(self)

    def doseresponseonecontrol(self, growthrate=0):
        self.data_plot = self.normalized_perc
        self.std_plot = self.stdnormalized_perc

        self.doseresponsecurve = DoseResponse(self)

        for i in range(len(self.comalone_list_group)):
            for j in range(len(self.comalone_list_group[i])):

                self.control_info = self.control[i]

                grup = [Groupping(c, self.compound_alone[i][j], self.celine[i], self.seeding,
                                  self.condition) for c in self.comalone_list_group[i][j]]
                if len(grup[0].concentration) > 2:
                    self.concentration_range = self.doseresponsecurve.get_concentration_range(grup[0])
                    datanormperce, std = self.doseresponsecurve.get_confluency_range(self, grup[0])
                    self.doseresponsecurve.get_curve_specifictimepoints(self, datanormperce, grup[0])
                    self.curves = self.doseresponsecurve.curves
                    DoseResponseSinglePlotting(self, grup[0], datanormperce, std)
                    # self.get_concentration_range(grup)
                    # self.get_confluency_range(grup)
                    # DoseResponseSingle(self.data_plot, self.std_range, grup[0], self.row_num,
                    #                    self.concentration_range,
                    #                    self.confluency_range, self.time_position, self.time_selected,
                    #                    self.txt_path,
                    #                    self.pdf_pages, self.readout, self.readout_units)

                    # Growth rate
                    if growthrate == 0:
                        datagr, stdgr = self.doseresponsecurve.get_confluency_rangeGR(self, grup[0], i, j)
                        self.doseresponsecurve.get_curve_specifictimepointsGR(self, datagr, grup[0])
                        self.curvesgr = self.doseresponsecurve.curves
                        DoseResponseGrowthRate(self, grup[0], i, j, datagr, stdgr)


    def inhibitiononecontrolovertime(self, data, std, info1=0, info2=0):
        self.data_plot = data
        self.std_plot = std
        for i in range(len(self.comalone_list_group)):
            # cell line
            # # plot all compounds alone
            for j in range(len(self.comalone_list_group[i])):
                self.control_info = self.control[i]

                grup = [Groupping(c, self.compound_alone[i][j], self.celine[i], self.seeding[i],
                                  self.condition[i]) for c in self.comalone_list_group[i][j]]
            #
                self.plotting.plot_single_compond(grup[0], self, info1, info2)
                # Growth Rate
                if info1 == 0:
                    self.plotting.plot_grwothrate(self, grup[0], i, j, 3, 2)

