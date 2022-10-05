import numpy as np
import os
from grupping import Groupping, Data_group, Grouppingzip
from synergism import Blissmethod, Hsamethod, LoeweMethod, ZipMethod, ZipMethomultidimen
from plotting import Overtime
from timepointselection import Timepointselection
from itertools import combinations
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm
from matplotlib.backends.backend_pdf import PdfPages
import math
import scipy.interpolate as inter
from doseresponse import DoseResponse
from synergyprocessing import Commonprocessing
from commonfunctions import CommonFunctions, OpenPdf, SaveTXTinfo, AdjustData, SynergyScoreSaveTXT, GrwothRateSaveTXT
from doseresponseplots import DoseResponsePlot, DoseResponseSinglePlotting, DoseResponseGrowthRate
from growthrate import GrowthRate, GrowthRateCompoundscreenSeveral, GrowthRateCombination
import time

import time as ti
class ExperimentCombination:

    def __init__(self, synergymethod, hdfcompoundcomb, file_info, branch, information_readout, readout_units,
                 fileobject, biologicalreplicate, i):

        self.synergy_method = synergymethod  # 0-->Bliss; 1-->HSA; 2-->ZIP
        self.readout = information_readout
        self.readout_units = readout_units
        self.elapse = file_info.elapsed

        self.date_info = file_info.date_info
        self.branch = branch
        self.experiment_name = i
        self.biologicalreplicateinfo = biologicalreplicate

        self.fileioriginaldatapath = fileobject.fileioriginaldatapath
        self.fileinhibiteddatapath = fileobject.fileinhibiteddatapath
        self.filenormalizedatapath = fileobject.filenormalizedatapath
        self.filegrowthrateresults = fileobject.filegrowthrateresults

        #  from hdfcompound object in main
        self.stdinfo = hdfcompoundcomb.std_info
        self.medium = hdfcompoundcomb.medium
        self.control = hdfcompoundcomb.control
        self.celine = hdfcompoundcomb.celline
        self.condition = hdfcompoundcomb.condition
        self.seeding = hdfcompoundcomb.seeding
        self.compound_alone = hdfcompoundcomb.compoundalone
        self.possible_comb = hdfcompoundcomb.possiblecombination
        self.comb_matrix = hdfcompoundcomb.possiblecombinationsize
        self.header = hdfcompoundcomb.fields

        # ### save original and processed data
        self.originaldata = hdfcompoundcomb.data
        self.stdinh = hdfcompoundcomb.std_inh

        # type of data from hdfcompound object in main used
        # for Dose-response curve
        self.normalized_perc = hdfcompoundcomb.normalized_perc
        self.inhib_data = hdfcompoundcomb.inhibited
        self.inhibitedperc = hdfcompoundcomb.inhibitedperc
        self.inhib_std = hdfcompoundcomb.std_inh
        self.stdnormalized_perc = hdfcompoundcomb.std
        self.normtozero = hdfcompoundcomb.normtozero

        self.branch_size = []
        self.control_info = []
        # for synergy
        self.curvessynergy = []
        self.concentration_rangeforzip = []

        #  ####
        self.folder = None
        self.pdf_pages = None
        self.plot_titel = None
        self.std_type = None
        self.time_position = None
        self.time_selected = None
        self.comb_list_group = None
        self.comalone_list_group = None
        self.conflue_data = None
        self.data_plot = None
        self.std_plot = None
        self.data_plotnorm = None
        self.std_plotnorm = None
        self.concentration_range = None
        self.confluency_range = None
        self.std_range = None
        self.plotting = None

        # for synergy
        self.comb_name_per_group = None
        self.synergy_score_per_group = None
        # self.synergy_score_per_groupdifference = None
        self.synergy_score_per_groupycyzip = None
        self.effect_per_group = None
        self.doseresponsecurve = None

        # ## this is defined from commonfunctions class AdjustData
        self.adjusteddata = None  # this is to set all values <0 or >100 as 0, important for synergism calculations
        self.startpdf = None

        #  from fileobject object in main
        # files names to save data
        self.txt_path = fileobject.fileictxtresultspath
        self.icgrowthratetxt_path = fileobject.fileicgrtxtresultspath
        self.pdf_path = fileobject.filepdfresultspath

        self.blisspredictpath = fileobject.filepredictedblisscorpath
        self.blissscorepath = fileobject.fileblisscorpath

        self.yLoewepath = fileobject.fileyloewetxtpath
        self.fileyloewedifferencetxtpath = fileobject.fileyloewedifferencetxtpath
        self.fileloeweCItxtpath = fileobject.fileloeweCItxtpath

        self.yZippath = fileobject.fileyZiptxtpath
        self.fileyZipdifferencetxtpath = fileobject.fileyZipdifferencetxtpath
        self.fileZipCItxtpath = fileobject.fileZipCItxtpath

        self.fileHsamaximumpath = fileobject.fileHsamaximumpath
        self.fileHsascorepath = fileobject.fileHsascorepath

        # get from doseresponse script
        self.curves = None
        self.curvesgr = None
        self.curveszip = None
        self.curves_specifictimepoints = None
        self.curvessynergynormal = None

        # Grwoth rate
        self.grresults = None
        self.time_selectedgr = None
        self.grheader = None
        self.position_gr = []
        self.position_elapsegr = []
        self.grselectedelap = []
        self.grselectedtimepoint = []
        if len(self.elapse) >1:
            self.extremepoints = [self.elapse[1], self.elapse[-2]]
            self.extremepointspositions = []

        # get main time point
        timeselect = Timepointselection(self.elapse)
        self.time_selected = timeselect.time_selected
        self.time_position = timeselect.time_position

        self.startcommonfunctions = CommonFunctions(self)

        self.row_num = self.startcommonfunctions.row_num
        self.colun_num = self.startcommonfunctions.colun_num

        # collect information of combination--> fiels and then groups

        self.make_empty_list(self.possible_comb)
        self.make_empty_list_single()
        self.get_headers_single()
        self.get_headers_comb()  # will get the list with possible combination

        self.comb_name_per_group = [[] for i in range(0, len(self.com_list_group))]
        self.synergy_score_per_group = [[] for i in range(0, len(self.com_list_group))]
        self.effect_per_group = [[] for i in range(0, len(self.com_list_group))]
        self.synergy_score_per_yzipref = [[] for i in range(0, len(self.com_list_group))]
        self.synergy_score_per_yzipfit = [[] for i in range(0, len(self.com_list_group))]
        self.combinationzipname = [[] for i in range(0, len(self.com_list_group))]

        # compute dose-response curve for synergy
        if len(self.elapse) > 1:
            gr = GrowthRateCombination(self, 1)
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
            self.grresultscombination = gr.grresultscombination
            self.grheadercombination = gr.grheadercombination
            txtgr = GrwothRateSaveTXT(self)
            txtgr.TXTcombination(self.grresultscombination, self.grheadercombination)

        self.doseresponsecurve = DoseResponse(self)

        for i in range(len(self.com_list_group)):
            for j in range(len(self.com_list_group[i])):
                # cell combination
                if self.synergy_method != 2:
                    self.adjusteddata = AdjustData(self)
                    self.normalized_perc = self.adjusteddata.normalized_perc
                    self.inhib_data = self.adjusteddata.inhib_data
                    self.inhibitedperc = self.adjusteddata.inhibitedperc

                    grup = [Groupping(c, self.possible_comb[i][j], self.branch[i][0], self.seeding,
                                      self.condition) for c in self.com_list_group[i][j][0]]
                    data_grup = [Data_group(self.inhib_data, i.concentration) for i in grup]

                    if len(grup) >= 3 and self.synergy_method == 0:
                        biscore = Blissmethod(grup, data_grup, i, j)
                        self.effect_per_group[i].append(biscore.predvalues)
                        self.synergy_score_per_group[i].append(biscore.bivalues)
                        self.comb_name_per_group[i].append(biscore.combname)

                    elif len(grup) >= 3 and self.synergy_method == 1:
                        print('HSA method')
                        hsascore = Hsamethod(grup, data_grup, i, j)
                        self.effect_per_group[i].append(np.zeros(hsascore.maximumeffect.shape))
                        self.synergy_score_per_group[i].append(hsascore.hsavalues)
                        self.comb_name_per_group[i].append(hsascore.combname)
                else:
                    if len(self.possible_comb[i][j]) > 3:
                        new = list(combinations(self.possible_comb[i][j][1:], 2))
                        new2 = ['_'.join(i) for i in new]
                        newcom = [self.possible_comb[i][j][0]]
                        newcom += new2
                        newcom += self.possible_comb[i][j][1:]
                        grup = [Grouppingzip(c, self.header, self.branch[i][0], self.seeding,
                                          self.condition) for c in newcom]
                        data_grup = [Data_group(self.inhib_data, i.concentration) for i in grup]
                    else:
                        grup = [Groupping(c, self.possible_comb[i][j], self.branch[i][0], self.seeding,
                                          self.condition) for c in self.com_list_group[i][j][0]]
                        data_grup = [Data_group(self.inhib_data, i.concentration) for i in grup]

                    if grup[-1].concentration.shape[0] > 2 and len(grup)>=3:
                        # get curves for all time points, and
                        zip = ZipMethomultidimen(self, grup, data_grup)
                        self.synergy_score_per_group[i].append(np.asarray(zip.zipvalues))
                        self.effect_per_group[i].append(np.asarray(zip.yzipvalues))
                        self.comb_name_per_group[i].append(zip.combname)
                        print('check after ZIP')
                    else:
                        self.effect_per_group[i].append([])
                        self.synergy_score_per_group[i].append([])
                        self.comb_name_per_group[i].append([])
                        self.combinationzipname[i].append([])

        print('check')
        SaveTXTinfo(self, self.fileioriginaldatapath, self.originaldata)
        SaveTXTinfo(self, self.filenormalizedatapath, self.normalized_perc)
        SaveTXTinfo(self, self.fileinhibiteddatapath, self.inhib_data)

        if self.synergy_method == 0:
            SynergyScoreSaveTXT(self, self.blisspredictpath, self.effect_per_group)
            SynergyScoreSaveTXT(self, self.blissscorepath, self.synergy_score_per_group)

        elif self.synergy_method == 1:
            SynergyScoreSaveTXT(self, self.fileHsamaximumpath, self.effect_per_group)
            SynergyScoreSaveTXT(self, self.fileHsascorepath, self.synergy_score_per_group)

        elif self.synergy_method == 2:
            SynergyScoreSaveTXT(self, self.yZippath, self.effect_per_group)
            SynergyScoreSaveTXT(self, self.fileyZipdifferencetxtpath, self.synergy_score_per_group)

        self.startpdf = OpenPdf(self)
        self.pdf_path = self.startpdf.pdf_path
        self.pdf_pages = self.startpdf.pdf_pages
        self.plotting = Overtime(self, self.synergy_method)

    def make_empty_list(self, possible_comb):
        header_count = len(possible_comb[0])
        list1 = []
        for i in range(len(possible_comb)):
            headers = [[] for i in range(0, header_count)]  # list_possible_combination
            list1.append(headers)

        self.com_list_group = list1

        return

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

    def get_headers_comb(self):

        for i in range(len(self.possible_comb)):
            # cell line
            for j in range(len(self.possible_comb[i])):
                # combination
                list2 = []
                for k in range(len(self.possible_comb[i][j])):
                    list1 = []
                    for u in self.header:
                        if u[0] == self.branch[i][0] and u[3] == self.possible_comb[i][j][k]:
                            list1.append(u)
                    list2.append(list1)

                self.com_list_group[i][j].append(list2)

    def close_pdf(self):
        self.startpdf.close_pdf()

    def doseresponse(self, growthrate=0):
        self.data_plot = self.normalized_perc
        for i in range(len(self.comalone_list_group)):
            self.control_info = self.control[i]
            for j in range(len(self.comalone_list_group[i])):
                grup = [Groupping(c, self.compound_alone[i][j], self.celine[i], self.seeding,
                                  self.condition) for c in self.comalone_list_group[i][j]]
                if len(grup[0].concentration) > 2:
                    print(grup[0].name)
                    self.concentration_range = self.doseresponsecurve.get_concentration_range(grup[0])
                    datanormperce, std = self.doseresponsecurve.get_confluency_range(self, grup[0])
                    self.doseresponsecurve.get_curve_specifictimepoints(self, datanormperce, grup[0])
                    self.curves = self.doseresponsecurve.curves
                    DoseResponseSinglePlotting(self, grup[0], datanormperce, std)

                    # # Growth rate
                    if growthrate ==0:
                        datagr, stdgr = self.doseresponsecurve.get_confluency_rangeGR(self, grup[0], i, j)
                        self.doseresponsecurve.get_curve_specifictimepointsGR(self, datagr, grup[0])
                        self.curvesgr = self.doseresponsecurve.curves
                        DoseResponseGrowthRate(self, grup[0], i, j, datagr, stdgr)

    def inhibition(self, data, std, info1=0, info2=0):
        self.data_plot = data
        self.data_plotgr = self.normtozero
        self.std_plot = std

        for i in range(len(self.comalone_list_group)):
            # cell line
            self.control_info = self.control[i]
            for j in range(len(self.comalone_list_group[i])):
                grup = [Groupping(c, self.compound_alone[i][j], self.celine[i], self.seeding,
                                  self.condition) for c in self.comalone_list_group[i][j]]

                self.plotting.plot_single_compond(grup[0], self, info1, info2)

            for k in range(len(self.com_list_group[i])):
                grup = [Groupping(c, self.possible_comb[i][k], self.celine[i], self.seeding,
                                  self.condition) for c in self.com_list_group[i][k][0]]

                grup_orig = grup.copy()
                grup_grwt = grup_orig.copy()
                grup_bar = grup.copy()
                data_grup = [Data_group(self.data_plot, m.concentration) for m in grup]
                data_orig = data_grup.copy()
                data_bar = data_grup.copy()
                std_data_grup = [Data_group(self.std_plot, m.concentration) for m in grup]
                std_orig = std_data_grup.copy()
                std_bar = std_data_grup.copy()

                if self.synergy_method == 2:
                    if len(grup[-1].concentration) > 1:
                        self.plotting.plot_combin_inhibition_recursive(self, grup, data_grup, std_data_grup, i, k,
                                                                       info1, info2)
                        if info1 != 0:
                            self.plotting.barplot_comb_inhibition_recursive(self, grup_bar, data_bar, std_bar, info1,
                                                                            info2, i, k)

                            self.plotting.heat_map_overtime_bi(self, i, k, grup_orig)
                            self.plotting.heat_map_selec_time_bi_biDim(self, grup_orig, i, k)
                else:
                    # self.plotting.plot_combin_inhibition_recursive(self, grup, data_grup, std_data_grup, i, k,
                    #                                                info1, info2)
                    if info1 == 0:
                        self.plotting.plot_grwothratecombination(self, grup_grwt, i, k)
                    if info1 != 0:
                        self.plotting.plot_combin_inhibition_recursive(self, grup, data_grup, std_data_grup, i, k,
                                                                       info1, info2)
                        self.plotting.heat_map_overtime_bi(self, i, k, grup_orig)
                        self.plotting.heat_map_selec_time_bi_biDim(self, grup_orig, i, k)
                        if info2 != 0:
                            self.plotting.barplot_comb_inhibition_recursive(self, grup_bar, data_bar, std_bar, info1,
                                                                            info2, i, k)


    def endpointinhibition(self, data, std, info1, info2):
        self.data_plot = data
        self.std_plot = std

        for i in range(len(self.comalone_list_group)):
            # cell line
            self.control_info = self.control[i]

            for k in range(len(self.com_list_group[i])):
                grup = [Groupping(c, self.possible_comb[i][k], self.celine[i], self.seeding,
                                  self.condition) for c in self.com_list_group[i][k][0]]

                grup_orig = grup.copy()
                grup_bar = grup.copy()
                data_grup = [Data_group(self.data_plot, m.concentration) for m in grup]
                data_orig = data_grup.copy()
                data_bar = data_grup.copy()
                std_data_grup = [Data_group(self.std_plot, m.concentration) for m in grup]
                std_orig = std_data_grup.copy()
                std_bar = std_data_grup.copy()

                # Commonprocessing.adjust_data_new(data_grup)
                if self.synergy_method == 2:
                    if len(grup[-1].concentration) > 1:
                        self.plotting.barplot_comb_inhibition_recursive(self, grup_bar, data_bar, std_bar, info1, info2,
                                                                        i, k)
                        if len(grup) == 3:
                            self.plotting.heat_map_selec_time_bi_biDim(self, grup_orig, i, k)
                else:
                    self.plotting.barplot_comb_inhibition_recursive(self, grup_bar, data_bar, std_bar, info1, info2,
                                                                    i, k)
                    if len(grup) == 3:
                        self.plotting.heat_map_selec_time_bi_biDim(self, grup_orig, i, k)

                # if self.comb_matrix[i][k][-1] == self.comb_matrix[i][k][-2] and len(grup_orig) == 3 and \
                #         len(grup_orig[1].concentration) > 2 and len(grup_orig[-1].concentration) > 2:
                # self.plotting.heat_map_selec_time_bi_biDim(grup_orig, self.time_selected, self.time_position,
                #                                            self.synergy_score_per_group[i][k])
                #
                # self.plotting.heat_map_overtime_bi(self, self.synergy_score_per_groupdifference, i, k, grup_orig)
                # if len(grup[-1].concentration) >1:
                #     self.plotting.heat_map_selec_time_bi_biDim(grup_orig, self.time_selected, self.time_position,
                #                                                self.synergy_score_per_groupdifference[i][k],
                #                                                np.round(self.synergy_score_per_groupdifference[i][k],
                #                                                         decimals=3),
                #                                                self.comb_name_per_group[i][k])

                    # self.plotting.heat_map_selec_time_bi_biDim(grup_orig, self.time_selected, self.time_position,
                    #                                            self.synergy_score_per_yzipref[i][k],
                    #                                            self.comb_name_per_group[i][k])
