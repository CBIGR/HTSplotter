import numpy as np
import os
from HTSplotter.grupping import Groupping, Data_group, GrouppingGrowthRate, Grouppingpertub
from HTSplotter.synergism import Blissmethod, Hsamethod
from HTSplotter.plotting import Overtime
from HTSplotter.timepointselection import Timepointselection
from HTSplotter.commonfunctions import CommonFunctions, OpenPdf, SaveTXTinfo, AdjustData, SynergyScoreSaveTXT, GrwothRateSaveTXT
from HTSplotter.doseresponse import DoseResponse
from copy import deepcopy
from HTSplotter.doseresponseplots import DoseResponseConditionPlotting, DoseResponseGrowthRate
from HTSplotter.growthrate import GrowthRate, GrowthRateGeneticChemical
# import matplotlib.pyplot as plt
# from matplotlib.pyplot import cm
# from matplotlib.backends.backend_pdf import PdfPages
# import math
# import scipy.interpolate as inter
#
# import time
#
# import time as ti


class GeneticChemicalPerturbagem:

    def __init__(self, synergymethod, hdfgeneticchemical, file_info, file_names, information_readout,
                 readout_units, i, biologicalreplicate):

        self.synergy_method = synergymethod
        self.readout = information_readout
        self.readout_units = readout_units
        self.experiment_name = i
        self.biologicalreplicateinfo = biologicalreplicate

        self.branch = hdfgeneticchemical.branch[::2]
        self.brachcomp = hdfgeneticchemical.branch
        self.control = hdfgeneticchemical.control
        self.celine = hdfgeneticchemical.celline
        self.condition = hdfgeneticchemical.condition
        self.seeding = hdfgeneticchemical.seeding
        self.compound_alone = hdfgeneticchemical.compoundalone
        self.possible_comb = hdfgeneticchemical.possiblecombination
        self.comb_matrix = hdfgeneticchemical.possiblecombinationsize
        self.header = hdfgeneticchemical.fields
        self.stdinfo = hdfgeneticchemical.std_info

        self.elapse = file_info.elapsed
        self.date_info = file_info.date_info

        self.originaldata = hdfgeneticchemical.data
        self.inhib_data = hdfgeneticchemical.inhibited
        self.stdinh = hdfgeneticchemical.std_inh
        self.normtozero = hdfgeneticchemical.normtozero

        # for Dose-response curve
        self.normalized_perc = hdfgeneticchemical.normalized_perc
        self.inhibitedperc = hdfgeneticchemical.inhibitedperc
        self.stdnormalized_perc = hdfgeneticchemical.std
        self.inhib_std = hdfgeneticchemical.std_inh

        self.medium = hdfgeneticchemical.fieldsmedium
        self.fieldsmediuminhib = hdfgeneticchemical.fieldsmediuminhibited
        self.mediumdata = hdfgeneticchemical.datamedium
        self.normalizedtranslationmedium = hdfgeneticchemical.normalizedtranslationmedium
        self.inhibitedmedium = hdfgeneticchemical.inhibitedmedium
        self.mediumstd = hdfgeneticchemical.stdmedium
        self.std_inhmedium = hdfgeneticchemical.std_inhmedium

        self.control_info = []
        self.branch_size = []
        #  from fileobject object in main
        # files names to save data
        self.fileioriginaldatapath = file_names.fileioriginaldatapath
        self.fileinhibiteddatapath = file_names.fileinhibiteddatapath
        self.filenormalizedatapath = file_names.filenormalizedatapath
        self.filegrowthrateresults = file_names.filegrowthrateresults
        self.txt_path = file_names.fileictxtresultspath
        self.icgrowthratetxt_path = file_names.fileicgrtxtresultspath
        self.pdf_path = file_names.filepdfresultspath

        self.blisspredictpath = file_names.filepredictedblisscorpath
        self.blissscorepath = file_names.fileblisscorpath

        self.fileHsamaximumpath = file_names.fileHsamaximumpath
        self.fileHsascorepath = file_names.fileHsascorepath

        self.std_type = None
        self.plot_titel = None
        self.onandoff = None
        self.comalone_list_group = None
        self.com_list_group = None

        # files names to save data

        self.pdf_pages = None
        self.conflue_data = None
        self.data_plot = None
        self.std_plot = None
        self.concentration_range = None
        self.confluency_range = None
        self.std_range = None

        self.compounds_name_on_off = None
        self.dose_curve_on_off = None
        self.data_curve_on_off = None
        self.std_curve_on_off = None
        self.ic_concentration_on_off = None
        self.concentration_value_on_off = None
        self.relative_IC50concentration_on_off = None

        self.compounds_name_on_off_gr = None
        self.dose_curve_on_off_gr = None
        self.data_curve_on_off_gr = None
        self.std_curve_on_off_gr = None
        self.ic_concentration_on_off_gr = None
        self.concentration_value_on_off_gr = None
        self.relative_IC50concentration_on_off_gr = None

        self.row_num = None
        self.colun_num = None
        self.plotting = None

        self.time_position = None
        self.time_selected = None

        # get from doseresponse script
        self.curves = None
        self.curves_specifictimepoints = None
        self.doseresponsecurve = None

        # growth rate
        self.grresults = None
        self.time_selectedgr = None
        self.grheader = None
        self.position_gr = []
        self.position_elapsegr = []
        self.grselectedelap = []
        self.grselectedtimepoint = []
        self.grresultscombination = None
        self.grheadercombination = None
        self.curvesgr = None
        self.conditionorder = [[] for i in range(len(self.possible_comb))]
        self.comb_list_group = [[] for i in range(len(self.possible_comb))]

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
        self.adjusteddata = AdjustData(self)
        # collect information of combination--> fiels and then groups
        self.get_control()
        self.make_empty_list(self.possible_comb)
        self.make_empty_list_single()
        self.get_headers_single()

        self.get_headers_comb()  # will get the list with possible combination

        self.comb_name_per_group = [[] for i in range(0, len(self.com_list_group))]
        self.synergy_score_per_group = [[] for i in range(0, len(self.com_list_group))]
        self.effect_per_group = [[] for i in range(0, len(self.com_list_group))]


        for i in range(len(self.com_list_group)):
            for j in range(len(self.com_list_group[i])):
                # cell combination
                print(self.possible_comb[i][j])
                grup = [Groupping(c, self.possible_comb[i][j], self.celine[i], self.seeding,
                                  self.condition) for c in self.com_list_group[i][j][0]]
                data_grup = [Data_group(self.inhib_data, m.concentration) for m in grup]
                if self.synergy_method == 0:
                    biscore = Blissmethod(grup, data_grup, i, j)
                    self.effect_per_group[i].append(biscore.predvalues)
                    self.synergy_score_per_group[i].append(biscore.bivalues)
                    self.comb_name_per_group[i].append(biscore.combname)

                elif self.synergy_method == 1:
                    print('HSA method')
                    hsascore = Hsamethod(grup, data_grup, i, j)
                    self.effect_per_group[i].append(np.zeros(hsascore.maximumeffect.shape))
                    self.synergy_score_per_group[i].append(hsascore.hsavalues)
                    self.comb_name_per_group[i].append(hsascore.combname)

        SaveTXTinfo(self, self.fileioriginaldatapath, self.originaldata)
        SaveTXTinfo(self, self.filenormalizedatapath, self.normalized_perc)
        SaveTXTinfo(self, self.fileinhibiteddatapath, self.inhib_data)

        if self.synergy_method == 0:
            SynergyScoreSaveTXT(self, self.blisspredictpath, self.effect_per_group)
            SynergyScoreSaveTXT(self, self.blissscorepath, self.synergy_score_per_group)

        elif self.synergy_method == 1:
            SynergyScoreSaveTXT(self, self.fileHsamaximumpath, self.effect_per_group)
            SynergyScoreSaveTXT(self, self.fileHsascorepath, self.synergy_score_per_group)
            # compute dose-response curve for synergy
        if len(self.elapse) > 1:
            gr = GrowthRateGeneticChemical(self)
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
        self.startpdf = OpenPdf(self)
        self.pdf_path = self.startpdf.pdf_path
        self.pdf_pages = self.startpdf.pdf_pages
        self.plotting = Overtime(self, self.synergy_method)

    def get_control(self):
        for j in self.header:
            for i in range(len(self.control)):
                for k in range(len(self.control[i])):
                    # per cell line
                    if self.control[i][k][0] == j[0] and self.control[i][k][1] == j[1] and \
                            self.control[i][k][2] == j[2] and self.control[i][k][3] == j[3]:
                        self.control[i][k] = j
        # self.combcontrol = self.control[::2]
        # self.doxcontrol = self.control[1::2]

    def make_empty_list(self, possible_comb):
        header_count = len(possible_comb[0])
        list1 = []
        for i in range(len(possible_comb)):
            headers = [[] for i in range(0, header_count)]  # list_possible_combination
            list1.append(headers)
        self.com_list_group = list1
        self.onandoff = deepcopy(list1)
        self.conditionorder.append([])
        for k in range(len(possible_comb)):
            for h in range(len(possible_comb[k])):
                first = possible_comb[k][h][0].split('_')
                order = [first[-1], 'Condition']
                self.comb_list_group[k].append(possible_comb[k][h][0:-1])
                self.conditionorder[k].append(order)
                # self.conditionorder.append('Condition')
        return

    def make_empty_list_single(self):
        # for j in range(len(self.compound_alone)):
        header_count = len(self.compound_alone)
        headers = [[] for i in range(0, header_count)]  # list_possible_combination
        list1 = headers
        for k in range(len(self.compound_alone)):
            for i in range(len(self.compound_alone[k])):
                headers2 = [[] for i in range(0, len(self.compound_alone[k][i]))]
                list1[k].append(headers2)
        self.comalone_list_group = list1

    def get_headers_single(self):
        for i in range(len(self.compound_alone)):
            for k in range(len(self.compound_alone[i])):
                for t in range(len(self.compound_alone[i][k])):
                    list1 = []
                    for u in self.header:
                        if u[0] == self.branch[i][0] and u[3] == self.compound_alone[i][k][t]:
                            list1.append(u)
                    self.comalone_list_group[i][k][t].append(list1)

    def get_headers_comb(self):
        for i in range(len(self.possible_comb)):
            # cell line
            for j in range(len(self.possible_comb[i])):
                # combination
                list2 = []
                # get the genetic perturbagen condition with solvent
                for u in self.control[i]:
                    nam = u[3].split("_")[1:]
                    if self.possible_comb[i][j][1] in u[3] and u[0] == self.celine[i] and u[1] == self.seeding[i]:
                        nam.append(self.condition[i][-1])
                        self.possible_comb[i][j][-1] = "_".join(nam)

                for k in range(len(self.possible_comb[i][j])):
                    list1 = []
                    for u in self.header:
                        if u[0] == self.celine[i] and u[3] == self.possible_comb[i][j][k]:
                            list1.append(u)
                    list2.append(list1)
                self.com_list_group[i][j].append(list2)

        for h in range(len(self.possible_comb)):
            for d in range(len(self.possible_comb[h])):
                lista4 = []
                for t in range(len(self.possible_comb[h][d][:-1])):
                    lista3 = []
                    for g in self.header:
                        if g[0] == self.celine[h] and g[3] == self.possible_comb[h][d][t]:
                            lista3.append(g)
                    lista4.append(lista3)
                self.onandoff[h][d] = lista4

    def close_pdf(self):
        self.startpdf.close_pdf()

    def confluencyovertime(self, data, std, info1, info2):
        self.data_plot = data
        # self.data_plotgrwot = self.normtozero
        self.std_plot = std
        # # first plot compounds alone
        for i in range(len(self.comalone_list_group)):
            # cell line
            # plot all compounds alone
            for j in range(len(self.comalone_list_group[i])):
                for k in range(len(self.comalone_list_group[i][j])):
                    for h in self.control[i]:
                        if h[3].split("_")[-1] == self.branch[i][2] \
                                and h[0] == self.branch[i][0] and h[1] == self.branch[i][1]:
                            if h[3].split("_")[-2] in self.comalone_list_group[i][j][k][0][0][3]:
                                self.control_info = h

                                grup = [Groupping(c, self.compound_alone[i][j], self.branch[i][0], self.branch[i][1],
                                                  self.branch[i][2]) for c in self.comalone_list_group[i][j][k]]

                        else:
                            if h[3].split("_")[-1] in self.comalone_list_group[i][j][k][0][0][3]:
                                self.control_info = h
                                grup = [Groupping(c, self.compound_alone[i][j], self.branch[i][0], self.branch[i][1],
                                                  self.branch[i][2]) for c in self.comalone_list_group[i][j][k]]
                                self.plotting.plot_single_compond(grup[0], self, info1, info2)

        for g in range(len(self.com_list_group)):
            for k in range(len(self.com_list_group[g])):
                for h in self.control[g]:
                    if h[3].split("_")[-1] != self.branch[g][2]:
                        if h[3].split("_")[-1] == self.com_list_group[g][k][0][1][0][3]:
                            self.control_info = h
                            grup = [Groupping(c, self.possible_comb[g][k], self.celine[g], self.seeding,
                                              self.condition) for c in self.com_list_group[g][k][0]]

                            grup_orig = grup.copy()
                            grup_grwt = grup.copy()
                            # grup_bar = grup.copy()
                            data_grup = [Data_group(self.data_plot, m.concentration) for m in grup]
                            data_orig = data_grup.copy()
                            data_grwt = data_grup.copy()
                            std_data_grup = [Data_group(self.std_plot, m.concentration) for m in grup]
                            std_orig = std_data_grup.copy()
                            std_grwt = std_data_grup.copy()
                            if info1 != 0:
                                self.plotting.plot_combin_inhibition_recursive(self, grup, data_grup, std_data_grup,
                                                                               g, k, info1, info2)
                                self.plotting.heat_map_overtime_bi(self, g, k,
                                                                   grup_orig)
                            if info1 == 0:
                                self.plotting.plot_grwothratecombination(self, grup_grwt, g, k)

        if info1 == 0 and len(self.medium) > 1:
            posi = 0
            mepos = 0
            for k in range(len(self.medium[::2])):
                grup = [Grouppingpertub(c, self.branch[k][0], self.branch[k][1],
                                        self.branch[k][2]) for c in self.medium[posi:posi+2]]

                posi +=2
                for j in range(len(grup)):
                    self.plotting.plot_control_medium(grup[j], self, mepos, k)
                    mepos += 1
                self.plotting.plot_medium(self, grup, 0, 0)

            # posi = 0
            # for k in range(len(self.medium[::2])):
            #     grup = [Grouppingpertub(c, self.branch[k][0], self.branch[k][1],
            #                             self.branch[k][2]) for c in self.medium[posi:posi + 2]]
            #     posi += 2
            #     for j in range(len(grup)):
            #         self.plotting.plot_control_medium(grup[j], self.data_plot, self.std_plot, self.datamedium[j],
            #                                           self.stdmedium[j], self.control[k], 0, 0, self.readout,
            #                                           self.readout_units)
            #     self.plotting.plot_medium(grup, self.datamedium, self.stdmedium, 0, 0, self.readout,
            #                               self.readout_units)

    def doseresponse(self):
        self.data_plot = self.normalized_perc
        self.std_plot = self.stdnormalized_perc
        self.doseresponsecurve = DoseResponse(self)

        for i in range(len(self.onandoff)):
            # self.control_info = self.control[i]
            for j in range(len(self.onandoff[i])):
                for u in self.control[i]:
                    if u[3].split("_")[-1] != self.branch[i][2]:
                        if u[3].split("_")[-1] in self.onandoff[i][j][0][0][3]:
                            self.control_info = u
                            grup = [Groupping(c, self.onandoff[i][j], self.celine[i], self.seeding,
                                              self.condition) for c in self.onandoff[i][j]]
                            self.compounds_name_on_off = []
                            self.dose_curve_on_off = []
                            self.data_curve_on_off = []
                            self.std_curve_on_off = []
                            self.ic_concentration_on_off = []
                            self.concentration_value_on_off = []
                            self.relative_IC50concentration_on_off = []
                            if len(grup[0].concentration) > 2:
                                for h in range(len(grup)):
                                    # self.get_row_plot_pag()
                                    # get x label => concentration range + lowest value for DMSO
                                    self.concentration_range = self.doseresponsecurve.get_concentration_range(grup[h])
                                    # compound, control_data, control_std
                                    datanormperce, std = self.doseresponsecurve.get_confluency_range(
                                        self, grup[h])
                                    # datainhib, stdinhi = self.get_confluency_range(self.normalized_perc, self.stdnormalized_perc, grup)
                                    self.doseresponsecurve.get_curve_specifictimepoints(self, datanormperce, grup[h])
                                    self.curves = self.doseresponsecurve.curves
                                    dr = DoseResponseConditionPlotting(self, grup[h], datanormperce, std)
                                    dr.singlecompoundplot(self, grup[h])
                                    if h == 1:
                                        dr.sigmoid_on_off_summary(grup[h])
        #
        # growth Rate
        # # Growth Rate
        normal = self.condition[0][0].index('Condition')
        if normal == 0:
            dox = normal+1
        if normal !=0:
            dox = normal-1
        for i in range(len(self.possible_comb)):
            # cell line
            # plot all compounds alone

            count = 0
            for j in range(len(self.possible_comb[i])):
                grup_gr = []
                self.compounds_name_on_off = []
                self.dose_curve_on_off = []
                self.data_curve_on_off = []
                self.std_curve_on_off = []
                self.ic_concentration_on_off = []
                self.concentration_value_on_off = []
                self.relative_IC50concentration_on_off = []
                for k in range(len(self.possible_comb[i][j][0:-1])):

                    for h in range(len(self.control[i])):
                        compounname = self.control[i][h][-3].split('_')[-1]
                        if compounname in self.possible_comb[i][j][k]:
                            self.control_info = self.control[i][h]
                            grup = GrouppingGrowthRate(self.possible_comb[i][j][k], self.header,
                                                       self.celine[i], self.conditionorder[i][j][k],
                                                       self.seeding[i])
                            grup_gr.append(grup)
                            if len(grup.concentration) > 2:
                                # if grup.concentration.shape[0] > 2:
                                self.concentration_range = self.doseresponsecurve.get_concentration_range(grup)
                                print('check ants', grup.name, k, i, j, count)
                                datagr, stdgr = self.doseresponsecurve.get_confluency_rangegeneticchemicalGR(self, grup,
                                                                                                             i, j, k)
                                self.doseresponsecurve.get_curve_specifictimepointsGR(self, datagr, grup)
                                self.curvesgr = self.doseresponsecurve.curves
                                # dr = DoseResponseGrowthRate(self, grup, datagr, stdgr)
                                dr = DoseResponseGrowthRate(self, grup, i, k, datagr, stdgr, 1)
                                print('check', count, grup.name, len(grup_gr))
                                if len(grup_gr) == 2:
                                    print('check after if', count, grup.name, len(grup_gr))
                                    dr.sigmoid_on_off_summarygr(grup_gr[0])
                                count += 1

    def endpointinhibition(self, data, std, info1, info2):
        self.data_plot = data
        self.std_plot = std
        # first plot compounds alone
        for g in range(len(self.com_list_group)):
            for k in range(len(self.com_list_group[g])):
                for h in self.control[g]:
                    if h[3].split("_")[-1] != self.branch[g][2]:
                        if h[3].split("_")[-1] == self.com_list_group[g][k][0][1][0][3]:
                            self.control_info = h
                            grup = [Groupping(c, self.possible_comb[g][k], self.celine[g], self.seeding,
                                              self.condition) for c in self.com_list_group[g][k][0]]

                            # data_grup = [Data_group(self.data_plot, m.concentration) for m in grup]
                            # std_data_grup = [Data_group(self.std_plot, m.concentration) for m in grup]
                            grup_orig = grup.copy()
                            grup_bar = grup.copy()
                            data_grup = [Data_group(self.data_plot, m.concentration) for m in grup]
                            data_orig = data_grup.copy()
                            data_bar = data_grup.copy()
                            std_data_grup = [Data_group(self.std_plot, m.concentration) for m in grup]
                            std_orig = std_data_grup.copy()
                            std_bar = std_data_grup.copy()
                            # self.adjust_data(data_grup)

                            self.plotting.barplot_comb_inhibition_recursive(self, grup_bar, data_bar, std_bar, info1,
                                                                            info2, g, k)
                            self.plotting.heat_map_selec_time_bi_biDim(self, grup_orig, g, k)

        # for cel in range(len(self.onandoff)):
        #     for cond in range(len(self.onandoff)):
        #         grup = [Groupping(c, self.onandoff[cel][cond], self.celine[cel], self.seeding,
        #                           self.condition) for c in self.onandoff[cel][cond]]
        #         self.plotting.barplot_geneticchemicalperturbgen_recursive(grup, data, std, 2, 1,
        #                                                                   self.time_selected, self.time_position,
        #                                                                   self.readout,
        #                                                                   self.readout_units)

