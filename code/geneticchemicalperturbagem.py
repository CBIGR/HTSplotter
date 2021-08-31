import numpy as np
import os
from grupping import Groupping, Data_group, Grouppingmedium, Grouppingpertub
from synergism import Blissmethod
from plotting import Overtime
from timepointselection import Timepointselection
from copy import deepcopy
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm
from matplotlib.backends.backend_pdf import PdfPages


import math
import scipy.interpolate as inter
from doseresponse import DoseResponseSingle, DoseResponseCondition
import time

import time as ti
class GeneticChemicalPerturbagem:

    def __init__(self, fields, elapse, branch, control, compound_alone,
                 possible_combinations, data_range_cpossi_combin, inhib, std_inh, celline, seeding, condition,
                 std_type, fieldsmedium, fieldsmediuminhib, datamedium,normalizedtranslationmedium, inhibitedmedium, stdmedium,
                 std_inhmedium, information_readout, readout_units):

        self.experiment_name = None
        self.folder = None
        self.pdf_path = None
        self.pdf_pages = None
        self.txt_path = None
        self.readout = information_readout
        self.readout_units = readout_units

        self.plot_titel = None
        self.elapse = elapse
        self.std_type = None
        self.stdinfo = std_type
        self.medium = fieldsmedium
        self.fieldsmediuminhib = fieldsmediuminhib
        self.datamedium = datamedium
        self.normalizedtranslationmedium = normalizedtranslationmedium
        self.inhibitedmedium = inhibitedmedium
        self.stdmedium = stdmedium
        self.std_inhmedium = std_inhmedium
        self.time_position = None
        self.time_selected = None

        self.branch = branch[::2]
        self.brachcomp = branch
        self.branch_size = []
        self.control = control
        self.control_info = []
        self.celine = celline
        self.condition = condition
        self.seeding = seeding
        self.compound_alone = compound_alone
        self.possible_comb = possible_combinations
        self.onandoff = None
        self.comb_matrix = data_range_cpossi_combin
        self.header = fields
        self.comb_list_group = None
        self.comalone_list_group = None
        # files names to save data
        self.txt_path = None
        self.experiment_name = None
        self.pdf_path = None

        self.inhib_data = inhib
        self.inhib_std = std_inh
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

        self.row_num = None
        self.colun_num = None
        self.plotting = None

        # get main time point
        timeselect = Timepointselection(self.elapse)
        self.time_selected = timeselect.time_selected
        self.time_position = timeselect.time_position

        # collect information of combination--> fiels and then groups
        self.get_control()
        self.make_empty_list(self.possible_comb)
        self.make_empty_list_single()
        self.get_headers_single()
        self.get_headers_comb()  # will get the list with possible combination

        self.comb_name_per_group = [[] for i in range(0, len(self.com_list_group))]
        self.bi_score_per_group = [[] for i in range(0, len(self.com_list_group))]
        self.predicted_per_group = [[] for i in range(0, len(self.com_list_group))]

        for i in range(len(self.com_list_group)):
            for j in range(len(self.com_list_group[i])):
                # cell combination
                grup = [Groupping(c, self.possible_comb[i][j], self.celine[i], self.seeding,
                                  self.condition) for c in self.com_list_group[i][j][0]]
                data_grup = [Data_group(self.inhib_data, m.concentration) for m in grup]
                if len(grup) >= 3:
                    biscore = Blissmethod(grup, data_grup, i, j)
                    self.predicted_per_group[i].append(biscore.predvalues)
                    self.bi_score_per_group[i].append(biscore.bivalues)
                    self.comb_name_per_group[i].append(biscore.combname)

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

                #self.onandoff[i][j].append(deepcopy(list2))
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

    def get_row_plot_pag(self):
        self.row_num = math.ceil(len(self.time_selected)/3)
        self.colun_num = 3

        return

    def get_concentration_range(self, grup):
        concentration_range = [0.001]  # this value is for DMSO control

        for j in range(len(grup.concentration)):
            concentration_range.append(grup.concentration[j][0])

        self.concentration_range = concentration_range
        return

    def get_confluency_range(self, grup):
        confluency_range = []  # this value is for DMSO control
        std_range = []
        for k in range(len(self.time_position)):
            tmp1 = [self.data_plot[self.control_info[-1]][self.time_position[k]]]  # this value is for DMSO control
            tmp2 = [self.std_plot[self.control_info[-1]][self.time_position[k]]]
            for j in range(len(grup.concentration)):
                pos = int(grup.concentration[j][-1])
                tmp1.append(self.data_plot[pos][self.time_position[k]])
                tmp2.append(self.std_plot[pos][self.time_position[k]])
            confluency_range.append(tmp1)
            std_range.append(tmp2)
        self.confluency_range = np.asarray(confluency_range)
        self.std_range = np.asarray(std_range)

        return

    def get_ic_txt_path(self, experiment_name):
        self.txt_path = experiment_name
        f = open(experiment_name, 'w')
        f.write('This is 1 biological replicate. Be careful drawing any conclusions.' + '\n')
        f.close()

    def bi_score_save_txt(self, name, score):
        bi_txt_path = name
        f = open(bi_txt_path, 'w')

        for i in range(len(self.possible_comb)):
            for k in range(len(self.possible_comb[i])):
                ela = np.asarray(["Compound=>" + str(self.possible_comb[i][k])])
                elapsed_cabe = np.hstack((ela, self.elapse))
                temp_name = np.asarray(self.comb_name_per_group[i][k])

                name_data = np.concatenate((temp_name[None, :], score[i][k]))
                new = np.hstack((elapsed_cabe[:, None], name_data))
                np.savetxt(f, new, delimiter="\t", fmt="%s")

        f.close()

    def open_pdf(self, pdfpathname, experiment_name):
        self.experiment_name = experiment_name
        self.pdf_path = pdfpathname
        self.pdf_pages = PdfPages(self.pdf_path)
        self.plotting = Overtime(self.pdf_pages, self.experiment_name, self.elapse, self.stdinfo)

    def close_pdf(self):
        self.pdf_pages.close()

    def confluencyovertime(self, data, std, info1, info2):
        self.data_plot = data
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
                                self.plotting.plot_single_compond(grup[0], self.control_info, self.data_plot,
                                                                  self.std_plot, info1, info2,
                                                                  self.readout, self.readout_units)
                        else:
                            if h[3].split("_")[-1] in self.comalone_list_group[i][j][k][0][0][3]:
                                self.control_info = h
                                grup = [Groupping(c, self.compound_alone[i][j], self.branch[i][0], self.branch[i][1],
                                                  self.branch[i][2]) for c in self.comalone_list_group[i][j][k]]
                                self.plotting.plot_single_compond(grup[0], self.control_info, self.data_plot,
                                                                  self.std_plot, info1, info2,
                                                                  self.readout, self.readout_units)
        for g in range(len(self.com_list_group)):
            for k in range(len(self.com_list_group[g])):
                for h in self.control[g]:
                    if h[3].split("_")[-1] != self.branch[g][2]:
                        if h[3].split("_")[-1] == self.com_list_group[g][k][0][1][0][3]:
                            self.control_info = h
                            grup = [Groupping(c, self.possible_comb[g][k], self.celine[g], self.seeding,
                                              self.condition) for c in self.com_list_group[g][k][0]]

                            grup_orig = grup.copy()
                            grup_bar = grup.copy()
                            data_grup = [Data_group(self.data_plot, m.concentration) for m in grup]
                            data_orig = data_grup.copy()
                            data_bar = data_grup.copy()
                            std_data_grup = [Data_group(self.std_plot, m.concentration) for m in grup]
                            std_orig = std_data_grup.copy()
                            std_bar = std_data_grup.copy()
                            self.plotting.plot_combin_inhibition_recursive(grup, self.control_info, self.data_plot,
                                                                           data_grup,
                                                                           std_data_grup, info1, info2,
                                                                           self.comb_name_per_group[g][k],
                                                                           self.predicted_per_group[g][k],
                                                                           self.bi_score_per_group[g][k],
                                                                           self.time_selected,
                                                                           self.time_position, self.readout,
                                                                           self.readout_units)
                            if info1 != 0:
                                self.plotting.heat_map_overtime_bi(grup_orig, self.bi_score_per_group[g][k],
                                                                   self.comb_name_per_group[g][k])

        if info1 == 0:
            posi = 0
            for k in range(len(self.medium[::2])):
                grup = [Grouppingpertub(c, self.branch[k][0], self.branch[k][1],
                                        self.branch[k][2]) for c in self.medium[posi:posi+2]]
                posi +=2
                for j in range(len(grup)):
                    self.plotting.plot_control_medium(grup[j], self.data_plot, self.std_plot, self.datamedium,
                                                      self.stdmedium, self.control[k], 0, 0, self.readout,
                                                      self.readout_units)
                self.plotting.plot_medium(grup, self.datamedium, self.stdmedium, 0, 0, self.readout,
                                          self.readout_units)

    def doseresponse(self, data, std):
        self.data_plot = data
        self.std_plot = std

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
                            if len(grup[0].concentration) > 2:
                                for h in range(len(grup)):
                                    self.get_row_plot_pag()
                                    # get x label => concentration range + lowest value for DMSO
                                    self.get_concentration_range(grup[h])
                                    self.get_confluency_range(grup[h])
                                    dr = DoseResponseCondition(self.data_plot, self.std_range, grup[h], self.row_num,
                                                               self.concentration_range,
                                                               self.confluency_range, self.time_position,
                                                               self.time_selected, self.txt_path, self.pdf_pages,
                                                               self.readout, self.readout_units)
                                    dr.singlecompoundplot(grup[h], self.compounds_name_on_off, self.dose_curve_on_off,
                                                          self.data_curve_on_off, self.std_curve_on_off,
                                                          self.ic_concentration_on_off,
                                                          self.concentration_value_on_off)
                                    if h == 1:
                                        dr.sigmoid_on_off_summary(grup[h], self.compounds_name_on_off,
                                                                  self.dose_curve_on_off,
                                                                  self.data_curve_on_off, self.std_curve_on_off,
                                                                  self.ic_concentration_on_off,
                                                                  self.concentration_value_on_off,
                                                                  self.concentration_range)

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

                            self.plotting.barplot_comb_inhibition_recursive(grup_bar, data_bar, std_bar, info1, info2,
                                                                            self.comb_name_per_group[g][k],
                                                                            self.predicted_per_group[g][k],
                                                                            self.bi_score_per_group[g][k],
                                                                            self.time_selected, self.time_position,
                                                                            self.readout, self.readout_units)
                            self.plotting.heat_map_selec_time_bi_biDim(grup_orig, self.time_selected,
                                                                       self.time_position,
                                                                       self.bi_score_per_group[g][k])


        # for cel in range(len(self.onandoff)):
        #     for cond in range(len(self.onandoff)):
        #         grup = [Groupping(c, self.onandoff[cel][cond], self.celine[cel], self.seeding,
        #                           self.condition) for c in self.onandoff[cel][cond]]
        #         self.plotting.barplot_geneticchemicalperturbgen_recursive(grup, data, std, 2, 1,
        #                                                                   self.time_selected, self.time_position,
        #                                                                   self.readout,
        #                                                                   self.readout_units)

