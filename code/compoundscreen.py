import numpy as np

import os
from grupping import Groupping
from grupping import Groupping, Data_group, Grouppingmedium
import matplotlib.pyplot as plt
from doseresponse import DoseResponseSingle, DoseResponseCondition
from plotting import Overtime
from timepointselection import Timepointselection

from matplotlib.pyplot import cm
from matplotlib.backends.backend_pdf import PdfPages
import math


class ExperimentSingleCompound:
    def __init__(self, branch, celline, seeding, fields, elapsed, control_name,
                 compound_alone, condition, std_type, medium, mediumdata, mediumstd, mediuminhfiels, mediuminhib,
                 mediuminhibstd, mediumitranslation, information_readout, readout_units):

        # ##Usefule information for each experiment, information gave by the calls==>new_save_data_base:

        self.experiment_name = None
        self.folder = None
        self.pdf_path = None
        self.pdf_pages = None
        self.txt_path = None
        self.readout = information_readout
        self.readout_units = readout_units

        self.plot_titel = None
        self.elapse = elapsed
        self.std_type = None
        self.stdinfo = std_type
        self.time_position = None
        self.time_selected = None
        self.header = fields
        self.medium = medium
        self.mediuminhfiels = mediuminhfiels
        self.mediumdata = mediumdata
        self.mediumstd = mediumstd
        self.inhibitemediumdata = mediuminhib
        self.translationmediumdata = mediumitranslation
        self.inhibitemediumstd = mediuminhibstd
        self.comalone_list_group = None

        self.compound_alone = compound_alone

        self.data_plot = None
        self.std_plot = None

        self.branch = branch
        self.seeding = seeding
        self.celine = celline
        self.condition = condition
        self.control = control_name
        self.plotting = None

        self.concentration_range = None

        # get main time point
        timeselect = Timepointselection(self.elapse)
        self.time_selected = timeselect.time_selected
        self.time_position = timeselect.time_position

        self.make_empty_list_single()
        self.get_headers_single()

    def get_control(self):
        for j in self.header:
            for i in range(len(self.control)):
                for k in range(len(self.control[i])):
                    # per cell line
                    if self.control[i][k][0] == j[0] and self.control[i][k][1] == j[1] and \
                            self.control[i][k][2] == j[2] and self.control[i][k][3] == j[3]:
                        self.control[i][k] = j

    def get_ic_txt_path(self, experiment_name):
        self.txt_path = experiment_name
        f = open(experiment_name, 'w')
        f.write('This is 1 biological replicate. ' + '\n')
        f.close()

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

    def get_row_plot_pag(self):
        self.row_num = math.ceil(len(self.time_selected) / 3)
        self.colun_num = 3

        return

    def get_concentration_range(self, grup):
        concentration_range = [0.001]  # this value is for DMSO control

        for j in range(len(grup[0].concentration)):
            concentration_range.append(grup[0].concentration[j][0])

        self.concentration_range = concentration_range
        return

    def get_confluency_range(self, grup):
        confluency_range = []  # this value is for DMSO control
        std_range = []
        for k in range(len(self.time_position)):
            tmp1 = [self.data_plot[self.control_info[-1]][self.time_position[k]]]  # this value is for DMSO control
            tmp2 = [self.std_plot[self.control_info[-1]][self.time_position[k]]]
            for j in range(len(grup[0].concentration)):
                pos = int(grup[0].concentration[j][-1])
                tmp1.append(self.data_plot[pos][self.time_position[k]])
                tmp2.append(self.std_plot[pos][self.time_position[k]])
            confluency_range.append(tmp1)
            std_range.append(tmp2)
        self.confluency_range = np.asarray(confluency_range)
        self.std_range = np.asarray(std_range)

        return

    # Methods for plotting main open and close pdf file

    def open_pdf(self, pdfpathname, experiment_name):
        self.experiment_name = experiment_name
        self.pdf_path = pdfpathname
        self.pdf_pages = PdfPages(self.pdf_path)
        self.plotting = Overtime(self.pdf_pages, self.experiment_name, self.elapse, self.stdinfo)

    def close_pdf(self):
        self.pdf_pages.close()


class SingleCompound(ExperimentSingleCompound):
    def __init__(self, branch, celline, seeding, fields, elapsed, control_name,
                 compound_alone, condition, std_type, medium, mediumdata, mediumstd, mediuminhfiels,
                 mediuminhib, mediuminhibstd, mediumitranslation, information_readout, readout_units):
        super().__init__(branch, celline, seeding, fields, elapsed, control_name,
                         compound_alone, condition, std_type, medium, mediumdata, mediumstd, mediuminhfiels,
                         mediuminhib, mediuminhibstd, mediumitranslation, information_readout, readout_units)
        self.get_control()
        # this subclass is for single compound screen with several controls

    def confluencyovertime(self, data, std):
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

                self.plotting.plot_single_compond(grup[0], self.control_info, self.data_plot, self.std_plot, 0, 0,
                                                  self.readout, self.readout_units)

            grup = [Groupping(self.medium[i], self.medium[i][0][3], self.celine[i],
                              self.seeding[i], self.condition[i])]

            self.plotting.plot_control_medium(grup[0], self.data_plot, self.std_plot, self.mediumdata,
                                              self.mediumstd, self.control[i], 0, 0, self.readout,
                                              self.readout_units)

    def doseresponse(self, data, std):
        self.data_plot = data
        self.std_plot = std

        for i in range(len(self.comalone_list_group)):
            for j in range(len(self.comalone_list_group[i])):
                self.control_info = self.control[i][j]
                grup = [Groupping(c, self.compound_alone[i][j], self.celine[i], self.seeding,
                                  self.condition) for c in self.comalone_list_group[i][j]]
                if len(grup[0].concentration) > 2:
                    self.get_row_plot_pag()
                    # # get x label => concentration range + lowest value for DMSO
                    self.get_concentration_range(grup)
                    self.get_confluency_range(grup)
                    DoseResponseSingle(self.data_plot, self.std_range, grup[0], self.row_num,
                                       self.concentration_range,
                                       self.confluency_range, self.time_position, self.time_selected,
                                       self.txt_path,
                                       self.pdf_pages, self.readout, self.readout_units)

    def inhibitionovertime(self, data, std, info1, info2):
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
                self.plotting.plot_single_compond(grup[0], self.control_info, self.data_plot, self.std_plot, info1,
                                                  info2, self.readout, self.readout_units)
            try:
                for j in range(len(self.mediuminhfiels[i])):
                    self.control_info = self.control[i][j]
                    grup = [Grouppingmedium(self.mediuminhfiels[i][j],  self.mediuminhfiels[i][j][3], self.celine[i],
                                            self.seeding[i], self.condition[i])]
                    if info1 == 2:
                        self.plotting.plot_control_mediuminhibited(grup[0], self.data_plot, self.std_plot,
                                                                   self.inhibitemediumdata[self.mediuminhfiels[i][j]
                                                                   [-1]],
                                                                   self.inhibitemediumstd[self.mediuminhfiels[i][j]
                                                                   [-1]], self.control_info, info1, info2, self.readout,
                                                                   self.readout_units)
                    elif info1 == 1:
                        self.plotting.plot_control_mediuminhibited(grup[0], self.data_plot, self.std_plot,
                                                                   self.translationmediumdata[self.mediuminhfiels[i][j]
                                                                   [-1]],
                                                                   self.inhibitemediumstd[self.mediuminhfiels[i][j]
                                                                   [-1]], self.control_info, info1, info2, self.readout,
                                                                   self.readout_units)

            except IndexError:
                pass

class SingleCompoundonecontrol(ExperimentSingleCompound):
    def __init__(self, branch, celline, seeding, fields, elapsed, control_name,
                 compound_alone, condition, std_type, medium, mediumdata, mediumstd, mediuminhfiels,
                 mediuminhib, mediuminhibstd, mediumitranslation, information_readout, readout_units):
        super().__init__(branch, celline, seeding, fields, elapsed, control_name,
                         compound_alone, condition, std_type, medium, mediumdata, mediumstd, mediuminhfiels,
                         mediuminhib, mediuminhibstd, mediumitranslation, information_readout, readout_units)

    def confluencyonecontrolovertime(self, data, std):
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
                self.plotting.plot_single_compond(grup[0], self.control_info, self.data_plot, self.std_plot, 0, 0,
                                                  self.readout, self.readout_units)

    def doseresponseonecontrol(self, data, std):
        self.data_plot = data
        self.std_plot = std

        for i in range(len(self.comalone_list_group)):
            for j in range(len(self.comalone_list_group[i])):
                self.control_info = self.control[i]
                grup = [Groupping(c, self.compound_alone[i][j], self.celine[i], self.seeding,
                                  self.condition) for c in self.comalone_list_group[i][j]]
                if len(grup[0].concentration) > 2:
                    self.get_row_plot_pag()
                    # # get x label => concentration range + lowest value for DMSO
                    self.get_concentration_range(grup)
                    self.get_confluency_range(grup)
                    DoseResponseSingle(self.data_plot, self.std_range, grup[0], self.row_num,
                                       self.concentration_range,
                                       self.confluency_range, self.time_position, self.time_selected,
                                       self.txt_path,
                                       self.pdf_pages, self.readout, self.readout_units)

    def inhibitiononecontrolovertime(self, data, std, info1, info2):
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
                self.plotting.plot_single_compond(grup[0], self.control_info, self.data_plot, self.std_plot,
                                                  info1, info2, self.readout, self.readout_units)
