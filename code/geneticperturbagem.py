import numpy as np

import os
from grupping import Groupping, Data_group, Grouppingmedium
import matplotlib.pyplot as plt
from doseresponse import DoseResponseSingle, DoseResponseCondition
from plotting import Overtime
from timepointselection import Timepointselection

from matplotlib.pyplot import cm
from matplotlib.backends.backend_pdf import PdfPages
import math
from copy import deepcopy


class ExperimentGeneticPerturbagem:
    def __init__(self, branch, celline, seeding, fields, elapsed, control_name,
                 compound_alone, condition, std_type, medium, mediumdata, mediumstd, mediuminhfiels, mediuminhib,
                 mediuminhibstd, information_readout, readout_units):
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
        self.perturbagem = None
        self.plotting = None

        self.concentration_range = None

        # get main time point
        timeselect = Timepointselection(self.elapse)
        self.time_selected = timeselect.time_selected
        self.time_position = timeselect.time_position

        self.make_empty_list_single()
        if len(self.elapse) == 1:
            self.get_headersendpoint()
        else:
            self.get_headers_single()

    def make_empty_list_single(self):

        list1 = []
        for j in range(len(self.compound_alone)):
            header_count = len(self.compound_alone[j])
            headers = [[] for i in range(0, header_count)]  # list_possible_combination
            list1.append(headers)
        self.comalone_list_group = deepcopy(list1)
        self.perturbagem = deepcopy(list1)

    def get_headers_single(self):
        for i in range(len(self.compound_alone)):
            for k in range(len(self.compound_alone[i])):
                list1 = []
                for u in self.header:
                    if u[0] == self.celine[i] and u[3] == self.compound_alone[i][k]:
                        list1.append(u)
                        self.perturbagem[i][k].append(u)
                self.comalone_list_group[i][k].append(list1)

    def get_headersendpoint(self):
        for i in range(len(self.compound_alone)):
            for k in range(len(self.compound_alone[i])):
                list1 = []
                for u in self.header:
                    if u[0] == self.celine[i] and u[3] == self.compound_alone[i][k]:
                        list1.append(u)
                        self.comalone_list_group[i][k].append(u)

    def open_pdf(self, pdfpathname, experiment_name):
        self.experiment_name = experiment_name
        self.pdf_path = pdfpathname
        self.pdf_pages = PdfPages(self.pdf_path)
        self.plotting = Overtime(self.pdf_pages, self.experiment_name, self.elapse, self.stdinfo)


    def close_pdf(self):
        self.pdf_pages.close()

    def perturbagemovertime(self, data, std, info1, info2, info3):
        self.data_plot = data
        self.std_plot = std
        for i in range(len(self.compound_alone)):
            # cell line
            # # plot all compounds alone
            for j in range(len(self.compound_alone[i])):
                self.control_info = self.control[i]
                grup = [Groupping(c, self.compound_alone[i][j], self.celine[i], self.seeding[i],
                                  self.condition[i]) for c in self.comalone_list_group[i][j]]
                self.plotting.plot_perturbagem(grup[0], self.control_info, self.data_plot, self.std_plot, info1,
                                               info2, info3, self.time_selected, self.time_position, self.readout,
                                               self.readout_units)
            if info3 == 1:
                self.plotting.heat_map_perturbation(grup[0], self.data_plot, self.perturbagem[i], self.control_info)

    def perturbagemendpoint(self, data, std, info1, info2, info3):
        self.data_plot = data
        self.std_plot = std
        for i in range(len(self.compound_alone)):
            # cell line
            # # plot all compounds alone
            grup = [Groupping(c, self.compound_alone[i], self.celine[i], self.seeding[i],
                              self.condition[i]) for c in self.comalone_list_group[i]]

            self.plotting.barplot_geneticperturbgen_recursive(grup, self.data_plot, self.std_plot, info1, info2,
                                                              self.time_selected, self.time_position, self.readout,
                                                              self.readout_units)


