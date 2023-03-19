import numpy as np

# import os
# from grupping import Groupping
from HTSplotter.grupping import Groupping, Data_group, Grouppingmedium
# import matplotlib.pyplot as plt
# from doseresponse import DoseResponseCondition
from HTSplotter.plotting import Overtime
from HTSplotter.timepointselection import Timepointselection
from HTSplotter.commonfunctions import CommonFunctions, OpenPdf, SaveTXTinfo, GrwothRateSaveTXT
from HTSplotter.growthrate import GrowthRate, GrowthRateCompoundscreenSeveral
# from matplotlib.pyplot import cm
# from matplotlib.backends.backend_pdf import PdfPages
# import math
from copy import deepcopy


class ExperimentGeneticPerturbagem:
    def __init__(self, header_info, hdfgenetic, file_info, catego,  file_names,
                 information_readout, readout_units, i, biologicalreplicate):

        # ##Usefule information for each experiment, information gave by the calls==>new_save_data_base:
        self.experiment_name = i
        self.stdinfo = catego.stdinfo
        self.elapse = file_info.elapsed
        self.date_info = file_info.date_info
        self.readout = information_readout
        self.readout_units = readout_units
        self.biologicalreplicateinfo = biologicalreplicate

        self.header = hdfgenetic.fields
        self.originaldata = hdfgenetic.data
        self.stdinh = hdfgenetic.std_inh
        self.normalized_perc = hdfgenetic.normalized_perc
        self.inhib_data = hdfgenetic.inhibited
        self.normtozero = hdfgenetic.normtozero

        self.mediuminhfiels = hdfgenetic.fieldsmedium
        self.mediumdata = hdfgenetic.datamedium
        self.mediumstd = hdfgenetic.stdmedium
        self.inhibitemediumdata = hdfgenetic.inhibitedmedium
        self.inhibitemediumstd = hdfgenetic.std_inhmedium

        self.compound_alone = hdfgenetic.compoundalone
        self.branch = header_info.branch
        self.seeding = hdfgenetic.seeding
        self.celine = hdfgenetic.celline
        self.condition = hdfgenetic.condition
        self.control = hdfgenetic.control

        self.fileioriginaldatapath = file_names.fileioriginaldatapath
        self.fileinhibiteddatapath = file_names.fileinhibiteddatapath
        self.filenormalizedatapath = file_names.filenormalizedatapath
        self.filegrowthrateresults = file_names.filegrowthrateresults

        self.std_type = None
        self.plot_titel = None
        self.time_position = None
        self.time_selected = None
        self.comalone_list_group = None
        self.data_plot = None
        self.std_plot = None
        self.perturbagem = None
        self.plotting = None
        self.concentration_range = None

        # files names to save data
        self.pdf_path = None
        self.folder = None
        self.pdf_pages = None

        # growth rate
        self.grresults = None
        self.time_selectedgr = None
        self.grheader = None
        self.position_gr = []
        self.position_elapsegr = []
        self.grselectedelap = []
        self.grselectedtimepoint = []
        self.icgrowthratetxt_path = file_names.fileicgrtxtresultspath
        if len(self.elapse) > 1:
            self.extremepoints = [self.elapse[1], self.elapse[-2]]
            self.extremepointspositions = []

        # get main time point
        timeselect = Timepointselection(self.elapse)
        self.time_selected = timeselect.time_selected
        self.time_position = timeselect.time_position

        if len(self.elapse) > 1:
            self.position_gr = []
            gr = GrowthRateCompoundscreenSeveral(self, 1)
            self.grresults = gr.grresults
            # self.extremepoints = [self.elapse[1], self.elapse[-2]]
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
        # files names to save data
        self.pdf_path = file_names.filepdfresultspath

        self.startcommonfunctions = CommonFunctions(self)

        self.row_num = self.startcommonfunctions.row_num
        self.colun_num = self.startcommonfunctions.colun_num

        self.make_empty_list_single()
        if len(self.elapse) == 1:
            self.get_headersendpoint()
        else:
            self.get_headers_single()
        SaveTXTinfo(self, self.fileioriginaldatapath, self.originaldata)
        SaveTXTinfo(self, self.filenormalizedatapath, self.normalized_perc)
        SaveTXTinfo(self, self.fileinhibiteddatapath, self.inhib_data)

        self.startpdf = OpenPdf(self)
        self.pdf_path = self.startpdf.pdf_path
        self.pdf_pages = self.startpdf.pdf_pages
        self.plotting = Overtime(self)

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

    # def open_pdf(self, pdfpathname, experiment_name):
    #     self.experiment_name = experiment_name
    #     self.pdf_path = pdfpathname
    #     self.pdf_pages = PdfPages(self.pdf_path)
    #     self.plotting = Overtime(self.pdf_pages, self.experiment_name, self.elapse, self.stdinfo)
    #     print(self.pdf_path)

    def close_pdf(self):
        self.startpdf.close_pdf()

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
                if info1 == 0:
                    self.plotting.plot_grwothrate(self, grup[0], i, j, 3, 2)
            if info3 == 1:
                self.plotting.heat_map_perturbation(grup[0], self.data_plot, self.perturbagem[i], self.control_info)

    def perturbagemendpoint(self, data, std, info1, info2):
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
            # for j in range(len(self.compound_alone[i])):
            #     self.control_info = self.control[i]
            #     grup = [Groupping(c, self.compound_alone[i][j], self.celine[i], self.seeding[i],
            #                       self.condition[i]) for c in self.comalone_list_group[i][j]]
            #     print(grup[0].name)
            #     if info3 == 1:
            #         self.plotting.heat_map_perturbation(grup[0], self.data_plot, self.perturbagem[i], self.control_info)



