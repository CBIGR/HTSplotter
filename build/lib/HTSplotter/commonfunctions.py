import numpy as np
import os
from grupping import Groupping
from txtsavedata import Savetxt
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm
from matplotlib import ticker
from matplotlib.backends.backend_pdf import PdfPages
import math
import scipy.interpolate as inter
from mpl_toolkits.mplot3d import Axes3D
from scipy.optimize import curve_fit

from plotting import Overtime

class CommonFunctions:

    def __init__(self, analysistype):

        self.analysistype = analysistype
        self.time_selected = analysistype.time_selected
        self.elapsed = None

        # commonfunctions methods
        self.row_num = None
        self.colun_num = None

        # ###OpenPdf subclass
        self.pdf_pages = None
        self.experiment_name = None
        self.pdf_pages = None

        # ###StartTXT subclass
        self.txt_path = None
        self.icgrowthratetxt_path = None
        self.experiment_name = None
        self.biologicalreplicateinfo = None

        self.date_info = None
        self.fields = None
        self.data = None

        self.possible_comb = None
        self.comb_matrix = None
        self.comb_name_per_group = None

        self.get_row_plot_pag()

    def get_row_plot_pag(self):
        self.row_num = math.ceil(len(self.time_selected) / 3)
        self.colun_num = 3

        return


class OpenPdf(CommonFunctions):
    def __init__(self, analysistype):
        super().__init__(analysistype)

        self.pdf_path = analysistype.pdf_path
        self.experiment_name = analysistype.experiment_name
        self.stdinfo = analysistype.stdinfo
        self.elapse = analysistype.elapse

        self.pdf_pages = PdfPages(self.pdf_path)
        # self.plotting = Overtime(self.pdf_pages, self.experiment_name, self.elapse, self.stdinfo)

    def close_pdf(self):
        self.pdf_pages.close()

class SaveTXTinfo(CommonFunctions):
    def __init__(self, analysistype, txt_path, data):
        super().__init__(analysistype)
        # this is to save all processed data
        self.txt_path = txt_path
        self.experiment_name = analysistype.experiment_name
        self.biologicalreplicateinfo = analysistype.biologicalreplicateinfo
        self.date_info = analysistype.date_info
        self.fields = analysistype.header
        self.data = data
        self.elapsed = analysistype.elapse

        Savetxt(self)

class SynergyScoreSaveTXT(CommonFunctions):
    def __init__(self, analysistype, name, score):
        super().__init__(analysistype)
        self.possible_comb = analysistype.possible_comb
        self.comb_matrix = analysistype.comb_matrix
        self.comb_name_per_group = analysistype.comb_name_per_group
        self.elapsed = analysistype.elapse

        bi_txt_path = name
        f = open(bi_txt_path, 'w')
        for i in range(len(self.possible_comb)):
            for k in range(len(self.possible_comb[i])):
                if self.comb_matrix[i][k][-1] > 1:  # Just for LOEWE testing
                    ela = np.asarray(["Compound=>" + str(self.possible_comb[i][k])])
                    elapsed_cabe = np.hstack((ela, self.elapsed))
                    temp_name = np.asarray(self.comb_name_per_group[i][k])
                    name_data = np.r_[np.asarray(temp_name)[None, :], np.asarray(score[i][k])]
                    new = np.c_[(elapsed_cabe[:, None], name_data)]
                    np.savetxt(f, new, delimiter="\t", fmt="%s")

        f.close()


class StartTXTic(CommonFunctions):
    def __init__(self, analysistype):
        super().__init__(analysistype)
        self.elapsed = analysistype.alltimepoints
        self.txt_path = analysistype.txt_path
        self.icgrowthratetxt_path = analysistype.icgrowthratetxt_path
        self.experiment_name = analysistype.experiment_name
        self.biologicalreplicateinfo = analysistype.biologicalreplicateinfo
        self.get_ic_txt_path()
        self.get_icgr_txt_path()

    def get_ic_txt_path(self):
        f = open(self.txt_path, 'w')
        if self.biologicalreplicateinfo == 0:
            f.write('This is 1 biological replicate. Be careful drawing any conclusions.' + '\n')
        else:
            f.write('Biological replicates.' + '\n')
        f.close()

    def get_icgr_txt_path(self):
        f = open(self.icgrowthratetxt_path, 'w')
        if self.biologicalreplicateinfo == 0:
            f.write('This is 1 biological replicate. Be careful drawing any conclusions.' + '\n')
        else:
            f.write('Biological replicates.' + '\n')
        f.close()


class AdjustData:
    def __init__(self, analysistype):
        self.normalized_perc = analysistype.normalized_perc
        self.inhib_data = analysistype.inhib_data
        self.inhibitedperc = analysistype.inhibitedperc

        for i in range(len(self.normalized_perc)):
            for j in range(len(self.normalized_perc[i])):
                if self.normalized_perc[i][j] >100:
                    self.normalized_perc[i][j] = 100
                if self.inhib_data[i][j] < 0:
                    self.inhib_data[i][j] = 0
                if self.inhibitedperc[i][j] < 0:
                    self.inhibitedperc[i][j] = 0

class GrwothRateSaveTXT(CommonFunctions):
    def __init__(self, analysistype):
        super().__init__(analysistype)
        self.headers = analysistype.grheader
        self.data = analysistype.grresults

        self.elapsed = analysistype.elapse
        self.growthpath = analysistype.filegrowthrateresults
        self.position_elapsegr = analysistype.position_elapsegr
        self.datacombination = None
        self.grselecteelasped = None
        self.grcombselecteelasped = None
        self.grheadercombination = None

        self.get_gr_size_elap()
        f = open(self.growthpath, 'w')
        for i in range(len(self.headers)):
            for k in range(len(self.headers[i])):
                ela = np.asarray("Growth Rate=>")
                elapsed_cabe = np.hstack((ela, self.elapsed[1:]))
                temp_name = ['_'.join(i[:-1]) for i in self.headers[i][k]]
                name_data = np.c_[np.asarray(temp_name)[:, None], np.asarray(self.grselecteelasped[i][k])]
                new = np.r_[elapsed_cabe[None, :], name_data]
                np.savetxt(f, new, delimiter="\t", fmt="%s")

        f.close()

    def TXTcombination(self, grresultscombination, grheadercombination):
        self.datacombination = grresultscombination
        self.grheadercombination = grheadercombination
        self.get_gr_size_elap_comb()
        f = open(self.growthpath, 'a')
        for com in range(len(self.grheadercombination)):
            for c in range(len(self.grheadercombination[com])):
                ela = np.asarray("Growth Rate Combination=>")
                elapsed_cabe = np.hstack((ela, self.elapsed[1:]))
                for s in range(len(self.grheadercombination[com][c])):
                    if s == 0:
                        temp_name = ['_'.join(i[:-1]) for i in self.grheadercombination[com][c][s]]
                        data = np.asarray(self.grcombselecteelasped[com][c][s])
                    else:
                        temp_name += ['_'.join(i[:-1]) for i in self.grheadercombination[com][c][s][1:]]
                        data = np.r_[data, np.asarray(self.grcombselecteelasped[com][c][s])[1:, :]]
                name_data = np.c_[np.asarray(temp_name)[:, None], np.asarray(data)]
                new = np.r_[elapsed_cabe[None, :], name_data]
                np.savetxt(f, new, delimiter="\t", fmt="%s")
                print('ceom')
        f.close()

    def get_gr_size_elap(self):
        self.grselecteelasped = [[] for i in range(len(self.data))]
        self.grselectedtimepoint = [[] for i in range(len(self.data))]
        for k in range(len(self.data)):
            for u in range(len(self.data[k])):
                self.grselecteelasped[k].append([])
                for i in range(len(self.data[k][u])):
                    self.grselecteelasped[k][u].append([])
                    for j in range(len(self.data[k][u][i])):
                        for time in self.position_elapsegr:
                            self.grselecteelasped[k][u][i].append(self.data[k][u][i][j][time])

    def get_gr_size_elap_comb(self):
        self.grcombselecteelasped = [[] for i in range(len(self.data))]
        for k in range(len(self.datacombination)):
            for u in range(len(self.datacombination[k])):
                self.grcombselecteelasped[k].append([])
                for i in range(len(self.datacombination[k][u])):
                    self.grcombselecteelasped[k][u].append([])
                    for j in range(len(self.datacombination[k][u][i])):
                        self.grcombselecteelasped[k][u][i].append([])
                        for time in self.position_elapsegr:
                            self.grcombselecteelasped[k][u][i][j].append(self.datacombination[k][u][i][j][0][time])


