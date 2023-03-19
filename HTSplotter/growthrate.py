# from txtsavedata import Savetxt
# import matplotlib.pyplot as plt
# from matplotlib.pyplot import cm
# from matplotlib import ticker
# from matplotlib.backends.backend_pdf import PdfPages
# import math
# import scipy.interpolate as inter
# from mpl_toolkits.mplot3d import Axes3D
# from scipy.optimize import curve_fit
import numpy as np
from HTSplotter.grupping import GrouppingGrowthRate, GrouppingGrowthRateCombination
from scipy.interpolate import UnivariateSpline, PchipInterpolator, LSQUnivariateSpline, splrep, splev, make_lsq_spline,\
make_interp_spline
# from scipy.special import erfinv
from scipy.ndimage import gaussian_filter

# from plotting import Overtime

class GrowthRate:

    def __init__(self, analysistype):

        self.elapsed = analysistype.elapse
        self.compound_alone = analysistype.compound_alone

        self.header = analysistype.header
        self.grresultscombination = None
        self.grheadercombination = None
        self.grselecteelasped = None
        self.grselectedtimepoint = None
        self.onecontrol = None

        self.data_growth = analysistype.normtozero
        self.cell_name = analysistype.celine
        self.condition = analysistype.condition
        self.seeding = analysistype.seeding
        self.taxa = []

        print('sto')

    @staticmethod
    def cheb_nodes(N):
        jj = 2. * np.arange(N) + 1
        x = np.cos(np.pi * jj / 2 / N)[::-1]
        return x

    def get_derivatives(self, positions):
        auxlist2 = []
        for conce in range(len(positions) + 1):
            xs = np.linspace(self.elapsed[0], self.elapsed[-1], int(1e5))
            if conce == 0:

                aux_x = self.elapsed.copy()

                ys = self.data_growth[self.control_info[-1]]
                ys = gaussian_filter(ys, 2.5)

                nodes = np.linspace(self.elapsed[0], self.elapsed[-1], int(round(self.elapsed.size*0.1)))
                nodes = np.r_[(nodes[0],) * (4+1), nodes, (nodes[-1],) * (4+1)]
                spl = make_lsq_spline(aux_x, ys, k=5, t=nodes)

                ys_aux = spl(xs)
                ysc_aux = ys_aux.copy()

            else:
                aux_x = self.elapsed.copy()

                ys = self.data_growth[positions[conce - 1]]
                ys = gaussian_filter(ys, 2.5)

                nodes = np.linspace(self.elapsed[0], self.elapsed[-1], int(round(self.elapsed.size * 0.1)))
                nodes = np.r_[(nodes[0],)*(4+1), nodes, (nodes[-1],)*(4+1)]
                spl = make_lsq_spline(aux_x, ys, k=5, t=nodes)

                ys_aux = spl(xs)

            condi = np.log2(ys/ys[0])

            if conce == 0:
                contro = condi.copy()

            auxlist = []

            gr_mean = 2**(condi/contro) - 1
            mask = ~ np.isnan(gr_mean)
            aux_x = self.elapsed[mask]
            gr_mean = gr_mean[mask]
            nodes = np.linspace(aux_x[0], aux_x[-1], int(round(aux_x.size * 0.1)))
            nodes = np.r_[(nodes[0],)*(4+1), nodes, (nodes[-1],)*(4+1)]
            spl = make_lsq_spline(aux_x, gr_mean, k=5, t=nodes)
            gr_mean = spl(xs)

            dys_dt = (ys_aux[2:] - ys_aux[:-2])/(xs[2:] - xs[:-2])
            dysc_dt = (ysc_aux[2:] - ysc_aux[:-2])/(xs[2:] - xs[:-2])

            df_dt = (1/np.log(2))*((ys_aux[0]/ys_aux[1:-1])*dys_dt*np.log2(ysc_aux[1:-1]/ysc_aux[0]) -
                                   (ysc_aux[0]/ysc_aux[1:-1])*dysc_dt*np.log2(ys_aux[1:-1]/ys_aux[0]))/\
                    (np.log2(ysc_aux[1:-1]/ysc_aux[0]))**2

            dgr_mean = gr_mean[1:-1]*np.log(2)*df_dt

            # gr = (dgr_mean[1:] + dgr_mean[:-1])*(xs[3:] - xs[:-3]/2) + (gr_mean[3:] + gr_mean[:-3])/2
            gr = (dgr_mean[1:] + dgr_mean[:-1])/2 * (xs[1:-2]) + (gr_mean[3:] + gr_mean[:-3]) / 2

            # auxlist.append(2 ** (condi / contro) - 1)
            auxlist.append(gr)
            # auxlist.append(condi)
            # auxlist.append((ys[2:] - ys[:-2])/(xs[2:] - xs[:-2]))
            # auxlist.append(ys[1:-1])
            pass
            # for tim in range(len(self.elapsed)):
            #     if conce == 0:
            #         # control data
            #
            #         if tim > 0 and tim < self.elapsed.size - 1:
            #             print('ied', tim, self.elapsed.size)
            #             condi = np.log2(self.data_growth[self.control_info[-1]][tim + 1]) - \
            #                      np.log2(self.data_growth[self.control_info[-1]][tim - 1])
            #             tem = self.elapsed[tim + 1] - self.elapsed[tim - 1]
            #             contro = np.log2(self.data_growth[self.control_info[-1]][tim + 1]) - \
            #                      np.log2(self.data_growth[self.control_info[-1]][tim - 1])
            #             auxlist.append(2 ** ((condi / tem) / (contro / tem)) - 1)
            #             print('oewf')
            #     # if tim == 0:
            #     else:
            #
            #         if tim > 0 and tim < self.elapsed.size-1:
            #             print('ied', tim, self.elapsed.size, positions[conce-1])
            #             condi = np.log2(self.data_growth[positions[conce-1]][tim + 1]) - \
            #                     np.log2(self.data_growth[positions[conce-1]][tim - 1])
            #             tem = self.elapsed[tim+1] - self.elapsed[tim-1]
            #             contro = np.log2(self.data_growth[self.control_info[-1]][tim + 1]) - \
            #                      np.log2(self.data_growth[self.control_info[-1]][tim - 1])
            #             auxlist.append(2**((condi/tem)/(contro/tem))-1)
            #             print('oewf')

            # if tim == len(self.elapsed):
            auxlist2.append(auxlist)

        return xs[3:-3], auxlist2

    def get_headers(self, grup, positions, i, j):
        for conce in range(len(positions) + 1):
            if conce == 0:
                if self.onecontrol != 0:
                    headerlist = [self.control[i].copy()]
                else:
                    headerlist = [self.control[i][j].copy()]

                # posi = self.control_info[-1]
            else:
                headerlist += [fi for fi in grup.headergroupping if fi[-1] == positions[conce - 1]]
                # posi = grup.headergroupping[positions[conce-1][-1]]
        return headerlist


class GrowthRateCompoundscreenSeveral(GrowthRate):
    def __init__(self, analysistype, onecontrol=0):
        super().__init__(analysistype)

        self.data = None
        self.data_control = None
        self.group = None
        self.control = analysistype.control
        self.time_selected = analysistype.time_selected
        self.position_gr = []
        self.grresults = [[] for i in range(len(self.compound_alone))]
        self.grresultsderivadascemtrais = [[] for i in range(len(self.compound_alone))]
        self.grheader = [[] for i in range(len(self.compound_alone))]
        self.onecontrol = onecontrol
        # self.derivadascentrais()
        for i in range(len(self.compound_alone)):
            # cell line
            # # plot all compounds alone
            for j in range(len(self.compound_alone[i])):
                if self.onecontrol != 0:
                    self.control_info = self.control[i]
                else:
                    self.control_info = self.control[i][j]
                print("is there a control", self.control_info)

                grup = GrouppingGrowthRate(self.compound_alone[i][j], self.header, self.cell_name[i],
                                           self.condition[i][0], self.seeding[i])

                positions = [int(c) for c in grup.concentration[:, -1]]

                headerlist = self.get_headers(grup, positions, i, j)

                xs, auxlist2 = self.get_derivatives(positions)

                self.xs = xs
                self.grresults[i].append(auxlist2)
                self.grheader[i].append(headerlist)

                # self.data = datagroupes.data_grup

        # for h in range(len())

class GrowthRateGeneticChemical(GrowthRate):
    def __init__(self, analysistype, onecontrol=0):
        super().__init__(analysistype)

        self.data = None
        self.data_control = None
        self.group = None
        self.com_list_group = analysistype.possible_comb
        self.simplecomb_list_group = analysistype.comb_list_group
        # self.comb_list_group
        self.conditionorder= analysistype.conditionorder
        self.control = analysistype.control
        self.time_selected = analysistype.time_selected
        self.branch = analysistype.branch
        self.position_gr = []
        self.grresults = [[] for i in range(len(self.compound_alone))]
        self.grresultsderivadascemtrais = [[] for i in range(len(self.compound_alone))]
        self.grresultscombination = [[] for i in range(len(self.com_list_group))]
        self.grheadercombination = [[] for i in range(len(self.com_list_group))]
        self.grheader = [[] for i in range(len(self.compound_alone))]
        self.onecontrol = onecontrol
        # self.derivadascentrais()
        for i in range(len(self.compound_alone)):
            # cell line
            # # plot all compounds alone

            for j in range(len(self.compound_alone[i])):
                for k in range(len(self.compound_alone[i][j])):
                    for h in range(len(self.control[i])):
                        compounname = self.control[i][h][-3].split('_')[-1]
                        if compounname in self.compound_alone[i][j][k]:
                            self.control_info = self.control[i][h]
                            print('confirm', self.control_info, self.compound_alone[i][j][k])
                            grup = GrouppingGrowthRate(self.compound_alone[i][j][k], self.header,
                                                       self.cell_name[i], self.condition[i][j],
                                                       self.seeding[i])
                            positions = [int(c) for c in grup.concentration[:, -1]]
                            # print('hey')

                            headerlist = self.get_headers(grup, positions, i, h)

                            xs, auxlist2 = self.get_derivatives(positions)


                            self.xs = xs
                            self.grresults[i].append(auxlist2)
                            self.grheader[i].append(headerlist)

                # self.data = datagroupes.data_grup

        for i in range(len(self.com_list_group)):
            # cell line
            # # plot all compounds alone
            self.grresultscombination[i] = [[] for k in range(len(self.com_list_group[i]))]
            self.grheadercombination[i] = [[] for k in range(len(self.com_list_group[i]))]
            for j in range(len(self.com_list_group[i])):
                for h in range(len(self.control[i])):
                    compounname = self.control[i][h][-3].split('_')[-1]
                    if compounname in self.com_list_group[i][j]:
                        self.control_info = self.control[i][h]
                        print('confirm', self.control_info, self.com_list_group[i][j])
                        grup = [GrouppingGrowthRateCombination(c, self.header,
                                                               self.cell_name[i],
                                                               self.seeding[i])for c in self.com_list_group[i][j]]
                        print('fom')
                        coun = 0
                        for grp in grup:
                            positions = [int(c) for c in grp.concentration[:, -1]]

                            headerlist = self.get_headers(grp, positions, i, h)

                            xs, auxlist2 = self.get_derivatives(positions)

                            self.xs = xs
                            # if coun == 0:
                            self.grresultscombination[i][j].append(auxlist2)
                            self.grheadercombination[i][j].append(headerlist)

                            coun +=1
                        print('to cehc')
    # def derivadascentrais(self):
    #     for i in range(len(self.compound_alone)):
    #         # cell line
    #         # # plot all compounds alone
    #         # self.grheader[i] = []
    #         for j in range(len(self.compound_alone[i])):
    #             self.control_info = self.control[i][j]
    #             print(self.control_info)
    #             grup = GrouppingGrowthRate(self.compound_alone[i][j], self.header, self.cell_name[i],
    #                                        self.condition[i][0], self.seeding[i])
    #
    #             positions = [int(c) for c in grup.concentration[:, -1]]
    #
    #             auxlist2 = []
    #             for conce in range(len(positions)+1):
    #
    #                 # if conce == 0:
    #                 #     spl = UnivariateSpline(self.elapsed, self.data_growth[self.control_info[-1]],
    #                 #                            s=self.elapsed.size*5.0, check_finite=True)
    #                 #
    #                 # else:
    #                 #     spl = UnivariateSpline(self.elapsed, self.data_growth[positions[conce-1]],
    #                 #                            s=self.elapsed.size*5.0, check_finite=True)
    #                 #
    #                 # xs = np.linspace(self.elapsed[0], self.elapsed[-1], int(1e5))
    #                 # ys = spl(xs)
    #                 #
    #                 # condi = ((np.log2(ys[2:]) - np.log2(ys[:-2]))/(xs[2:] - xs[:-2])).flatten()
    #                 # if conce == 0:
    #                 #     contro = condi.copy()
    #
    #                 auxlist = []
    #                 # if conce == 0:
    #                 #     headerlist = [self.control[i][j].copy()]
    #                 #     # posi = self.control_info[-1]
    #                 # else:
    #                 #     headerlist += [fi for fi in grup.headergroupping if fi[-1] == positions[conce-1]]
    #                     # posi = grup.headergroupping[positions[conce-1][-1]]
    #
    #                 # auxlist.append(2 ** (contro - condi) - 1)
    #                 # auxlist.append(ys[1:-1])
    #                 for tim in range(len(self.elapsed)):
    #                     if conce == 0:
    #                         # control data
    #
    #                         if tim > 0 and tim < self.elapsed.size - 1:
    #                             print('ied', tim, self.elapsed.size)
    #                             condi = np.log2(self.data_growth[self.control_info[-1]][tim + 1]) - \
    #                                      np.log2(self.data_growth[self.control_info[-1]][tim - 1])
    #                             tem = self.elapsed[tim + 1] - self.elapsed[tim - 1]
    #                             contro = np.log2(self.data_growth[self.control_info[-1]][tim + 1]) - \
    #                                      np.log2(self.data_growth[self.control_info[-1]][tim - 1])
    #                             auxlist.append(2 ** ((condi / tem) / (contro / tem)) - 1)
    #                             print('oewf')
    #                     # if tim == 0:
    #                     else:
    #                         if tim > 0 and tim < self.elapsed.size-1:
    #                             print('ied', tim, self.elapsed.size, positions[conce-1])
    #                             condi = np.log2(self.data_growth[positions[conce-1]][tim + 1]) - \
    #                                     np.log2(self.data_growth[positions[conce-1]][tim - 1])
    #                             tem = self.elapsed[tim+1] - self.elapsed[tim-1]
    #                             contro = np.log2(self.data_growth[self.control_info[-1]][tim + 1]) - \
    #                                      np.log2(self.data_growth[self.control_info[-1]][tim - 1])
    #                             auxlist.append(2**((condi/tem)/(contro/tem))-1)
    #                             print('oewf')
    #
    #                     # if tim == len(self.elapsed):
    #                 auxlist2.append(np.asarray(auxlist))
    #
    #             self.grresultsderivadascemtrais[i].append(auxlist2)
    #             # self.grheader[i].append(headerlist)
    #
    #             # self.data = datagroupes.data_grup
    # def growthrate(self):
    #     for i in range(len(self.compound_alone)):
    #         # cell line
    #         # # plot all compounds alone
    #         # self.grheader[i] = []
    #         for j in range(len(self.compound_alone[i])):
    #             conditionposition = 0
    #             if onecontrol != 0:
    #                 self.control_info = self.control[i]
    #             else:
    #                 self.control_info = self.control[i][j]
    #             print(self.control_info)
    #             # if len(self.condition[i]) >1:
    #             #     grup = GrouppingGrowthRate(self.compound_alone[i][j], self.header, self.cell_name[i],
    #             #                                self.condition[i][conditionposition], self.seeding[i])
    #             #     conditionposition +=1
    #             # else:
    #             grup = GrouppingGrowthRate(self.compound_alone[i][j], self.header, self.cell_name[i],
    #                                        self.condition[i], self.seeding[i])
    #
    #             positions = [int(c) for c in grup.concentration[:, -1]]
    #
    #             auxlist2 = []
    #             for conce in range(len(positions) + 1):
    #                 va = 3
    #                 xs = np.linspace(self.elapsed[0], self.elapsed[-1], int(1e5))
    #                 if conce == 0:
    #
    #                     aux_x = self.elapsed.copy()
    #
    #                     ys = self.data_growth[self.control_info[-1]]
    #                     # m_end = np.mean((ys[-1] - ys[-va:-1])/(self.elapsed[-1] - self.elapsed[-va:-1]))
    #                     # m_ini = np.mean((ys[1:va] - ys[0])/(self.elapsed[1:va] - self.elapsed[0]))
    #                     # dx = self.elapsed[-1] - self.elapsed[-2]
    #                     # y_end = m_end * np.asarray([a*dx for a in range(va-1)]) + ys[-1]
    #                     # dx = self.elapsed[1] - self.elapsed[0]
    #                     # y_ini = m_ini * np.asarray([-a * dx for a in range(va-1, 0, -1)]) + ys[0]
    #                     #
    #                     # aux_x = np.r_[np.asarray([-a * dx for a in range(va-1, 0, -1)]),
    #                     #               aux_x, np.asarray([a*dx + aux_x[-1] for a in range(1, va)])]
    #
    #                     # ys = gaussian_filter(ys, 2.5)
    #                     # ys = np.r_[y_ini, ys, y_end]
    #
    #                     # spl = PchipInterpolator(aux_x, ys)  # ,
    #                     # s=self.elapsed.size*0.0, check_finite=True)
    #                     # ys = spl(aux_x)
    #                     spl = splrep(aux_x, ys, k=5)
    #                     ys = splev(xs, spl)
    #                     # ys = gaussian_filter(ys, 2.5)
    #                     # ys = ys[va:-va]
    #                     # spl = splrep(aux_x, ys)
    #                     # ys = splev(xs, spl)
    #
    #                 else:
    #
    #                     aux_x = self.elapsed.copy()
    #
    #                     ys = self.data_growth[positions[conce - 1]]
    #
    #                     m_end = np.mean((ys[-1] - ys[-va:-1]) / (self.elapsed[-1] - self.elapsed[-va:-1]))
    #                     m_ini = np.mean((ys[1:va] - ys[0]) / (self.elapsed[1:va] - self.elapsed[0]))
    #                     dx = self.elapsed[-1] - self.elapsed[-2]
    #                     y_end = m_end * np.asarray([a * dx for a in range(va - 1)]) + ys[-1]
    #                     dx = self.elapsed[1] - self.elapsed[0]
    #                     y_ini = m_ini * np.asarray([-a * dx for a in range(va - 1, 0, -1)]) + ys[0]
    #
    #                     ys = gaussian_filter(ys, 2.5)
    #                     ys = np.r_[y_ini, ys, y_end]
    #
    #                     aux_x = np.r_[np.asarray([-a * dx for a in range(va - 1, 0, -1)]),
    #                                   aux_x, np.asarray([a * dx + aux_x[-1] for a in range(1, va)])]
    #
    #                     ys = gaussian_filter(ys, 2.5)
    #
    #                     for w in range(1):
    #                         spl = splrep(aux_x, ys, k=5)
    #                         if w < -1:
    #                             aux_x = np.r_[np.asarray([-a * dx for a in range(va - 1, 0, -1)]),
    #                                           np.random.choice(aux_x, len(self.elapsed)),
    #                                           np.asarray([a * dx + aux_x[-1] for a in range(1, va)])]
    #                             aux_x = np.unique(aux_x)
    #                             ys = splev(np.unique(aux_x), spl)
    #
    #                     # spl = UnivariateSpline(self.elapsed, self.data_growth[positions[conce-1]],
    #                     #                        s=self.elapsed.size*5.0, check_finite=True)
    #
    #                     ys = splev(xs, spl)
    #                 # ys = gaussian_filter(ys, 5000)
    #                 # condi = ((np.log2(ys) - np.log2(ys)) / (xs - xs)).flatten()
    #                 condi = ((np.log2(ys[2:]) - np.log2(ys[:-2])) / (xs[2:] - xs[:-2])).flatten()
    #                 if conce == 0:
    #                     condi = gaussian_filter(condi, 10)
    #                     contro = condi.copy()
    #
    #                 auxlist = []
    #                 if conce == 0:
    #                     if onecontrol != 0:
    #                         headerlist = [self.control[i].copy()]
    #                     else:
    #                         headerlist = [self.control[i][j].copy()]
    #                     # posi = self.control_info[-1]
    #                 else:
    #                     headerlist += [fi for fi in grup.headergroupping if fi[-1] == positions[conce - 1]]
    #                     # posi = grup.headergroupping[positions[conce-1][-1]]
    #
    #                 auxlist.append(2 ** (condi / contro) - 1)
    #                 # auxlist.append(condi)
    #                 # auxlist.append(ys[1:-1])
    #                 pass
    #                 # for tim in range(len(self.elapsed)):
    #                 #     if conce == 0:
    #                 #         # control data
    #                 #
    #                 #         if tim > 0 and tim < self.elapsed.size - 1:
    #                 #             print('ied', tim, self.elapsed.size)
    #                 #             condi = np.log2(self.data_growth[self.control_info[-1]][tim + 1]) - \
    #                 #                      np.log2(self.data_growth[self.control_info[-1]][tim - 1])
    #                 #             tem = self.elapsed[tim + 1] - self.elapsed[tim - 1]
    #                 #             contro = np.log2(self.data_growth[self.control_info[-1]][tim + 1]) - \
    #                 #                      np.log2(self.data_growth[self.control_info[-1]][tim - 1])
    #                 #             auxlist.append(2 ** ((condi / tem) / (contro / tem)) - 1)
    #                 #             print('oewf')
    #                 #     # if tim == 0:
    #                 #     else:
    #                 #
    #                 #         if tim > 0 and tim < self.elapsed.size-1:
    #                 #             print('ied', tim, self.elapsed.size, positions[conce-1])
    #                 #             condi = np.log2(self.data_growth[positions[conce-1]][tim + 1]) - \
    #                 #                     np.log2(self.data_growth[positions[conce-1]][tim - 1])
    #                 #             tem = self.elapsed[tim+1] - self.elapsed[tim-1]
    #                 #             contro = np.log2(self.data_growth[self.control_info[-1]][tim + 1]) - \
    #                 #                      np.log2(self.data_growth[self.control_info[-1]][tim - 1])
    #                 #             auxlist.append(2**((condi/tem)/(contro/tem))-1)
    #                 #             print('oewf')
    #
    #                 # if tim == len(self.elapsed):
    #                 auxlist2.append(auxlist)
    #             self.xs = xs
    #             self.grresults[i].append(auxlist2)
    #             self.grheader[i].append(headerlist)
    #
    #             # self.data = datagroupes.data_grup

class GrowthRateCombination(GrowthRate):
    def __init__(self, analysistype, onecontrol=0):
        super().__init__(analysistype)

        self.data = None
        self.data_control = None
        self.group = None
        self.com_list_group = analysistype.possible_comb
        self.control = analysistype.control
        self.time_selected = analysistype.time_selected
        self.branch = analysistype.branch
        self.position_gr = []
        self.grresults = [[] for i in range(len(self.compound_alone))]
        self.grresultsderivadascemtrais = [[] for i in range(len(self.compound_alone))]
        self.grresultscombination = [[] for i in range(len(self.com_list_group))]
        self.grheadercombination = [[] for i in range(len(self.com_list_group))]
        self.grheader = [[] for i in range(len(self.compound_alone))]
        self.onecontrol = onecontrol
        # self.derivadascentrais()
        for i in range(len(self.compound_alone)):
            # cell line
            # # plot all compounds alone
            for j in range(len(self.compound_alone[i])):
                    self.control_info = self.control[i]
                    print('confirm', self.control_info, self.compound_alone[i][j])
                    grup = GrouppingGrowthRate(self.compound_alone[i][j], self.header, self.cell_name[i],
                                               self.condition[i][0], self.seeding[i])
                    positions = [int(c) for c in grup.concentration[:, -1]]
                    # print('hey')

                    headerlist = self.get_headers(grup, positions, i, j)

                    xs, auxlist2 = self.get_derivatives(positions)


                    self.xs = xs
                    self.grresults[i].append(auxlist2)
                    self.grheader[i].append(headerlist)

                # self.data = datagroupes.data_grup

        for i in range(len(self.com_list_group)):
            # cell line
            # # plot all compounds alone
            self.grresultscombination[i] = [[] for k in range(len(self.com_list_group[i]))]
            self.grheadercombination[i] = [[] for k in range(len(self.com_list_group[i]))]
            for j in range(len(self.com_list_group[i])):
                self.control_info = self.control[i]
                print('confirm', self.control_info, self.com_list_group[i])
                grup = [GrouppingGrowthRateCombination(c, self.header,
                                                       self.cell_name[i],
                                                       self.seeding[i])for c in self.com_list_group[i][j]]
                print('fom')
                coun = 0
                for grp in grup:
                    positions = [int(c) for c in grp.concentration[:, -1]]

                    headerlist = self.get_headers(grp, positions, i, j)

                    xs, auxlist2 = self.get_derivatives(positions)

                    self.xs = xs
                    # if coun == 0:
                    self.grresultscombination[i][j].append(auxlist2)
                    self.grheadercombination[i][j].append(headerlist)
                    coun +=1
                    print('to cehc')
