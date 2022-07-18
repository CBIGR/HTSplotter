import numpy as np
import os
# from grupping import Groupping, Data_group
# from synergism import Blissmethod
os.environ['MPLCONFIGDIR'] = '/opt/HTSplotter/.config/matplotlib'
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm
from arrayRGB import ColorsHtsplots
import random
# from matplotlib.backends.backend_pdf import PdfPages
import math
from interactionsuser import Outputfile
import scipy.interpolate as inter
from scipy.interpolate import UnivariateSpline
from interactionsuser import Outputfile


# from doseresponse import DoseResponseSingle
# import time

# import time as ti


class Overtime:

    def __init__(self, analysis, synergy_method=False):
        # self.pdf_pages, self.experiment_name, self.elapse, self.stdinfo

        self.plot_titel = analysis.experiment_name
        self.elapse = analysis.elapse
        self.std_type = analysis.stdinfo[0]
        self.pdf_pages = analysis.pdf_pages
        self.time_selected = analysis.time_selected
        self.time_position = analysis.time_position
        self.time_selectedgr = analysis.time_selectedgr
        self.position_gr = analysis.position_gr
        # self.position_elapsegr = analysis.position_elapsegr
        self.readout = analysis.readout
        self.readout_units = analysis.readout_units

        self.datanormtozero = analysis.normtozero
        self.data = None
        self.grup = None
        self.control_info = None
        self.data_plot = None
        self.std_plot = None
        self.control_list = None

        self.plotdata = None
        self.header = None
        self.extremepointspositions = None
        self.elapsegrwothrate = None

        self.legtitel = None
        self.mediumdata = None
        self.mediumstd = None
        self.mediuminhfields = None
        self.colormap = None
        self.controldata = None

        #  atributes for combination
        self.comb_name_per_group = None
        self.effect_per_group = None
        self.synergy_score_per_group = None
        self.synergy_method = synergy_method

        # atributes for Growth Rate
        self.grselecteelasped = None
        self.grselectedtimepoint = None

        if synergy_method is False:
            self.synergy_text = ' BI score'
        else:
            self.synergy_method = analysis.synergy_method
            if self.synergy_method == 0:
                self.synergy_text = ' BI score'
            elif self.synergy_method == 1:
                self.synergy_text = ' HSA score'
            # elif self.synergy_method == 2:
            #     self.synergy_text = ' LOEWE method--> CI (combination index) score'
            elif self.synergy_method == 2:
                self.synergy_text = ' ZIP score'


    def set_y_axis_name(self, inf):
        y_axis_title = None
        if inf == 0:
            y_axis_title = self.readout + " " + self.readout_units
        if inf == 1:
            y_axis_title = "relative " + self.readout
        if inf == 2:
            y_axis_title = "relative " + self.readout + " inhbition"
        if inf == 3:
            y_axis_title = "Growth Rate "

        return y_axis_title

    def set_y_axis_limt(self, inf, maxim, minin):

        y_axis_min = None
        y_axis_max = None
        if minin == None:
            y_axis_min = 0
        # stepsize = round((maximoriginal-minim)/len(self.data_plot[0]), 1)
        if inf == 0:
            y_axis_min = 0
            if maxim > 10 and maxim < 100:
                y_axis_max = 120
            else:
                y_axis_max = maxim * 1.5  # 120

        if inf == 1:
            y_axis_min = minin*0.5
            if minin > 0:
                y_axis_min = -1.5
            if minin == 0:
                y_axis_min = -1.5
            if minin < 0 and minin > -1.5:
                y_axis_min = -1.5
            if maxim > 1.5:
                y_axis_max = maxim
            if maxim < 1.5:
                y_axis_max = 1.5
        # if inf == 2:
        #     y_axis_min = minin
        #     y_axis_max = maxim
        if inf == 3:
            if minin < -1:
                y_axis_min = minin
            else:
                y_axis_min = -1
            if maxim > 2:
                y_axis_max = maxim
            else:
                y_axis_max = 2
        if inf == 4:
            y_axis_min = 0
            y_axis_max = maxim+1

        return y_axis_min, y_axis_max

    def get_leg_title(self, grup):

        title_leg = grup.cell_name + " " + grup.seeding.replace("per", "/") + "\n" + self.std_type

        return title_leg

    def plot_info(self, fr, grup, counter, ind=0):

        if len(grup) > ind:
            if counter is None:
                counter = np.zeros(len(grup), dtype=np.int16)

            group = grup[ind]  # .pop(0)

            for i in range(ind, -1, -1):
                if counter[i] == grup[i].concentration.shape[0]:
                    counter[i] = 0
                    if i > 0:
                        counter[i - 1] += 1

            tem = str(group.name) + " " + str(group.concentration[counter[ind], 0]) + " " + str(group.unidade)
            fr.append(tem)
            fr, counter, counter_old = self.plot_info(fr, grup, counter, ind + 1)

            if ind == len(counter) - 1:
                counter_old = counter.copy()
                counter[ind] += 1

            return fr, counter, counter_old

        else:
            counter_old = counter

            return fr, counter, counter_old

    def plot_infogr(self, fr, grup, counter, ind=0):

        if len(grup) > ind:
            if counter is None:
                counter = np.zeros(len(grup), dtype=np.int16)

            group = grup[ind]  # .pop(0)

            for i in range(ind, -1, -1):
                if counter[i] == grup[i].concentration.shape[0]:
                    counter[i] = 0
                    if i > 0:
                        counter[i - 1] += 1

            tem = str(group.name) + " " + str(group.concentration[counter[ind], 0]) + " " + str(group.unidade)
            fr.append(tem)
            fr, counter, counter_old = self.plot_infogr(fr, grup, counter, ind + 1)

            if ind == len(counter) - 1:
                counter_old = counter.copy()
                counter[ind] += 1

            return fr, counter, counter_old

        else:
            counter_old = counter

            return fr, counter, counter_old

    def get_gr(self):
        self.grselecteelasped = [[] for i in range(len(self.plotdata))]
        self.grselectedtimepoint = [[] for i in range(len(self.plotdata))]
        for i in range(len(self.plotdata)):
            for j in range(len(self.plotdata[i])):
                self.grselecteelasped[i].append([])
                for time in self.position_gr:
                    self.grselecteelasped[i][j].append(self.plotdata[i][j][0][time])


    @staticmethod
    def get_addition_text(titel, biscore, timeselected, timeposition, cout):

        for a in range(len(timeselected)):
            posit = timeposition[a]
            titel.append(str(timeselected[a]) + " h=" +
                         str("{:.2f}".format(biscore[posit, cout]))
                         + " ")

        te = " ".join(titel)

        return te

    def get_addition_grtext(self, titel, textnamegr, position):
        text = []
        for condi in range(len(textnamegr)):
            sit = position[condi]
            for a in range(len(self.grselecteelasped)):
                # posit = self.time_selected[a]
                if a == 0:
                    text.append(textnamegr[condi] + ' ' +str(self.time_selected[a]) + " h=" +
                                 str("{:.2f}".format(self.grselecteelasped[sit[0]][sit[-1]][a]))
                                 + " ")
                else:
                    text[-1] += str(self.time_selected[a]) + " h=" + \
                        str("{:.2f}".format(self.grselecteelasped[sit[0]][sit[-1]][a])) + str(' ')
            # text = " ".join(text)
        te = "\n".join(text)

        return te

    @staticmethod
    def get_perturbagem_text(titel, biscore, timeselected, timeposition):

        for a in range(len(timeselected)):
            posit = timeposition[a]
            titel.append(str(timeselected[a]) + " h=" +
                         str("{:.2f}".format(biscore[posit]))
                         + " ")

        te = " ".join(titel)

        return te

    # XY plots overtime
    def plot_single_compond(self, grup, forplotting, info1=0, info2=0):

        self.control_info = forplotting.control_info
        self.grup = grup
        self.data_plot = forplotting.data_plot
        self.std_plot = forplotting.std_plot
        self.readout = forplotting.readout
        self.readout_units = forplotting.readout_units

        y_axis_name = self.set_y_axis_name(info1)
        y_axis_min, y_axis_max = self.set_y_axis_limt(info2, np.max(self.data_plot), np.min(self.data_plot))
        fig = plt.figure(figsize=(8, 8), dpi=80)
        ax = plt.axes()
        label_titel = "compound: " + self.grup.name
        label_b = str(self.control_info[3]) + " " + self.control_info[-2]
        titel_page = str(label_titel)
        post_control = self.control_info[-1]
        plt.plot(self.elapse, self.data_plot[post_control], label=label_b, color="black")
        plt.fill_between(self.elapse, self.data_plot[post_control] - self.std_plot[post_control],
                         self.data_plot[post_control] + self.std_plot[post_control],
                         alpha=0.5, edgecolor='lightgray', facecolor='lightgray')

        self.colormap = list(iter(cm.tab20b(np.linspace(0, 1, len(grup.concentration)))))
        cout = 0
        if len(self.colormap) < len(grup.concentration):
            col = ColorsHtsplots(len(grup.concentration))
            self.colormap = col.map

        for i in range(len(grup.concentration)):
            post1 = int(grup.concentration[i][-1])
            if len(grup.concentration[i]) == 1:
                label_a = str(grup.name)
            else:
                label_a = str(grup.name) + " " + str(grup.concentration[i][0]) + " " + grup.unidade

            plt.plot(self.elapse, self.data_plot[post1], label=label_a, color=self.colormap[i])
            plt.fill_between(self.elapse, self.data_plot[post1] - self.std_plot[post1],
                             self.data_plot[post1] + self.std_plot[post1],
                             alpha=0.5, edgecolor="lightgray", facecolor='lightgray')
            cout +=1
        plt.xlabel("time (h)")
        plt.ylabel(y_axis_name)
        plt.ylim(y_axis_min, y_axis_max)
        # plt.yticks(np.arange(y_axis_min, y_axis_max, stepsize))
        sub = plt.suptitle(self.plot_titel, fontsize=12)
        plt.title(titel_page, fontsize=8,
                  fontweight='bold')
        title_leg = self.get_leg_title(self.grup)
        lgd = plt.legend(title=title_leg, loc='center left', bbox_to_anchor=(1, 0.5))
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

        self.pdf_pages.savefig(fig, bbox_extra_artists=(lgd, sub), bbox_inches='tight')
        del fig

    def plot_combin_inhibition_recursive(self, combination, grup, data_grup, std_data_grup, poistioni, positionk,
                                         info1=0, info2=0):

        self.data_plot = combination.data_plot
        self.comb_name_per_group = combination.comb_name_per_group[poistioni][positionk]
        self.effect_per_group = combination.effect_per_group[poistioni][positionk]
        self.synergy_score_per_group = combination.synergy_score_per_group[poistioni][positionk]

        self.control_info = combination.control_info

        group = grup.pop(0)
        data_sel = data_grup.pop(0)
        std_sel = std_data_grup.pop(0)

        y_axis_name = self.set_y_axis_name(info1)
        y_axis_min, y_axis_max = self.set_y_axis_limt(info2, np.max(data_sel.data_grup), np.min(data_sel.data_grup))
        counter = None

        for i in range(len(group.concentration)):
            grup_new = grup.copy()
            fig = plt.figure(figsize=(8, 8), dpi=80)
            label_titel = "compound: " + group.name
            titel_page = str(label_titel)
            ax = plt.axes()

            label_d = str(self.control_info[3]) + " " + self.control_info[-2]  # control
            post_control = self.control_info[-1]
            plt.plot(self.elapse, self.data_plot[post_control], label=label_d, color="black")
            label_a = str(group.name) + " " + ' '.join(map(str, group.concentration[i][:-1])) + " " + str(group.unidade)
            label_b = str(self.comb_name_per_group[i])
            plt.plot(self.elapse, data_sel.data_grup[:, i], label=label_a, color="darkblue")
            plt.fill_between(self.elapse, data_sel.data_grup[:, i] - std_sel.data_grup[:, i],
                             data_sel.data_grup[:, i] + std_sel.data_grup[:, i],
                             alpha=0.5, edgecolor="lightgray", facecolor='lightgray')

            fr = []

            fr, counter, counter_old = self.plot_info(fr, grup_new, counter)

            col = ColorsHtsplots(len(fr))
            self.colormap = col.map
            for j in range(len(fr)):
                label_tem = str(grup[j].name) + " " + \
                            str(grup[j].concentration[counter_old[j]][:-1][0]) + " " + str(grup[j].unidade)

                plt.plot(self.elapse, data_grup[j].data_grup[:, counter_old[j]],
                         label=label_tem, color=self.colormap[j])
                plt.fill_between(self.elapse, data_grup[j].data_grup[:, counter_old[j]] -
                                 std_data_grup[j].data_grup[:, counter_old[j]],
                                 data_grup[j].data_grup[:, counter_old[j]] +
                                 std_data_grup[j].data_grup[:, counter_old[j]],
                                 alpha=0.5, edgecolor="lightgray", facecolor='lightgray')

            if len(self.effect_per_group[:, i]) == 1 and info1 == 2 and self.synergy_method == 0:
                # this is only for Bliss method
                # a = [element for tupl in self.predicted_per_group[k][matrix_count] for element in tupl]
                plt.plot(self.elapse, np.reshape(self.effect_per_group,
                                                 (self.effect_per_group.size, 1)),
                         label=label_b, color="black",
                         linestyle="--")
            elif self.synergy_method == 0:
                if info1 == 2:
                    plt.plot(self.elapse, self.effect_per_group[:, i], label=label_b, color="black",
                             linestyle="--")
                elif info1 == 1:
                    plt.plot(self.elapse, 0-self.effect_per_group[:, i], label=label_b, color="black",
                             linestyle="--")
            if info1 != 0:

                text_titel_bi = [self.synergy_text + ' \n']
                te = self.get_addition_text(text_titel_bi,  self.synergy_score_per_group,
                                            self.time_selected, self.time_position, i)

                tr = plt.text(self.elapse[-1], y_axis_min + 0.5, te, bbox=dict(facecolor='none', edgecolor='silver'))
            # text_titel_hsa = ["HSA_score :\n"]

            # hsa_text = plt.text(80, -0.4, text_titel_hsa, bbox=dict(facecolor='none', edgecolor='silver'))
            plt.xlabel("time (h)")
            plt.ylabel(y_axis_name)
            plt.ylim(y_axis_min, y_axis_max)
            sub = plt.suptitle(self.plot_titel, fontsize=12)
            plt.title(titel_page, fontsize=8, fontweight='bold')
            title_leg = self.get_leg_title(grup[0])
            # title_leg = grup[0].cell_name + "  " + self.std_type
            lgd = plt.legend(title=title_leg, loc='center left', bbox_to_anchor=(1, 0.5))  # , borderaxespad=0)
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            if info1 != 0:
                self.pdf_pages.savefig(fig, bbox_extra_artists=(lgd, sub, tr), bbox_inches='tight')
            else:
                self.pdf_pages.savefig(fig, bbox_extra_artists=(lgd, sub), bbox_inches='tight')

            del fig

    def plot_perturbagem(self, grup, control, data, std, info1, info2, info3, timeselected, timeposition, readout,
                         readout_units):

        self.controldata = control
        self.grup = grup
        self.data_plot = data
        self.std_plot = std
        self.readout = readout
        self.readout_units = readout_units

        y_axis_name = self.set_y_axis_name(info1)
        y_axis_min, y_axis_max = self.set_y_axis_limt(info2, np.max(self.data_plot), np.min(self.data_plot))
        fig = plt.figure(figsize=(8, 8), dpi=80)
        ax = plt.axes()
        # label_titel = "Cell lines: " + grup.cell_name

        label_b = str(self.controldata[3])  # + "_" + self.control_info[-2]
        # titel_page = str(label_titel)

        plt.plot(self.elapse, self.data_plot[self.controldata[-1]], label=label_b, color="black")
        if len(self.std_plot[0].shape) == 2:
            leg_titel = grup.cell_name + " " + grup.seeding.replace("per", "/") + "\n" + \
                        self.std_info_plot
            plt.fill_between(self.elapse, self.controlstd[0],
                             self.controlstd[1], alpha=0.5, edgecolor='black', facecolor='silver')

        else:
            # plt.fill_between(self.elapse, self.controldata - self.controlstd,
            #                  self.controldata + self.controlstd, alpha=0.5, edgecolor='black', facecolor='silver')
            leg_titel = grup.cell_name + " " + grup.seeding.replace("per", "/") + "\n" + \
                        self.std_type
            plt.fill_between(self.elapse, self.data_plot[self.controldata[-1]] - self.std_plot[self.controldata[-1]],
                             self.data_plot[self.controldata[-1]] + self.std_plot[self.controldata[-1]],
                             alpha=0.5, edgecolor='black', facecolor='silver')
        col = ColorsHtsplots(len(grup.concentration)*5)
        self.colormap = col.map

        for i in range(len(grup.concentration)):
            post1 = int(grup.concentration[i][-1])
            label_a = str(grup.name) + " " + str(grup.concentration[i][0]) + " " + grup.unidade.replace("per", "/")

            plt.plot(self.elapse, self.data_plot[post1], label=label_a, color=self.colormap[i])

            if len(self.std_plot[0].shape) == 2:
                plt.fill_between(self.elapse, self.std_plot[post1][0],
                                 self.std_plot[post1][1],
                                 alpha=0.5, edgecolor=self.colormap[i], facecolor='silver')
                # if info3 == 0:
                #     text_titel_bi = ["perturbagem effect (confluency):\n"]
                #     te = self.get_perturbagem_text(text_titel_bi, self.data_plot[post1], timeselected, timeposition)
                #     tr = plt.text(self.elapse[-1], 5, te, bbox=dict(facecolor='none', edgecolor='silver'))
            else:
                plt.fill_between(self.elapse, self.data_plot[post1] - self.std_plot[post1],
                                 self.data_plot[post1] + self.std_plot[post1],
                                 alpha=0.5, edgecolor=self.colormap[i], facecolor='silver')
                # if info3 == 0:
                #     # text_titel_bi = ["condition effect " + self.readout + self.readout_units + "\n"]
                #     # te = self.get_perturbagem_text(text_titel_bi, self.data_plot[post1], timeselected, timeposition)
                #     # tr = plt.text(self.elapse[-1], y_axis_min + 5, te, bbox=dict(facecolor='none', edgecolor='silver'))
                #     leg = ()
                if info3 != 0:
                    text_titel_bi = [y_axis_name + "\n"]
                    te = self.get_perturbagem_text(text_titel_bi, self.data_plot[post1], timeselected, timeposition)
                    tr = plt.text(self.elapse[-1], y_axis_min - 0.2, te, bbox=dict(facecolor='none',
                                                                                   edgecolor='silver'))

        plt.xlabel("time (h)")
        plt.ylabel(y_axis_name)
        plt.ylim(y_axis_min, y_axis_max)
        sub = plt.suptitle(self.plot_titel, fontsize=12)
        # plt.title(titel_page, fontsize=8,
        #           fontweight='bold')

        lgd = plt.legend(title=leg_titel, loc='center left', bbox_to_anchor=(1, 0.5))
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        if info3 == 1:
            self.pdf_pages.savefig(fig, bbox_extra_artists=(lgd, sub, tr), bbox_inches='tight')
        else:
            self.pdf_pages.savefig(fig, bbox_extra_artists=(lgd, sub), bbox_inches='tight')

    # plotting grwoth rate
    def plot_grwothrate(self, forplotting, grup, i, k, info1=0, info2=0):

        self.plotdata = forplotting.grresults[i][k]
        self.header = forplotting.grheader[i][k]
        self.elapsegrwothrate = self.elapse
        self.extremepointspositions = forplotting.extremepointspositions
        self.control_info = forplotting.control_info
        # self.std_plot = forplotting.std_plot
        self.grup = grup
        # self.extremepointspositions = [0, -1]
        # self.std_plot = forplotting.std_plot
        # self.readout = 'grwoth rate'
        self.readout_units = forplotting.readout_units

        # y_axis_name = self.set_y_axis_name(info1)
        y_axis_min, y_axis_max = self.set_y_axis_limt(4, np.max(self.datanormtozero),
                                                      np.min(self.datanormtozero))
        y_axis_min2, y_axis_max2 = self.set_y_axis_limt(3, -1, 2)

        # label_titel = "GR compound: " + self.header[0][3]
        xs = np.linspace(self.elapsegrwothrate[1], self.elapsegrwothrate[-1], int(1e5))
        linetype = [0 for i in range(len(xs[self.extremepointspositions[0]:self.extremepointspositions[-1]]))]
        self.colormap = list(iter(cm.tab20b(np.linspace(0, 1, len(self.plotdata)))))
        cout = 0
        if len(self.colormap) < len(self.plotdata):
            col = ColorsHtsplots(len(self.plotdata))
            self.colormap = col.map
        for i in range(len(self.plotdata)-1):
            fig, axs = plt.subplots(1, 2, figsize=(15, 5))
            for plo in range(2):
                if plo == 0:
                    post_control = self.control_info[-1]
                    label_b = str(self.control_info[3]) + " " + self.control_info[-2]
                    axs[plo].plot(self.elapse, self.datanormtozero[post_control], label=label_b, color="black")
                    post1 = int(grup.concentration[i][-1])
                    label_a = str(grup.name) + " " + str(grup.concentration[i][0]) + " " + grup.unidade
                    axs[plo].plot(self.elapse, self.datanormtozero[post1], label=label_a, color=self.colormap[i])
                    # axs[plo].set_title('growth normal ' + label_a, fontsize=8)
                    axs[plo].set_ylim(y_axis_min, y_axis_max)
                    lgd2 = axs[plo].legend(title='growth normal',
                                           loc='center left', bbox_to_anchor=(1, 0.5))
                else:

                    ind = np.where(xs == 1.0)

                    label_b = '_'.join(self.header[0][-3:-1])
                    axs[plo].plot(xs[self.extremepointspositions[0]:self.extremepointspositions[-1]],
                                    self.plotdata[0][0][self.extremepointspositions[0]:self.extremepointspositions[-1]],
                                    label=label_b, color='black')

                    label_a = '_'.join(self.header[i+1][-3:-1])
                    axs[plo].plot(xs[self.extremepointspositions[0]:self.extremepointspositions[-1]],
                                    self.plotdata[i+1][0][self.extremepointspositions[0]:self.extremepointspositions[-1]],
                                    label=label_a, color=self.colormap[i])
                    axs[plo].plot(xs[self.extremepointspositions[0]:self.extremepointspositions[-1]],
                                     linetype, '--', color="silver")
                    # axs[plo].set_title('growth rate ' + label_a, fontsize=8)
                    axs[plo].set_ylim(y_axis_min2, y_axis_max2)
                    # axs[plo].yticks(np.arange(y_axis_min, y_axis_max))

                    lgd = axs[plo].legend(title='growth rate ',
                                          loc='center left', bbox_to_anchor=(1, 0.5))
            cout += 1
            sub = plt.suptitle(self.plot_titel, fontsize=12)
            fig.tight_layout()
            fig.subplots_adjust(top=0.88)
            # plt.show()
            self.pdf_pages.savefig(fig, bbox_extra_artists=(lgd2, lgd, sub), bbox_inches='tight')
        del fig

    def plot_grwothratecombination(self, forplotting, grup, i, k, info1=0, info2=0):
        self.control_info = forplotting.control_info
        self.plotdata = forplotting.grresultscombination[i][k]
        self.header = forplotting.grheadercombination[i][k]
        self.elapsegrwothrate = self.elapse
        self.synergy_score_per_group = forplotting.synergy_score_per_group[i][k]
        xs = np.linspace(self.elapsegrwothrate[1], self.elapsegrwothrate[-1], int(1e5))
        self.extremepointspositions = forplotting.extremepointspositions

        self.grup = grup
        self.get_gr()
        counter = None

        self.readout_units = forplotting.readout_units

        y_axis_min, y_axis_max = self.set_y_axis_limt(4, np.max(self.datanormtozero),
                                                      np.min(self.datanormtozero))
        y_axis_min2, y_axis_max2 = self.set_y_axis_limt(3, -1, 2)

        label_titel = self.grup[0].cell_name + " combination " # + self.header[0][0][3]
        linetype = [0 for i in range(len(xs[self.extremepointspositions[0]:self.extremepointspositions[-1]]))]
        self.colormap = list(iter(cm.tab20b(np.linspace(0, 1, len(self.plotdata[0])))))

        if len(self.colormap) < len(self.plotdata):
            col = ColorsHtsplots(len(self.plotdata))
            self.colormap = col.map
        grup_new = self.grup.copy()
        group = grup.pop(0)

        count = 0
        for i in range(len(self.plotdata[0])-1):
            fr = []
            fr, counter, counter_old = self.plot_infogr(fr, grup, counter)
            fig, axs = plt.subplots(2, 2, figsize=(15, 5), gridspec_kw={'height_ratios': [2, 1]}) #, tight_layout=True
            axs[1, 0].axis('off')  # , 20,5
            axs[1, 1].axis('off')  # , 20,5
            countplotleft = 0
            countplotright = 0
            for plo in range(2):
                if plo == 0:
                    post_control = self.control_info[-1]
                    label_b = str(self.control_info[3]) + " " + self.control_info[-2]
                    axs[0, plo].plot(self.elapse, self.datanormtozero[post_control], label=label_b, color="black")
                    post1 = int(grup_new[0].concentration[i][-1])
                    label_a = str(grup_new[0].name) + " " + ' '.join(map(str, grup_new[0].concentration[i][:-1])) +\
                              " " + grup_new[0].unidade
                    axs[0, plo].plot(self.elapse, self.datanormtozero[post1], label=label_a, color="darkblue")
                    countplotright += 1
                    col = ColorsHtsplots(len(fr))
                    self.colormap = col.map
                    for each in range(len(fr)):
                        print(counter, counter_old)
                        post2 = int(grup[each].concentration[counter_old[each]][-1])
                        label_c = str(grup[each].name) + " " + \
                                  str(grup[each].concentration[counter_old[each]][0]) + " " + \
                                  grup[each].unidade
                        axs[0, plo].plot(self.elapse, self.datanormtozero[post2], label=label_c,
                                         color=self.colormap[each])
                        countplotright +=1
                    text_titel_bi = [self.synergy_text + ' \n']
                    te1 = self.get_addition_text(text_titel_bi, self.synergy_score_per_group,
                                                 self.time_selected, self.time_position, i)
                    axs[0, plo].set_title(self.readout + 'normalized to 0 time point', fontsize=8)  # + label_a
                    lgd2 = axs[0, plo].legend(title='Norm to 0 time point', loc='center left',
                                              bbox_to_anchor=(1, 0.5), fontsize=8)
                    axs[0, plo].set_ylim(y_axis_min, y_axis_max)

                    axs[1, plo].text(0, 0.5, te1, bbox=dict(facecolor='none', edgecolor='silver'), fontsize=8)

                else:
                    textnamegr = []
                    position = []
                    label_b = '_'.join(self.header[0][0][-3:-1])
                    axs[0, plo].plot(xs[self.extremepointspositions[0]:self.extremepointspositions[-1]],
                                     self.plotdata[0][0][0][self.extremepointspositions[0]:
                                                            self.extremepointspositions[-1]],
                                     label=label_b, color='black')
                    label_a = '_'.join(self.header[0][i + 1][-3:-1])
                    textnamegr.append(label_a)
                    position.append([0, i + 1])
                    axs[0, plo].plot(xs[self.extremepointspositions[0]:self.extremepointspositions[-1]],
                                     self.plotdata[0][i + 1][0][
                                     self.extremepointspositions[0]:self.extremepointspositions[-1]],
                                     label=label_a, color="darkblue")
                    axs[0, plo].plot(xs[self.extremepointspositions[0]:self.extremepointspositions[-1]],
                                     linetype, '--', color="silver")
                    countplotleft += 1
                    if count == len(grup[-1].concentration):
                        count = 0

                    col = ColorsHtsplots(len(fr))
                    self.colormap = col.map
                    for each in range(len(fr)):
                        label_c = '_'.join(self.header[each + 1][counter_old[each] + 1][-3:-1])
                        textnamegr.append(label_c)
                        position.append([each + 1, counter_old[each]+1])
                        axs[0, plo].plot(xs[self.extremepointspositions[0]:self.extremepointspositions[-1]],
                                         self.plotdata[each + 1][counter_old[each] + 1][0][
                                         self.extremepointspositions[0]:self.extremepointspositions[-1]],
                                         label=label_c, color=self.colormap[each])

                    te = self.get_addition_grtext('grwoth rate', textnamegr, position)

                    axs[1, plo].text(0, 0, te, bbox=dict(facecolor='none', edgecolor='silver'), fontsize=8)

                    axs[0, plo].set_title('growth rate ', fontsize=8)
                    lgd = axs[0, plo].legend(title='growth rate', loc='center left', bbox_to_anchor=(1, 0.5),
                                             fontsize=8)
                    axs[0, plo].set_ylim(y_axis_min2, y_axis_max2)
                    axs[0, plo].set_yticks(np.arange(y_axis_min2, y_axis_max2, 0.5))

                    sub = plt.suptitle(label_titel, fontsize=12)
                    fig.tight_layout()
                    fig.subplots_adjust(top=0.88)
                    self.pdf_pages.savefig(fig, bbox_extra_artists=(lgd2, lgd, sub), bbox_inches='tight')
                    # del fig

    # ####

    # Bar plot methods

    def barplot_comb_inhibition_recursive(self, combination, grup, data_grup, std_data_grup, info1, info2, i, k):

        self.comb_name_per_group = combination.comb_name_per_group[i][k]
        self.effect_per_group = combination.effect_per_group[i][k]
        self.synergy_score_per_group = combination.synergy_score_per_group[i][k]

        self.data_plot = data_grup
        self.std_plot = std_data_grup

        group = grup.pop(0)
        data_sel = data_grup.pop(0)
        std_sel = std_data_grup.pop(0)

        y_axis_name = self.set_y_axis_name(info1)
        y_axis_min, y_axis_max = self.set_y_axis_limt(info2, np.max(data_sel.data_grup), np.min(data_sel.data_grup))

        col = ColorsHtsplots(len(grup[0].concentration))
        self.colormap = col.map

        for tim in range(len(self.time_selected)):
            combination_counter = 0
            counter = None

            if len(grup[0].concentration) <= 4:
                rs, cols = (1, 4)
                compound_info = "compound: " + grup[0].name + "\n" + " cell line: " + grup[0].cell_name + "\n" + \
                                "Time point" + " " + str(self.time_selected[tim]) + " h"

                fig, axs = plt.subplots(int(rs), int(cols), figsize=(30, 20))  # sharey='row',
                fig.suptitle(compound_info, fontsize=26, fontweight='bold', ha="center")

                cont = 0

                self.barplot_1row(grup, cols, tim, cont, group, combination_counter, axs, data_sel,
                                  std_sel, data_grup, std_data_grup, y_axis_name, y_axis_max, y_axis_min, counter,
                                  info1)
                self.pdf_pages.savefig(fig, bbox_inches='tight')
                del fig

            else:
                for i in range(len(grup[-1].concentration)):
                    grup = grup.copy()
                    data_grup = data_grup.copy()

                    rows = math.ceil(len(grup[0].concentration) / 4)
                    rs, cols = (rows, 4)
                    compound_info = "compound: " + group.name + "\n" + " cell line: " + group.cell_name + "\n" + \
                                    "Time point" + " " + str(int(self.time_selected[tim])) + " h"
                    if rs >= 4:
                        fig, axs = plt.subplots(int(rs), int(cols), figsize=(90, 50))  # sharey='row',
                        fig.suptitle(compound_info, fontsize=42, fontweight='bold', ha="center")
                    else:
                        fig, axs = plt.subplots(int(rs), int(cols), figsize=(45, 20))  # sharey='row',
                        fig.suptitle(compound_info, fontsize=26, fontweight='bold', ha="center")

                    cont = 0
                    combination_counter, counter, cont = self.barplot_rows(grup, rs, cols, tim, cont,
                                                                           group, combination_counter, axs, data_sel,
                                                                           std_sel, data_grup, std_data_grup,
                                                                           y_axis_name, y_axis_max, y_axis_min, counter,
                                                                           info1)
                    # plt.savefig('C:/Users/cdcarval/Dropbox (speleman lab)/Personal Lab/'
                    #             'Carol personal Lab/In silico lab/Drugging_script/paperfigures/HTSplotter_2/' +
                    #             compound_info + '.pdf')
                    self.pdf_pages.savefig(fig, bbox_inches='tight')
                    del fig

    def barplot_1row(self, grup, cols, tim, cont, group, combination_counter, axs, data_sel, std_sel, data_grup,
                     std_data_grup, y_axis_name, y_axis_max, y_axis_min, counter, info1):
        # self.adjust_combinationvalues(data_grup)
        for col in range(cols):
            if cont >= len(grup[0].concentration):
                axs[col].axis('off')
                pass
            else:
                grup_new = grup.copy()
                label_a = str(group.name) + " " + ' '.join(map(str, group.concentration[combination_counter][:-1])) \
                          + " " + str(group.unidade)
                label_b = str(self.comb_name_per_group[combination_counter])
                if info1 != 0:
                    if self.synergy_method == 0:
                        te = "Predicted value = " + \
                             str("{:.2f}".format(self.effect_per_group
                                                 [self.time_position[tim], combination_counter])) + \
                             self.synergy_text + ' = ' + str("{:.2f}".format(self.synergy_score_per_group[
                                                                                 self.time_position[tim],
                                                                                 combination_counter]))
                    else:
                        te = self.synergy_text + ' = ' + str("{:.2f}".format(self.synergy_score_per_group[
                                                                                 self.time_position[tim],
                                                                                 combination_counter]))
                fr = []

                fr, counter, counter_old = self.plot_info_barplot(fr, grup_new, counter)

                x_barplot = np.arange(len(fr) + 1)

                axs[col].bar(x_barplot[0],
                             data_sel.data_grup[self.time_position[tim], combination_counter],
                             yerr=std_sel.data_grup[self.time_position[tim], combination_counter],
                             color="darkblue")
                # horizontal line indicating the threshold
                if info1 == 2:
                    if self.synergy_method == 0:
                        #  horizontal line indicating the threshold
                        axs[col].axhline(y=self.effect_per_group[self.time_position[tim], combination_counter],
                                         linewidth=2, ls="--", color="black")
                        leg = [te, label_a]
                    else:
                        axs[col].axhline(y=self.effect_per_group[self.time_position[tim], combination_counter],
                                         linewidth=2, alpha=0.1, color='white')
                        leg = [te, label_a]

                elif info1 == 1:
                    if self.synergy_method == 0:
                        axs[col].axhline(y=0-self.effect_per_group[self.time_position[tim], combination_counter],
                                         linewidth=1.5, ls="--", color="black")
                        leg = [te, label_a]
                    else:
                        axs[col].axhline(y=0-self.effect_per_group[self.time_position[tim], combination_counter],
                                         linewidth=2, alpha=0.1, color='white')
                        leg = [te, label_a]

                for j in range(len(fr)):

                    label_tem = str(grup[j].name) + " " + \
                                str(grup[j].concentration[counter_old[j]][:-1][0]) + " " + str(grup[j].unidade)
                    leg.append(label_tem)

                    if data_grup[j].data_grup[self.time_position[tim], counter_old[j]] < 0 and info1 == 2:
                        data_grup[j].data_grup[self.time_position[tim], counter_old[j]] = 0

                    axs[col].bar(x_barplot[j + 1],
                                 data_grup[j].data_grup[self.time_position[tim], counter_old[j]],
                                 yerr=std_data_grup[j].data_grup[self.time_position[tim],
                                                                 counter_old[j]], color=self.colormap[j])

                axs[col].set_ylabel(y_axis_name, fontsize=18)
                axs[col].set_ylim(y_axis_min - 0.75, y_axis_max)
                axs[col].tick_params(axis="y", labelsize=18)
                axs[col].tick_params(axis="x", labelsize=18)

                title_leg = "  " + self.std_type

                axs[col].legend(leg, title=title_leg, fontsize=16, loc='lower center', frameon=False,
                                title_fontsize=16)

                combination_counter += 1

                cont += 1

    def barplot_rows(self, grup, rs, cols, tim, cont, group, combination_counter, axs, data_sel, std_sel, data_grup,
                     std_data_grup, y_axis_name, y_axis_max, y_axis_min, counter, info1):

        for r in range(int(rs)):
            for col in range(cols):
                if cont >= len(grup[0].concentration):
                    axs[r, col].axis('off')
                    pass
                else:
                    grup_new = grup.copy()
                    label_a = str(group.name) + " " + ' '.join(
                        map(str, group.concentration[combination_counter][:-1])) + " " + str(
                        group.unidade)
                    label_b = str(self.comb_name_per_group[combination_counter])
                    if info1 !=0:
                        if self.synergy_method == 0:
                            te = "Predicted value = " + \
                                 str("{:.2f}".format(self.effect_per_group
                                                     [self.time_position[tim], combination_counter])) + \
                                 self.synergy_text + " = " + \
                                 str("{:.2f}".format(self.synergy_score_per_group[self.time_position[tim],
                                                                                  combination_counter]))
                        else:
                            te = self.synergy_text + " = " + str("{:.2f}".format(self.synergy_score_per_group[
                                                                                     self.time_position[tim],
                                                                                     combination_counter]))
                    fr = []

                    fr, counter, counter_old = self.plot_info_barplot(fr, grup_new, counter)

                    x_barplot = np.arange(len(fr) + 1)

                    axs[r, col].bar(x_barplot[0],
                                    data_sel.data_grup[self.time_position[tim], combination_counter],
                                    yerr=std_sel.data_grup[self.time_position[tim], combination_counter],
                                    color="darkblue")

                    if info1 == 2:
                        if self.synergy_method == 0:
                        #  horizontal line indicating the threshold
                            axs[r, col].axhline(y=self.effect_per_group[self.time_position[tim], combination_counter],
                                                linewidth=2, ls="--", color="black")
                            leg = [te, label_a]
                        else:
                            axs[r, col].axhline(y=self.effect_per_group[self.time_position[tim], combination_counter],
                                                linewidth=2, alpha=0.1, color='white')
                            leg = [te, label_a]

                    elif info1 == 1:
                        if self.synergy_method == 0 or self.synergy_method ==2:
                            axs[r, col].axhline(y=0-self.effect_per_group[self.time_position[tim], combination_counter],
                                                linewidth=1.5, ls="--", color="black")
                            leg = [te, label_a]
                        else:
                            axs[r, col].axhline(y=0 - self.effect_per_group[self.time_position[tim],
                                                                            combination_counter],
                                                linewidth=1.5, alpha=0.1, color="white")
                            leg = [te, label_a]

                    else:
                        leg = [label_a]
                    # color = iter(cm.tab20b(np.linspace(0, 1, len(fr))))
                    for j in range(len(fr)):
                        # c = next(color)
                        label_tem = str(grup[j].name) + " " + \
                                    str(grup[j].concentration[counter_old[j]][:-1][0]) + " " + str(grup[j].unidade)
                        leg.append(label_tem)

                        if data_grup[j].data_grup[self.time_position[tim], counter_old[j]] < 0 and info1 == 2:
                            data_grup[j].data_grup[self.time_position[tim], counter_old[j]] = 0

                        axs[r, col].bar(x_barplot[j + 1],
                                        data_grup[j].data_grup[self.time_position[tim], counter_old[j]],
                                        yerr=std_data_grup[j].data_grup[self.time_position[tim],
                                                                        counter_old[j]], color=self.colormap[j])

                    axs[r, col].set_ylabel(y_axis_name, fontsize=18)
                    axs[r, col].set_ylim(y_axis_min - 0.75, y_axis_max)
                    axs[r, col].tick_params(axis="y", labelsize=18)
                    axs[r, col].tick_params(axis="x", labelsize=18)

                    title_leg = "  " + self.std_type

                    axs[r, col].legend(leg, title=title_leg, fontsize=16, loc='lower left', frameon=False,
                                       title_fontsize=16)

                    combination_counter += 1

                    cont += 1

        return combination_counter, counter, cont

    def barplot_geneticperturbgen_recursive(self, grup, data_grup, std_data_grup, info1, info2, timeselected,
                                            timeposition, readout, readout_units):
        self.readout = readout
        self.readout_units = readout_units
        y_axis_name = self.set_y_axis_name(info1)
        y_axis_min, y_axis_max = self.set_y_axis_limt(info2, np.max(data_grup), np.min(data_grup))
        # self.adjust_combinationvalues(data_grup)

        counter = None

        compound_info = " cell line: " + grup[0].cell_name + "\n" + \
                        "Time point" + " " + str(timeselected[0]) + " h"

        fig, axs = plt.subplots(1, 1, figsize=(30, 20))  # sharey='row',
        fig.suptitle(compound_info, fontsize=26, fontweight='bold', ha="center")

        cont = 0

        fr = []

        fr, counter, counter_old = self.plot_info_barplot(fr, grup, counter)
        data = []
        std = []
        leg = []
        x_barplot = np.arange(len(grup))
        grup_valore = []
        col = ColorsHtsplots(len(grup))
        self.colormap = col.map
        for j in range(len(grup)):
            te = 'effect = ' + (str("{:.2f}".format(data_grup[int(grup[j].concentration[0][-1])][0]))) + \
                 "_" + grup[j].name

            axs.bar(x_barplot[j], data_grup[int(grup[j].concentration[0][-1])][0],
                    yerr=std_data_grup[int(grup[j].concentration[0][-1])][0], color=self.colormap[j])
            leg.append(te)

        axs.set_ylabel(y_axis_name, fontsize=18)
        axs.set_ylim(-0.75, y_axis_max)
        axs.tick_params(axis="y", labelsize=18)
        axs.tick_params(axis="x", labelsize=18)

        title_leg = "  " + self.std_type

        axs.legend(leg, ncol=6, title=title_leg, fontsize=16, loc=2, frameon=False,
                   title_fontsize=16)

        self.pdf_pages.savefig(fig, bbox_inches='tight')
        del fig

    def barplot_geneticchemicalperturbgen_recursive(self, grup, data_grup, std_data_grup, info1, info2, timeselected,
                                                    timeposition, readout, readout_units):
        self.readout = readout
        self.readout_units = readout_units
        y_axis_name = self.set_y_axis_name(info1)
        y_axis_min, y_axis_max = self.set_y_axis_limt(info2, np.max(data_grup), np.min(data_grup))
        # self.adjust_combinationvalues(data_grup)

        # for tim in range(len(timeselected)):
        combination_counter = 0
        counter = None

        compound_info = " cell line: " + grup[0].cell_name + "\n" + \
                        "Time point" + " " + str(timeselected[0]) + " h"

        fig, axs = plt.subplots(1, 1, figsize=(30, 20))  # sharey='row',
        fig.suptitle(compound_info, fontsize=26, fontweight='bold', ha="center")

        cont = 0
        col = ColorsHtsplots(len(grup[0].concentration))
        self.colormap = col.map
        # color = iter(cm.Set2(np.linspace(0, 2, len(grup) * 5)))

        fr = []

        fr, counter, counter_old = self.plot_info_barplot(fr, grup, counter)
        data = []
        std = []
        leg = []
        x_barplot = np.arange(len(grup))
        for j in range(len(grup)):
            # c = next(color)

            te = 'effect = ' + (str("{:.2f}".format(data_grup[int(grup[j].concentration[0][-1])][0]))) + \
                 "_" + grup[j].name

            # data.append(data_grup[grup[i].concentration[0][-1][0]])
            # std.append(data_grup[grup[i].concentration[0][-1][0]])
            axs.bar(x_barplot[j], data_grup[int(grup[j].concentration[0][-1])][0],
                    yerr=std_data_grup[int(grup[j].concentration[0][-1])][0], color=self.colormap[j])

            leg.append(te)
        axs.set_ylabel(y_axis_name, fontsize=18)
        axs.set_ylim(-0.75, y_axis_max)
        axs.tick_params(axis="y", labelsize=18)
        axs.tick_params(axis="x", labelsize=18)

        title_leg = "  " + self.std_type

        axs.legend(leg, ncol=4, title=title_leg, fontsize=16, loc=2, frameon=False,
                   title_fontsize=16)

        # self.barplot_geneticperturbagen(grup, branch, tim, cont, data_grup,
        #                                 std_data_grup, y_axis_name, y_axis_max, counter, timeposition)
        self.pdf_pages.savefig(fig, bbox_inches='tight')
        del fig

    def plot_info_barplot(self, fr, grup, counter, ind=0):

        if len(grup) > ind:

            if counter is None:
                counter = np.zeros(len(grup), dtype=np.int16)

            group = grup[ind]

            for i in range(ind, -1, -1):
                if counter[i] == grup[i].concentration.shape[0]:
                    counter[i] = 0
                    if i > 0:
                        counter[i - 1] += 1

            tem = str(group.name) + " " + str(group.concentration[counter[ind], 0]) + " " + str(group.unidade)
            fr.append(tem)
            fr, counter, counter_old = self.plot_info(fr, grup, counter, ind + 1)

            if ind == len(counter) - 1:
                counter_old = counter.copy()
                counter[ind] += 1

            return fr, counter, counter_old

        else:
            counter_old = counter

            return fr, counter, counter_old

    @staticmethod
    def adjust_combinationvalues(data_grup):
        for i in range(len(data_grup)):
            if data_grup[i] < 0:
                data_grup[i] = 0

    # ####

    # # HEATMAP
    def heat_map_overtime_bi(self, combination, positinioni, positionk, grup):

        self.synergy_score_per_group = combination.synergy_score_per_group[positinioni][positionk]
        self.comb_name_per_group = combination.comb_name_per_group[positinioni][positionk]

        title = grup[0].cell_name + self.synergy_text + " " + grup[0].name

        a_ff = "{:.2f}".format(np.min(self.synergy_score_per_group))
        min_string_a = str(a_ff)
        b_ff = "{:.2f}".format(np.max(self.synergy_score_per_group))
        max_string_a = str(b_ff)
        string_completa_a = '>0 syne; <0- anta' + " " + "min = " + min_string_a + " & " + "max = " + max_string_a

        fig, ax = plt.subplots(figsize=(40, 20))

        ax.set_xticks(np.arange(self.elapse.shape[0]))
        ax.set_yticks(np.arange(len(self.comb_name_per_group)))
        ax.set_xticklabels(self.elapse, fontsize=15)
        ax.set_yticklabels(self.comb_name_per_group, fontsize=15)

        # if self.synergy_method == 2 or self.synergy_method == 3:
        #     heatplot = ax.imshow(self.synergy_score_per_group.T, vmin=a_ff, vmax=b_ff, cmap='RdBu',
        #                          label="teste")
        # else:
        heatplot = ax.imshow(self.synergy_score_per_group.T, vmin=-1, vmax=1, cmap='RdBu',
                             label="teste")

        plt.setp(ax.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor")
        cbar = fig.colorbar(heatplot, ticks=[-1, 0.5, 0, -0.5, 1],
                            fraction=0.15, pad=0.05, shrink=0.5, aspect=20)
        # cbar = fig.colorbar(heatplot, fraction=0.15, pad=0.05, shrink=0.5, aspect=20)

        cbar.ax.tick_params(labelsize=20)
        cbar.set_label(string_completa_a, labelpad=100, rotation=270,
                       fontsize=20)

        ax.set_title(title, fontsize=20)
        fig.tight_layout()

        self.pdf_pages.savefig(fig)

        del fig

    def heat_map_selec_time_bi_biDim(self, combination, grup, i, k):
        combname = np.round(combination.synergy_score_per_group[i][k], decimals=3)
        self.synergy_score_per_group = combination.synergy_score_per_group[i][k]
        self.comb_name_per_group = combination.comb_name_per_group[i][k]

        for j in range(len(self.time_selected)):
            fig = plt.figure(figsize=(25, 15))
            ti_posit = self.time_position[j]

            title = grup[0].cell_name + self.synergy_text + " " \
                    + grup[0].name + " " + str(self.time_selected[j])

            tem3d, tem, name_anotation = self.get_temp_bi_score(grup, self.synergy_score_per_group[ti_posit, :],
                                                                combname[ti_posit, :])
            a_ff = "{:.2f}".format(np.min(tem))
            min_string_a = str(a_ff)
            b_ff = "{:.2f}".format(np.max(tem))
            max_string_a = str(b_ff)
            string_completa_a = '>0 syne; <0- anta' + " " + "min = " + min_string_a + " & " + "max = " + max_string_a

            name_x, name_y = self.get_name(self.comb_name_per_group, 1)

            ax = fig.add_subplot(1, 2, 1)

            ax.set_xticks(np.arange(name_x.shape[0]))
            ax.set_yticks(np.arange(name_y.shape[0]))
            ax.set_xticklabels(name_x, fontsize=15)
            ax.set_yticklabels(name_y, fontsize=15)

            heatplot = ax.imshow(tem, vmin=-1, vmax=1, cmap='RdBu')

            plt.setp(ax.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor")

            cbar = fig.colorbar(heatplot, ticks=[-1, 0.5, 0, -0.5, 1],
                                fraction=0.15, pad=0.05, shrink=0.5, aspect=20)

            for i in range(tem.shape[0]):
                for j in range(tem.shape[1]):
                    text = ax.text(j, i, "{:.2f}".format(tem[i, j]),
                                   ha="center", va="center", color="black")
            cbar.ax.tick_params(labelsize=18)
            cbar.set_label(string_completa_a, labelpad=100, rotation=-90, fontsize=20)
            fig.tight_layout()
            # plt.show()
            ax.set_title(title, fontsize=20)
            # if self.synergy_method != 2 and self.synergy_method != 0:
            if len(grup[1].concentration[:, 0]) != 1 and len(grup[2].concentration[:, 0]) != 1:

                ax = fig.add_subplot(1, 2, 2, projection="3d")

                xmax = max(grup[1].concentration[:, 0])
                xmin = min(grup[1].concentration[:, 0])
                ymax = max(grup[0].concentration[:, 0])
                ymin = min(grup[0].concentration[:, 0])

                conc1 = np.zeros(len(grup[1].concentration[:, 0]) * len(grup[2].concentration[:, 0]))
                conc2 = np.zeros(len(grup[1].concentration[:, 0]) * len(grup[2].concentration[:, 0]))
                if len(grup[1].concentration) > len(grup[2].concentration):
                    for p in range(len(grup[2].concentration)):
                        conc1[p * len(grup[1].concentration[:, 0]):(p + 1) * len(grup[1].concentration[:, 0])] = \
                            grup[1].concentration[:, 0]
                        conc2[p * len(grup[1].concentration[:, 0]):(p + 1) * len(grup[1].concentration[:, 0])] = \
                            grup[1].concentration[p, 0]
                elif len(grup[1].concentration) < len(grup[2].concentration):
                    for p in range(len(grup[1].concentration)):
                        conc1[p * len(grup[2].concentration[:, 0]):(p + 1) * len(grup[2].concentration[:, 0])] = \
                            grup[2].concentration[:, 0]
                        conc2[p * len(grup[2].concentration[:, 0]):(p + 1) * len(grup[2].concentration[:, 0])] = \
                            grup[2].concentration[p, 0]
                else:
                    for p in range(len(grup[1].concentration)):
                        conc1[p * len(grup[1].concentration[:, 0]):(p + 1) * len(grup[1].concentration[:, 0])] = \
                            grup[1].concentration[:, 0]
                        conc2[p * len(grup[1].concentration[:, 0]):(p + 1) * len(grup[1].concentration[:, 0])] = \
                            grup[1].concentration[p, 0]

                rbf = inter.Rbf(np.log10(conc1), np.log10(conc2), tem3d)

                x = np.linspace(np.log10(xmin), np.log10(xmax), 100)
                y = np.linspace(np.log10(ymin), np.log10(ymax), 100)

                arr2 = np.zeros((100 * 100, 3))
                # make interpolation
                for u in range(len(x)):
                    for lin in range(len(y)):
                        arr2[u * len(x) + lin, :2] = [x[u], y[lin]]

                arr2[:, 2] = rbf(arr2[:, 0], arr2[:, 1])

                xx, yy = np.meshgrid(x, y)
                surf = ax.plot_surface(xx, yy, np.reshape(arr2[:, 2], (100, 100)), vmin=-1, vmax=1,
                                       cmap="RdBu", antialiased=True)

                cbar = fig.colorbar(heatplot, ticks=[-1, 0.5, 0, -0.5, 1], fraction=0.15,
                                    pad=0.08, shrink=0.5, aspect=20)

                cbar.ax.tick_params(labelsize=18)
                cbar.set_label(string_completa_a, labelpad=100, rotation=-90,
                               fontsize=20)

                ax.grid(b=None)
                ax.xaxis.pane.fill = False
                ax.yaxis.pane.fill = False
                ax.zaxis.pane.fill = False

                name_xlabel = str(grup[1].name) + "\n" + " Concentration range log$_{10}$ from: " + "\n" + "[" \
                              + str(min(grup[1].concentration[:, 0])) \
                              + " " + str(max(grup[1].concentration[:, 0])) + "] " + str(grup[1].unidade)

                name_ylabel = str(grup[2].name) + "\n" + " Concentration range log$_{10}$ from: " + "\n" + "[" \
                              + str(min(grup[2].concentration[:, 0])) \
                              + " " + str(max(grup[2].concentration[:, 0])) + "] " + str(grup[2].unidade)

                ax.set_xlabel(name_xlabel, fontsize=18, labelpad=30)

                ax.set_ylabel(name_ylabel, fontsize=18, labelpad=30)

                ax.set_zlabel('Synergism Score', fontsize=18, labelpad=50)
                ax.view_init(elev=15, azim=-70)

                ax.set_zlim(-1, 1)
                ax.tick_params(axis='z', which='major', pad=18)
                ax.tick_params(axis='y', which='major', pad=8)
                ax.zaxis.set_tick_params(labelsize=16)
                ax.yaxis.set_tick_params(labelsize=16)
                ax.xaxis.set_tick_params(labelsize=16)

                fig.tight_layout()
                fig.subplots_adjust(top=0.88)
            # else:
            #     ax = fig.add_subplot(1, 2, 2)
            #
            #     ax.set_xticks(np.arange(name_x.shape[0]))
            #     ax.set_yticks(np.arange(name_y.shape[0]))
            #     ax.set_xticklabels(name_x, fontsize=15)
            #     ax.set_yticklabels(name_y, fontsize=15)
            #
            #     # heatplot = ax.imshow(tem, vmin=-1, vmax=1, cmap='RdBu')
            #     heatplot = ax.imshow(tem, vmin=-1, vmax=1, cmap='RdBu')
            #     # for ynum in range(len(name_x)):
            #     #     for xnum in range(len(name_y)):
            #     #         text = ax.text(ynum, xnum, np.round(tem[xnum, ynum], 2),
            #     #                        ha="center", va="center", color="black")
            #
            #     plt.setp(ax.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor")
            #     cbar = fig.colorbar(heatplot, ticks=[-1, 0, 1],
            #                         fraction=0.15, pad=0.05, shrink=0.5, aspect=20)
            #     # cbar = fig.colorbar(heatplot, ticks=[-1, 0.5, 0, -0.5, 1], fraction=0.15, pad=0.05, shrink=0.5, aspect=20)
            #     # cbar = fig.colorbar(heatplot, ticks=[-1, 0, 1], fraction=0.15, pad=0.05, shrink=0.5, aspect=20)
            #     for i in range(name_anotation.shape[0]):
            #         for j in range(name_anotation.shape[1]):
            #             text = ax.text(j, i, name_anotation[i, j],
            #                            ha="center", va="center", color="black")
            #
            #     cbar.ax.tick_params(labelsize=18)
            #     cbar.set_label(string_completa_a, labelpad=100, rotation=-90, fontsize=20)
            #     fig.tight_layout()
            #     # plt.show()
            #     ax.set_title(title, fontsize=20)
            self.pdf_pages.savefig(fig)
            # plt.savefig('C:/Users/cdcarval/Dropbox (speleman lab)/Personal Lab/'
            #             'Carol personal Lab/In silico lab/Drugging_script/paperfigures/HTSplotter_2/' + title + '.pdf')
            del fig

    def get_temp_bi_score(self, grup, bi_score_per_group, combname):
        # i=cell line position
        # j=combination matrix position

        temp3d = bi_score_per_group

        concentrationposition = []
        concentrationposition2 = []
        # for i in range(len(grup[1].concentration)):
        # concentrationposition.append([grup[0].concentration[:, k] for k in range(grup[0].concentration.shape[1]-1)])
        # for i in range(len(concentrationposition[0][0])):
        #     string = str(concentrationposition[0][0][i]) + '_' + str(concentrationposition[0][1][i])
        #     concentrationposition2.append(string)
        # # concentrationposition = []
        # concentrationposition2 = np.asarray(concentrationposition2)

        if len(grup) >= 4:
            temp = np.reshape(temp3d, (len(grup[1].concentration), len(grup[2].concentration),
                                       len(grup[3].concentration))).T
        else:
            temp = np.reshape(temp3d, (len(grup[1].concentration), len(grup[2].concentration))).T
            # name_bi = np.reshape(concentrationposition2, (len(grup[1].concentration), len(grup[2].concentration))).T
            name_bi = np.reshape(np.asarray(combname), (len(grup[1].concentration), len(grup[2].concentration))).T
        return temp3d, temp, name_bi

    def heat_map_perturbation(self, grup, data, name, control):
        self.data_plot = data
        self.compounds_grups = name
        self.control_info = control

        title = "Genetic perturbagem" + "_relative to the control"

        min_string_a = str("{:.2f}".format(np.min(self.data_plot)))
        max_string_a = str("{:.2f}".format(np.max(self.data_plot)))
        string_completa_a = "min effect = " + min_string_a + " & " + "max effect= " + max_string_a

        fig, ax = plt.subplots(figsize=(40, 20))

        lis = np.zeros((len(self.compounds_grups), 2))
        for k in range(lis.shape[0]):
            lis[k, :] = [k, self.compounds_grups[k][0][-1]]
        lis.view("float, float").sort(axis=0, kind="mergesort", order=['f1'])
        self.compounds_grups = [self.compounds_grups[int(k)] for k in lis[:, 0]]
        com_name = []
        for u in self.compounds_grups:
            com_name.append(u[0][-3])

        ax.set_xticks(np.arange(self.elapse.shape[0]))
        ax.set_yticks(np.arange(len(self.compounds_grups)))
        ax.set_xticklabels(self.elapse, fontsize=15)
        ax.set_yticklabels(com_name, fontsize=15)
        data_plot_h = self.data_plot.copy()
        del data_plot_h[self.control_info[-1]]
        heatplot = ax.imshow(np.asarray(data_plot_h), vmin=-1, vmax=1, cmap='RdBu',
                             label="teste")

        plt.setp(ax.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor")

        cbar = fig.colorbar(heatplot, fraction=0.15, pad=0.05, shrink=0.5, aspect=20)  # fraction = 0.0046, pad = 0.04

        cbar.ax.tick_params(labelsize=20)
        cbar.set_label(string_completa_a, labelpad=100, rotation=270,
                       fontsize=20)  # labelpad >0 fica mais long, <0 fica mais perto da barra

        ax.set_title(title, fontsize=20)
        fig.tight_layout()

        self.pdf_pages.savefig(fig)
        del fig

    # on an off plots
    def plot_inhibition_on_off(self, grup, data, std, info1, info2, readout, readout_units):
        self.data_plot = data
        self.std_plot = std
        self.readout = readout
        self.readout_units = readout_units
        y_axis_name = self.set_y_axis_name(info1)
        y_axis_min, y_axis_max = self.set_y_axis_limt(info2, np.max(data), np.min(data))

        for i in range(len(grup[0].concentration)):
            color = iter(cm.tab20b(np.linspace(0, 1, len(grup[0].concentration) + 1)))
            # plot mehtod is only for on and off cases
            fig = plt.figure(figsize=(8, 8), dpi=80)
            ax = plt.axes()
            label_titel = "Cell lines: " + grup[0].cell_name + "compound: " + grup[0].name + " & " + grup[1].name
            label_a = grup[0].name + "_ " + str(grup[0].concentration[i][0]) + "_" + grup[0].unidade
            label_b = grup[1].name + "_ " + str(grup[1].concentration[i][0]) + "_" + grup[1].unidade

            titel_page = str(label_titel)
            post_a = int(grup[0].concentration[i, -1])
            post_b = int(grup[1].concentration[i, -1])

            c = next(color)
            plt.plot(self.elapse, self.data_plot[post_a], label=label_a, color=c)
            plt.fill_between(self.elapse, self.data_plot[post_a] - self.std_plot[post_a],
                             self.data_plot[post_a] + self.std_plot[post_a],
                             alpha=0.5, edgecolor="lightgray", facecolor='lightgray')
            c = next(color)
            plt.plot(self.elapse, self.data_plot[post_b], label=label_b, color=c)
            plt.fill_between(self.elapse, self.data_plot[post_b] - self.std_plot[post_b],
                             self.data_plot[post_b] + self.std_plot[post_b],
                             alpha=0.5, edgecolor="lightgray", facecolor='lightgray')
            plt.xlabel("time (h)")
            plt.ylabel(y_axis_name)
            plt.ylim(y_axis_min, y_axis_max)
            sub = plt.suptitle(self.plot_titel, fontsize=12)
            plt.title(titel_page, fontsize=8,
                      fontweight='bold')
            title_leg = grup[0].cell_name + " " + self.std_type
            lgd = plt.legend(title=title_leg, loc='center left', bbox_to_anchor=(1, 0.5))
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            self.pdf_pages.savefig(fig, bbox_extra_artists=(lgd, sub), bbox_inches='tight')

    @staticmethod
    def get_name(grup, comp):
        x = []
        y = []
        for i in grup:
            a = i.split('_')
            if a[0] not in x:
                x.append(a[0])
            if a[1] not in y:
                y.append(a[1])

        uniquex = np.asarray(x)
        uniquey = np.asarray(y)
        # string = []
        #
        # for i in range(len(grup[comp].concentration)):
        #     n = grup[comp].concentration[i][:-1]
        #     name1 = grup[comp].name + " " + str(n[0])
        #
        #     string.append(name1)

        return uniquex, uniquey

    # #####
    # Plot medium and all the controls
    def plot_control_medium(self, grup, forplotting, mepos, positmedium, info1=0, info2=0):

        self.control_info = forplotting.control_info
        self.control_list = forplotting.control[positmedium]
        self.grup = grup
        self.data_plot = forplotting.data_plot
        self.std_plot = forplotting.std_plot
        self.mediumdata = forplotting.mediumdata[mepos]
        self.mediumstd = forplotting.mediumstd[mepos]
        self.readout = forplotting.readout
        self.readout_units = forplotting.readout_units

        y_axis_name = self.set_y_axis_name(info1)
        y_axis_min, y_axis_max = self.set_y_axis_limt(info2, np.max(self.data_plot), np.min(self.data_plot))

        # plot mehtod is only for cells with medium only and controls
        fig = plt.figure(figsize=(8, 8), dpi=80)
        ax = plt.axes()
        label_titel = "Cell lines: " + grup.cell_name + "\n Compound: " + grup.name + " & " + "Control"
        label_a = str(grup.name)
        titel_page = str(label_titel)
        post_a = int(grup.concentration[0][-1])
        plt.plot(self.elapse, self.mediumdata, label=label_a, color="black")
        if len(self.std_plot[post_a]) == 2:
            plt.fill_between(self.elapse, self.mediumstd[0],
                             self.mediumstd[-1],
                             alpha=0.5, edgecolor='lightgray', facecolor='lightgray')
        else:
            plt.fill_between(self.elapse, self.mediumdata - self.mediumstd,
                             self.mediumdata + self.mediumstd,
                             alpha=0.5, edgecolor='lightgray', facecolor='lightgray')
        color = iter(cm.tab20c(np.linspace(0, 1, len(self.control_list))))

        for j in range(len(self.control_list)):

            c = next(color)
            label_b = str(self.control_list[j][3]) + "_" + self.control_list[j][-2]  # control
            post_b = self.control_list[j][-1]

            plt.plot(self.elapse, self.data_plot[post_b], label=label_b, color=c)
            if len(self.std_plot[post_b]) == 2:
                plt.fill_between(self.elapse, self.std_plot[post_b][0],
                                 self.std_plot[post_b][-1],
                                 alpha=0.5, edgecolor='lightgray', facecolor='lightgray')
            else:
                plt.fill_between(self.elapse, self.data_plot[post_b] - self.std_plot[post_b],
                                 self.data_plot[post_b] + self.std_plot[post_b],
                                 alpha=0.5, edgecolor='lightgray', facecolor='lightgray')

        plt.xlabel("time (h)")
        plt.ylabel(y_axis_name)
        plt.ylim(y_axis_min, y_axis_max)
        sub = plt.suptitle(self.plot_titel, fontsize=12)
        plt.title(titel_page, fontsize=8,
                  fontweight='bold')
        title_leg = grup.cell_name + "  " + self.std_type
        lgd = plt.legend(title=title_leg, loc='center left', bbox_to_anchor=(1, 0.5))
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        self.pdf_pages.savefig(fig, bbox_extra_artists=(lgd, sub), bbox_inches='tight')

    def plot_medium(self, forplotting, grup, info1, info2):

        mediumdata = forplotting.mediumdata
        mediumstd = forplotting.mediumstd

        y_axis_name = self.set_y_axis_name(info1)
        y_axis_min, y_axis_max = self.set_y_axis_limt(info2, np.max(self.data_plot), np.min(self.data_plot))

        # plot mehtod is only for cells with medium only and controls
        fig = plt.figure(figsize=(8, 8), dpi=80)
        ax = plt.axes()
        name = [grup[i].name for i in range(len(grup))]
        label_titel = "Cell lines: " + grup[0].cell_name + "\n Compound: " + "&".join(name)

        titel_page = str(label_titel)

        color = iter(cm.tab20c(np.linspace(0, 1, len(grup))))

        for j in range(len(grup)):
            c = next(color)
            label_b = str(grup[j].name)   # control
            post_b = int(grup[j].concentration[0][-1])
            print(post_b)
            plt.plot(self.elapse, mediumdata[post_b], label=label_b, color=c)
            if len(mediumstd[post_b]) == 2:
                plt.fill_between(self.elapse, mediumstd[post_b][0],
                                 mediumstd[post_b][-1],
                                 alpha=0.5, edgecolor='lightgray', facecolor='lightgray')
            else:
                plt.fill_between(self.elapse, mediumdata[post_b] - mediumstd[post_b],
                                 mediumdata[post_b] + mediumstd[post_b],
                                 alpha=0.5, edgecolor='lightgray', facecolor='lightgray')

        plt.xlabel("time (h)")
        plt.ylabel(y_axis_name)
        plt.ylim(y_axis_min, y_axis_max)
        sub = plt.suptitle(self.plot_titel, fontsize=12)
        plt.title(titel_page, fontsize=8,
                  fontweight='bold')
        title_leg = grup[0].cell_name + "  " + self.std_type
        lgd = plt.legend(title=title_leg, loc='center left', bbox_to_anchor=(1, 0.5))
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        self.pdf_pages.savefig(fig, bbox_extra_artists=(lgd, sub), bbox_inches='tight')

    def plot_control_mediuminhibited(self, grup, forplotting, positioni, positionj, info1, info2):
        self.control_info = forplotting.control_info
        # self.grup = grup
        # self.data_plot = forplotting.data_plot
        # self.std_plot = forplotting.std_plot
        # self.readout = forplotting.readout
        # self.readout_units = forplotting.readout_units
        self.control_list = forplotting.control[positioni][positionj]

        self.grup = grup
        self.data_plot = forplotting.data_plot
        self.std_plot = forplotting.std_plot
        self.readout = forplotting.readout
        self.readout_units = forplotting.readout_units

        self.mediuminhfields = forplotting.mediuminhfiels[positioni][positionj][-1]
        if info1 == 1:
            self.mediumdata = forplotting.inhibitemediumdata[self.mediuminhfields]
            self.mediumstd = forplotting.inhibitemediumstd[self.mediuminhfields]

        elif info1 == 2:
            self.mediumdata = forplotting.translationmediumdata[self.mediuminhfields]
            self.mediumstd = forplotting.translationmediumdata[self.mediuminhfields]

        y_axis_name = self.set_y_axis_name(info1)
        y_axis_min, y_axis_max = self.set_y_axis_limt(info2, np.max(self.data_plot), np.min(self.data_plot))

        # plot mehtod is only for cells with medium only and controls
        fig = plt.figure(figsize=(8, 8), dpi=80)
        ax = plt.axes()
        label_titel = "Cell lines: " + grup.cell_name + "\n Compound: " + grup.name + " & " + "Control"
        label_a = str(grup.name)
        titel_page = str(label_titel)
        post_a = int(grup.concentration[0][-1])
        plt.plot(self.elapse, self.mediumdata, label=label_a, color="black")
        if len(self.std_plot[post_a]) == 2:
            plt.fill_between(self.elapse, self.mediumstd,
                             self.mediumstd[-1],
                             alpha=0.5, edgecolor='lightgray', facecolor='lightgray')
        else:
            plt.fill_between(self.elapse, self.mediumdata - self.mediumstd,
                             self.mediumdata + self.mediumstd,
                             alpha=0.5, edgecolor='lightgray', facecolor='lightgray')
        color = iter(cm.tab20c(np.linspace(0, 1, len(self.control_list))))
        c = next(color)

        label_b = str(self.control_list[3]) + "_" + self.control_list[-2]  # control
        post_b = self.control_list[-1]

        plt.plot(self.elapse, self.data_plot[post_b], label=label_b, color=c)
        if len(self.std_plot[post_b]) == 2:
            plt.fill_between(self.elapse, self.std_plot[post_b][0],
                             self.std_plot[post_b][-1],
                             alpha=0.5, edgecolor='lightgray', facecolor='lightgray')
        else:
            plt.fill_between(self.elapse, self.data_plot[post_b] - self.std_plot[post_b],
                             self.data_plot[post_b] + self.std_plot[post_b],
                             alpha=0.5, edgecolor='lightgray', facecolor='lightgray')

        plt.xlabel("time (h)")
        plt.ylabel(y_axis_name)
        plt.ylim(y_axis_min, y_axis_max)
        sub = plt.suptitle(self.plot_titel, fontsize=12)
        plt.title(titel_page, fontsize=8,
                  fontweight='bold')
        title_leg = grup.cell_name + "  " + self.std_type
        lgd = plt.legend(title=title_leg, loc='center left', bbox_to_anchor=(1, 0.5))
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        self.pdf_pages.savefig(fig, bbox_extra_artists=(lgd, sub), bbox_inches='tight')


