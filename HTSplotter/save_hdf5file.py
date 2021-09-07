import h5py
import os
from scipy import stats
from scipy.stats import norm
import numpy as np
from HTSplotter.hdf5functions import Hdf5functions
from copy import deepcopy

class Hdf5database:
    def __init__(self, experiment_type, branch,
                 file_name, header, elapse, date_info, date, data, std, std_info, medium, compound):

        self.experiment_name = file_name
        self.file_path = file_name
        self.file = None
        self.header = header
        self.branch = branch

        self.control_stru_level = 5
        self.experiment_type = experiment_type

        self.elapse = np.asarray(elapse)
        self.date_info = np.asarray(date_info)
        self.date = np.asarray(date)
        self.data_ini = data
        self.std_ini = std
        self.std_info = std_info

        self.control = []
        self.controlpath = []
        self.controlpathaverage = []
        self.controlpathdataaverage = []
        self.mediumpath = []
        self.controlpathdata = []
        self.medium = medium
        self.mediumword = None
        self.compound = compound

        # information obtained after HDF5 processing
        self.celline = []
        self.seeding = []
        self.condition = []
        self.conditionpath = []
        self.fields = []
        self.data = []
        self.inhibited = []
        self.normalized = []
        self.normalized_perc = []
        self.std = []
        self.std_inh = []
        self.confinterval = []
        # information for cells only:
        self.fieldsmedium = []
        self.fieldsmediuminhibited = []
        self.datamedium = []
        self.inhibitedmedium = []
        self.normalizedtranslation = []
        self.normalizedmedium = []
        self.normalizedtranslationmedium = []
        self.normalized_percmedium = []
        self.stdmedium = []
        self.std_inhmedium = []

        # combination information
        self.possiblecombination = []
        self.possiblecombinationsize = []
        # self.combination = []
        self.compoundalone = []
        self.hdfnorm = Hdf5functions()
        self.compoundalonearranged = []
        self.get_medium()
    def get_medium(self):
        if len(self.medium[0]) > 0:
            self.mediumword = self.medium[0].split("_")[0]

    # ##comomn functions!!!
    def open_file(self):
        self.file = h5py.File(self.file_path, "w")
        return

    def close_file(self):
        self.file.close()

    def structure_data(self):
        # create main information (elapse, complete date and date infor(first and last date)
        self.create_main_information()

        # create groups
        self.create_first_subgroup()

    def create_main_information(self):
        self.file.create_dataset("0Elapse", data=self.elapse, dtype="f")

        asciilistdate = [n.encode("ascii", "ignore") for n in self.date]
        self.file.create_dataset('000Date', (len(self.date), 1), 'S25', asciilistdate)

        asciilistdate_info = [n.encode("ascii", "ignore") for n in self.date_info]
        self.file.create_dataset('00Date_info', (len(self.date_info), 1), 'S25', asciilistdate_info)
        asciilistdate_info = [n.encode("ascii", "ignore") for n in self.std_info]
        self.file.create_dataset('00' + self.std_info[0], (len(self.std_info), 1), 'S25', asciilistdate_info)

    def create_first_subgroup(self):
        # build hdf5 according to the data structure, only 1 data contition
        if len(self.header[0][-1]) != 1:
            for j in self.header:
                try:
                    a = self.file.create_group("/".join(j[:-1]))
                    for k in j[-1]:
                        a.create_dataset('data_' + str(k), data=self.data_ini[:, k])
                except ValueError:
                    existpath = "/" + "/".join(j[:-1])
                    existpath2 = self.file[existpath]
                    for k in j[-1]:
                        existpath2.create_dataset('data_' + str(k), data=self.data_ini[:, k])
        else:
            for i in self.header:
                a = self.file.create_group("/".join(i[:-1]))
                a.create_dataset('data', data=self.data_ini[:, i[-1][0]])
                if len(self.std_ini) != 0:
                    a.create_dataset('std', data=self.std_ini[:, i[-1][0]])
                if len(self.std_ini) == 0:
                    a.create_dataset('std', data=np.zeros(self.data_ini[:, i[-1][0]].shape))
                if "Control" in a.name and "Condition" in a.name:
                    new = a.name.split("/")
                    compound = new[-2].split("_")
                    if len(compound) > 3:
                        com = "_".join(compound[-2:])
                    else:
                        com = compound[-1]
                    newname = "/".join(new[:4]) + "/" + com
                    self.controlpath.append(newname)
                    self.controlpathdata.append(self.data_ini[:, i[-1][0]])

    def comput_comb_average_std(self, f, level):
        if level > 0:
            try:
                keys = list(f.keys())
                for key in keys:
                    grp1 = f[key]
                    self.comput_comb_average_std(grp1, level - 1)

            except AttributeError:
                pass
        else:
            k = list(f.keys())
            for j in k:
                grp1 = f[j]
                n = list(grp1.keys())
                self.std_info = ['std', '95%CI']
                for c in range(len(n)):
                    grp2 = grp1[n[c]]
                    new = list(grp2.keys())
                    temp = [[] for i in range(0, len(new))]
                    for o in range(len(new)):
                        gr = grp2[new[o]]
                        temp[o] = gr[:]
                    df = len(temp) - 1
                    alpha = 1 - 0.95
                    grp2.create_dataset('data', data=np.mean(temp, axis=0))
                    grp2.create_dataset('std', data=np.std(temp, axis=0))
                    grp2.create_dataset('95%CI', data=stats.t.interval(1 - alpha, df, loc=np.mean(temp, axis=0),
                                                                       scale=np.std(temp, axis=0) / np.sqrt(
                                                                           len(temp))))
                    if "Control" in j:
                        self.controlpath.append(grp1.name)
                        self.controlpathdata.append(grp2['data'][:])

    def normalizeseveralcontrol(self):
        cell = list(self.file.keys())
        # skip always 0Elapse; 000Date and 00Date_info
        level = 0
        for k in cell[4:]:
            self.hdfnorm.get_normalizeseveralcontrol(self.file, level, 0, k, self.celline, self.seeding,
                                                     self.condition, self.controlpath, self.controlpathdata,
                                                     self.medium, self.mediumword)

    def normalizeonecontrol(self):
        cell = list(self.file.keys())
        # cellpath, seedingpath, conditionpath, compound = self.hdf5walk(self.file, cell)
        # skip always 0Elapse; 000Date and 00Date_info
        for k in cell[4:]:
            self.hdfnorm.get_normalizeonecontrol(self.file, 0, k, self.celline, self.seeding,
                                                 self.condition, self.conditionpath, self.controlpath,
                                                 self.controlpathdata, self.medium, self.mediumword)

    def get_medium_control(self, f, level):
        if level > 0:
            try:
                keys = list(f.keys())
                for key in keys:
                    grp1 = f[key]
                    self.get_medium_control(grp1, level - 1)
            except AttributeError:
                pass
        else:
            self.hdfnorm.medium_control(f, 0, self.mediumword, self.control, self.compoundalone,
                                        self.fieldsmedium, self.fieldsmediuminhibited)

    def get_combination(self, f, level):
        if level > 0:
            try:
                keys = list(f.keys())
                for key in keys:
                    grp1 = f[key]
                    self.get_combination(grp1, level-1)
            except AttributeError:
                pass
        else:
            self.hdfnorm.get_combination(f, self.compoundalone, self.possiblecombination,
                                         self.possiblecombinationsize)

    def get_compoundalone(self, f, level):
        if level > 0:
            try:
                keys = list(f.keys())
                for key in keys:
                    grp1 = f[key]
                    self.get_compoundalone(grp1, level - 1)
            except AttributeError:
                pass
        else:
            alone = []
            compound = list(f.keys())
            for i in compound:
                temp = []
                if "Control" not in i:
                    alone.append(i)
            self.compoundalone.append(alone)

    def add_predictedbiscore(self, groupname, biscore, predicted):
        self.file = h5py.File(self.file_path, "r+")
        self.hdfnorm.add_comboinformation(groupname, biscore, predicted, self.file,
                                          self.possiblecombination, self.conditionpath,
                                          self.possiblecombinationsize)
        self.close_file()

class Compoundscreenonecontrol(Hdf5database):
    def __init__(self, experiment_type, branch, file_name, header, elapse, date_info, date, data, std,
                 std_info, medium, compound):
        super().__init__(experiment_type, branch, file_name, header, elapse, date_info, date, data, std,
                         std_info, medium, compound)

        self.open_file()
        self.structure_data()
        if len(self.header[0][-1]) == 1 and len(self.std_ini) == 0:
            self.std_info = ["No STD"]
        elif len(self.header[0][-1]) > 1:
            self.std_info = ["STD computed by HTSplotter"]
            self.comput_comb_average_std(self.file, 3)
        self.normalizeonecontrol()
        self.hdfnorm.get_fieldsonecontrolmain(self.file, self.control_stru_level, 0, self.fields,
                                              self.control, self.data,
                                              self.inhibited, self.normalized_perc, self.normalized,
                                              self.std_inh, self.std, self.normalizedtranslation,
                                              self.datamedium, self.inhibitedmedium, self.normalized_percmedium,
                                              self.normalizedmedium, self.std_inhmedium, self.stdmedium,
                                              self.confinterval, self.normalizedtranslationmedium)

        self.get_compoundalone(self.file, 3)

        self.close_file()


class Compoundcombination(Hdf5database):
    def __init__(self, experiment_type, branch, file_name, header, elapse, date_info, date, data, std,
                 std_info, medium, compound):
        super().__init__(experiment_type, branch, file_name, header, elapse, date_info, date, data, std,
                         std_info, medium, compound)

        self.open_file()
        self.structure_data()
        if len(self.header[0][-1]) == 1 and len(self.std_ini) == 0:
            self.std_info = ["No STD"]

        elif len(self.header[0][-1]) > 1:
            self.std_info = ["STD computed by HTSplotter"]
            self.comput_comb_average_std(self.file, 3)

        self.normalizeonecontrol()

        self.hdfnorm.get_fieldsonecontrolmain(self.file, self.control_stru_level, 0, self.fields,
                                              self.control, self.data,
                                              self.inhibited, self.normalized_perc, self.normalized,
                                              self.std_inh, self.std, self.normalizedtranslation,
                                              self.datamedium, self.inhibitedmedium, self.normalized_percmedium,
                                              self.normalizedmedium, self.std_inhmedium, self.stdmedium,
                                              self.confinterval, self.normalizedtranslationmedium)
        self.get_combination(self.file, 3)
        self.close_file()


class Compoundscreen(Hdf5database):
    def __init__(self, experiment_type, branch, file_name, header, elapse, date_info, date, data, std,
                 std_info, medium, compound):
        super().__init__(experiment_type, branch, file_name, header, elapse, date_info, date, data, std,
                         std_info, medium, compound)

        self.open_file()
        self.structure_data()
        if len(self.header[0][-1]) == 1 and len(self.std_ini) == 0:
            self.std_info = ["No STD"]

        elif len(self.header[0][-1]) > 1:
            self.std_info = ["STD computed by HTSplotter"]
            self.comput_comb_average_std(self.file, 3)

        self.normalizeseveralcontrol()
        self.hdfnorm.get_fieldsseveralcontrolmain(self.file, self.control_stru_level, 0, self.fields, self.data,
                                                  self.inhibited,
                                                  self.normalized_perc, self.normalized, self.std_inh, self.std,
                                                  self.normalizedtranslation, self.datamedium, self.inhibitedmedium,
                                                  self.normalized_percmedium, self.normalizedmedium, self.std_inhmedium,
                                                  self.stdmedium, self.confinterval, self.normalizedtranslationmedium)
        self.get_medium_control(self.file, 3)
        self.close_file()


class Geneticperturbagem(Hdf5database):
    def __init__(self, experiment_type, branch, file_name, header, elapse, date_info, date, data, std,
                 std_info, medium, compound):
        super().__init__(experiment_type, branch, file_name, header, elapse, date_info, date, data, std,
                         std_info, medium, compound)
        # this is the only experimental condition that accepts more than 1 control
        self.controlnamenew = []
        self.controlnamenewdata = []
        for i in self.header:
            if "Control" in i[3] or "control" in i[3]:
                name = i[3].split("_")
                i[3] = name[-1]
                self.controlnamenew.append(name[-1])
                self.controlnamenewdata.append(i[-1])
        self.controlnamenew = "_".join(self.controlnamenew)

        self.open_file()
        self.structure_dataarray()

        if len(self.header[0][-1]) == 1 and len(self.std_ini) == 0:
            self.std_info = ["No STD"]

        elif len(self.header[0][-1]) > 1:
            self.std_info = ["STD computed by HTSplotter"]
            self.comput_comb_average_std(self.file, 3)


        self.normalizeonecontrol()
        self.hdfnorm.get_fieldsonecontrolmain(self.file, self.control_stru_level, 0, self.fields,
                                              self.control, self.data,
                                              self.inhibited, self.normalized_perc, self.normalized,
                                              self.std_inh, self.std, self.normalizedtranslation,
                                              self.datamedium, self.inhibitedmedium, self.normalized_percmedium,
                                              self.normalizedmedium, self.std_inhmedium, self.stdmedium,
                                              self.confinterval, self.normalizedtranslationmedium)

        self.get_compoundalone(self.file, 3)

        self.close_file()

    def structure_dataarray(self):
        # create main information (elapse, complete date and date infor(first and last date)
        self.create_main_information()

        # create groups
        self.create_first_subgroup_array()

    def create_first_subgroup_array(self):
        # build hdf5 according to the data structure, only 1 data contition
        if len(self.header[0][-1]) != 1:
            for j in self.header:
                try:
                    a = self.file.create_group("/".join(j[:-1]))
                    for k in j[-1]:
                        a.create_dataset('data_' + str(k), data=self.data_ini[:, k])
                except ValueError:
                    existpath = "/" + "/".join(j[:-1])
                    existpath2 = self.file[existpath]
                    for k in j[-1]:
                        existpath2.create_dataset('data_' + str(k), data=self.data_ini[:, k])
            controlmainpath = self.file.create_group("/" + "/".join(j[:-3]) + "/Control" + "/" + self.controlnamenew)
            for i in self.controlnamenewdata:
                for h in i:
                    controlmainpath.create_dataset('data_' + str(h), data=self.data_ini[:, h])

        else:
            for i in self.header:
                a = self.file.create_group("/".join(i[:-1]))
                a.create_dataset('data', data=self.data_ini[:, i[-1]])
                a.create_dataset('std', data=self.std_ini[:, i[-1]])


class Geneticchemicalperturbagem(Hdf5database):
    def __init__(self, experiment_type, branch, file_name, header, elapse, date_info, date, data, std,
                 std_info, medium, compound):
        super().__init__(experiment_type, branch, file_name, header, elapse, date_info, date, data, std,
                         std_info, medium, compound)
        self.branch = []
        self.open_file()
        self.structure_data()
        if len(self.header[0][-1]) == 1 and len(self.std_ini) == 0:
            self.std_info = ["No STD"]

        elif len(self.header[0][-1]) > 1:
            self.std_info = ["STD computed by HTSplotter"]
            self.comput_comb_average_std(self.file, 3)


        self.normalizeseveralcontrol()

        self.hdfnorm.get_fieldsseveralcontrolmain(self.file, self.control_stru_level, 0, self.fields, self.data, self.inhibited,
                                                  self.normalized_perc, self.normalized, self.std_inh, self.std,
                                                  self.normalizedtranslation, self.datamedium, self.inhibitedmedium,
                                                  self.normalized_percmedium, self.normalizedmedium, self.std_inhmedium,
                                                  self.stdmedium, self.confinterval, self.normalizedtranslationmedium)

        self.hdfnorm.get_genetcombmain(self.file, 2, self.mediumword, self.possiblecombination,
                                       self.possiblecombinationsize)

        self.get_medium_control(self.file, 3)
        self.compoundalonearranged = self.hdfnorm.get_branch(0, self.file, self.branch,
                                                             self.compoundalone)
        self.compoundalone = self.compoundalonearranged

        self.close_file()
