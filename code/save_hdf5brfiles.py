import h5py
import os
import numpy as np
from scipy import stats
from hdf5functions import Hdf5functions
from scipy.stats import norm

file_name_each = []
branch_BR = []
cell_line_BR = []
BiologicalReplicates = []
seeding_BR = []
condition_BR = []
compound_BR = []
concentration = []
data_BR = []
std_BR = []
compound_control_list = []


class BRHdf5database:

    def __init__(self, count, file_name_br, experiment_type, branch,
                 file_name, header, elapse, date_info, date, data, std, medium):

        self.file_path = file_name_br
        self.experiment_name = file_name
        self.individual_file_name = file_name
        self.count_br = count
        # self.file_path = os.path.join(folder, self.experiment_name)
        self.header = header
        self.branch = branch
        self.elapse = np.asarray(elapse)
        self.date_info = np.asarray(date_info)
        self.date = np.asarray(date)

        self.medium = medium
        self.mediumword = medium[0].split("_")[0]

        self.f = None
        self.main_grup = None

        # self.control_information = []
        # self.compound_information = []
        # self.control_name = None  # control name for Combination
        self.control_stru_level = 3
        self.experiment_type = experiment_type

        # ### Important information to generate the Data_base
        self.data_ini = data
        self.std_ini = std

        self.brpath = None
        self.eachpath = []  # path to compute technical replicates avarage

        # information that can be onbatin usin the method get_normalize
        self.fields = []
        self.data = []
        self.normalized = []
        self.normalized_perc = []
        self.inhibited = []
        self.std = []
        self.std_inh = []
        self.std_perc = []
        self.confinterval = []
        ##################

        # ### information for BR
        self.controlpath = []
        self.controlpathdata = []
        self.control = []

        # just for genetic perturbagem
        self.controlnamenew = None
        # information for cells only:
        self.fieldsmedium = []
        self.fieldsmediuminhibited = []
        self.datamedium = []
        self.inhibitedmedium = []
        self.normalizedmedium = []
        self.normalizedtranslation = []
        self.normalizedtranslationmedium = []
        self.normalized_percmedium = []
        self.stdmedium = []
        self.std_inhmedium = []

        # combination information
        self.possiblecombination = []
        self.possiblecombinationsize = []
        self.compoundalone = []
        self.compoundalonearranged = []

        # ##This part is given by the method synersm_groups

        self.celline = []  # in case of having more than 1 we need to know
        self.seeding = []  # in case the seeding is different
        self.condition = []
        self.conditionpath = []# in case condition is more than 1



    # ##comomn functions!!!
    def open_file(self):

        if self.count_br == 0:
            self.f = h5py.File(self.file_path,  "w")
        else:
            self.f = h5py.File(self.file_path, "a")

        return

    def read_file(self):
        self.f = h5py.File(self.file_path, "r+")

    def close_file(self):
        self.f.close()

    def structure_data(self):
        self.main_grup = self.f.create_group(self.individual_file_name)

        # create main information (elapse, complete date and date infor(first and last date)
        self.create_main_information()

        # create groups
        self.create_first_subgroup()

    def create_main_information(self):
        self.main_grup.create_dataset("0Elapse", data=self.elapse, dtype="f")

        asciilistdate = [n.encode("ascii", "ignore") for n in self.date]
        self.main_grup.create_dataset('000Date', (len(self.date), 1), 'S25', asciilistdate)

        asciilistdate_info = [n.encode("ascii", "ignore") for n in self.date_info]
        self.main_grup.create_dataset('00Date_info', (len(self.date_info), 1), 'S25', asciilistdate_info)

    def create_first_subgroup(self):
        if len(self.header[0][-1]) != 1:
            for j in self.header:
                try:
                    a = self.main_grup.create_group("/".join(j[:-1]))
                    self.eachpath.append(a)

                    for k in j[-1]:
                        a.create_dataset('data_' + str(k), data=self.data_ini[:, k])
                except ValueError:
                    existpath = "/".join(j[:-1])
                    existpath2 = self.main_grup[existpath]
                    for k in j[-1]:
                        existpath2.create_dataset('data_' + str(k), data=self.data_ini[:, k])
        else:
            k = 0
            for field in self.header:
                grp = self.main_grup
                for i in range(len(field)):
                    keys = list(grp.keys())
                    if field[i] in keys:
                        grp = grp[field[i]]
                    else:
                        break
                if i + 1 <= len(field):
                    for j in range(i, len(field) - 1):
                        grp.create_group(field[j])
                        grp = grp[field[j]]

                grp.create_dataset('data', data=self.data_ini[:, k])
                if len(self.std_ini) != 0:
                    grp.create_dataset('std', data=self.std_ini[:, k])
                if len(self.std_ini) == 0:
                    grp.create_dataset('std', data=np.zeros(self.data_ini[:, k].shape))

                k += 1

    def get_pathways(self):
        # cycle for all compounds:
        path_com_list = []
        compound_control_list.sort()

        for i in range(len(compound_control_list)):
            tem_list = []
            name_list = []
            mean_list = []
            std_list = []
            confinterval = []
            self.recursive_function(self.f, self.control_stru_level, i, tem_list, name_list)
        path_com_list, mean_list, std_list, confinterval = self.get_compound_concentration_eachbr(tem_list, name_list,
                                                                                    path_com_list, mean_list,
                                                                                    std_list, confinterval)
        self.save_br_info(path_com_list, mean_list, std_list, confinterval)

    def get_compound_concentration_eachbr(self, tem_list, name_list, path_com_list, mean_list,
                                          std_list, confinterval):
        name_p = None
        aux_list = None
        aux_list2 = None

        cicle_concentration = int(len(tem_list) / len(cell_line_BR))

        for k in range(int(len(name_list) / len(cell_line_BR))):
            aux_list = np.zeros(len(tem_list[k]))
            aux_list2 = np.zeros(len(tem_list[k]))

            for t in range(len(cell_line_BR)):
                name_p = name_list[k + t * cicle_concentration]
                aux_list += tem_list[k + t * cicle_concentration]  # each BR

                aux_list2 += tem_list[k + t * cicle_concentration] ** 2  # each BR Sum of squares

            # compute average and STD

            self.get_average_std(mean_list, std_list, aux_list, aux_list2, confinterval, cell_line_BR)

            path = self.get_path_br(name_p)
            path_com_list.append(path)

            # Save mean and STD
        return path_com_list, mean_list, std_list, confinterval

    def recursive_function(self, f, level, i, tem_list, name_list):
        if level > 0:
            try:
                keys = list(f.keys())
                for key in keys:
                    grp1 = f[key]
                    self.recursive_function(grp1, level - 1, i, tem_list, name_list)

            except AttributeError:
                pass
        else:
            self.processing(f, i, tem_list, name_list)

    @staticmethod
    def get_path_br(name_path):
        path = []
        name_path_split = name_path.split("/")
        name = name_path_split[1]
        subpath = "/".join(name_path_split[2:])
        name2 = name.split(" ")
        add = "_BR"
        path1 = "/" + add + "/" + subpath  # name2[0] + "_" +
        path.append(path1)

        return path

    @staticmethod
    def get_average_std(mean_list, std_list, aux_list, aux_list2, confinterval, cell_line_BR):

        mean_list.append(aux_list / len(cell_line_BR))

        std_list.append(((aux_list2 / len(cell_line_BR) - (mean_list[-1] ** 2)))**0.5)  # STD sample
        std_list[-1][np.isnan(std_list[-1])] = 0
        # calculate 95 % CI
        df = len(mean_list) - 1
        alpha = 1 - 0.95
        inter = (stats.t.interval(alpha=0.95, df=len(aux_list)-1, loc=mean_list,
                                  scale=stats.sem(aux_list)))
        confinterval.append(np.squeeze([entry for entry in inter]))

        return mean_list, std_list, confinterval

    @staticmethod
    def processing(f, i, tem_list, name_list):
        k = list(f.keys())
        for j in k:
            a = list(f[j].keys())
            for cada in a:
                nom = j + "/" + cada
                grp1 = f[nom]
                b = list(grp1.keys())
                # through all concentrations
                for u in b:
                    nom2 = nom + "/" + u
                    grp2 = f[nom2]
                    name_list.append(grp2.name)
                    tem_list.append(f[nom2]['data'][:])

    def save_br_info(self, name_path, mean_list, std_list, confinterval):
        g = self.f

        for i in range(len(name_path)):
            grp2 = g.create_group(name_path[i][0])
            grp2.create_dataset('mean', data=mean_list[i])
            grp2.create_dataset('std_BR', data=std_list[i])
            grp2.create_dataset('95%CI_BR', data=confinterval[i].T)
            if "Control" in name_path[i][0]:
                path = name_path[i][0].split("/")
                if "Condition" not in path[-3]:
                    controlname = path[-2].split("_")[-2:]
                    path[-2] = "_".join(controlname)
                    self.controlpath.append("/".join(path[:-1]))
                    self.controlpathdata.append(grp2['mean'][:])
                else:
                    controlname = path[-2].split("_")[-1]
                    path[-2] = controlname
                    self.controlpath.append("/".join(path[:-1]))
                    self.controlpathdata.append(grp2['mean'][:])

    # #####
    # Compute avarage for single data set
    def computaverage(self):
        for i in self.eachpath:
            k = list(i.keys())
            temp = [[] for r in range(0, len(k))]
            for ea in range(len(k)):
                grp = i[k[ea]]
                temp[ea] = grp[:]
            i.create_dataset('data', data=np.mean(temp, axis=0))

    ## Mean normalization only for the BR

    def normalizeseveralcontrol(self):
        experiments = list(self.f.keys())
        self.hdfnorm = Hdf5functions()
        # skip always 0Elapse; 000Date and 00Date_info
        level = 0
        for exp in experiments:
            if "_BR" == exp:
                self.brpath = self.f[exp]
                cell = list(self.f[exp].keys())
                for k in cell:
                    self.hdfnorm.get_normalizeseveralcontrol(self.f, level, exp, k, self.celline, self.seeding,
                                                             self.condition, self.controlpath, self.controlpathdata,
                                                             self.medium, self.mediumword)
                    # level += 1
                    # cellpath = self.f[exp][k]
                    # self.celline.append(k)
                    # seeding = list(cellpath.keys())
                    # for j in seeding:
                    #     self.seeding.append(j)
                    #     level += 1
                    #     seedingpath = cellpath[j]
                    #     condition = list(seedingpath.keys())
                    #     self.condition.append(condition)
                    #     for g in condition:
                    #         level += 1
                    #         conditionpath = seedingpath[g]
                    #         compound = list(conditionpath.keys())
                    #         for c in compound:
                    #             compoundpath = conditionpath[c]
                    #             for i in range(len(self.controlpath)):
                    #                 if compoundpath.name == self.controlpath[i]:
                    #                     dosage = list(compoundpath.keys())
                    #                     for d in dosage:
                    #                         grup = compoundpath[d]
                    #                         controlname = self.controlpath[i].split("/")
                    #                         self.hdf5normalization(grup, controlname[-1], self.controlpathdata[i])
                    #                 elif "Condition" not in compoundpath.name:
                    #                     splitpath = self.controlpath[i].split("/")
                    #                     pathmain = "/".join(splitpath[:-2])
                    #                     pathcontrolcompound = splitpath[-1]
                    #                     if pathmain in compoundpath.name and pathcontrolcompound in compoundpath.name:
                    #                         dosage = list(compoundpath.keys())
                    #                         for d in dosage:
                    #                             grup = compoundpath[d]
                    #                             controlname = self.controlpath[i].split("/")
                    #                             self.hdf5normalization(grup, controlname[-1], self.controlpathdata[i])
                    #                     elif self.mediumword in c:
                    #                         singledosage = list(compoundpath.keys())
                    #                         grup = compoundpath[singledosage[0]]
                    #                         new4 = self.controlpath[i].split("/")
                    #                         new3 = "/".join(new4[:-2])
                    #                         if new3 in conditionpath.name:
                    #                             self.hdf5normalization(grup, new4[-1], self.controlpathdata[i])
                    #                 if c in self.medium:
                    #                     singledosage = list(compoundpath.keys())
                    #                     grup = compoundpath[singledosage[0]]
                    #                     new = self.controlpath[i].split("/")
                    #                     new2 = "/".join(new[:-1])
                    #                     if new2 in conditionpath.name:
                    #                         self.hdf5normalization(grup, new[-1], self.controlpathdata[i])
                    #                 if "Control" in c:
                    #                     singledosage = list(compoundpath.keys())
                    #                     grup = compoundpath[singledosage[0]]
                    #                     new = compoundpath.name.split("/")
                    #                     new2 = new[-1].split("_")
                    #                     if len(new2) > 3:
                    #                         new3 = "_".join(new2[-2:])
                    #                         new4 = "/".join(new[:-1])
                    #                         complete = new4 + "/" + new3
                    #                         if complete == self.controlpath[i]:
                    #                             self.hdf5normalization(grup, new[-1], self.controlpathdata[i])
                    #                     else:
                    #                         new3 = "_".join(new2[-1:])
                    #                         new4 = "/".join(new[:-1])
                    #                         complete = new4 + "/" + new3
                    #                         if complete == self.controlpath[i]:
                    #                             self.hdf5normalization(grup, new[-1], self.controlpathdata[i])

    def normalizeonecontrol(self):
        experiments = list(self.f.keys())
        self.hdfnorm = Hdf5functions()
        # skip always 0Elapse; 000Date and 00Date_info
        level = 0
        for exp in experiments:
            if "_BR" == exp:
                self.brpath = self.f[exp]
                cell = list(self.f[exp].keys())
                for k in cell:
                    self.hdfnorm.get_normalizeonecontrol(self.f, exp, k, self.celline, self.seeding,
                                                         self.condition, self.conditionpath, self.controlpath,
                                                         self.controlpathdata, self.medium, self.mediumword)
                    # cellpath = self.f[exp][k]
                    # seeding = list(cellpath.keys())
                    # self.celline.append(k)
                    # for j in seeding:
                    #     seedingpath = cellpath[j]
                    #     condition = list(seedingpath.keys())
                    #     self.seeding.append(j)
                    #     for g in condition:
                    #         conditionpath = seedingpath[g]
                    #         compound = list(conditionpath.keys())
                    #         self.condition.append(g)
                    #         self.conditionpath.append(conditionpath.name)
                    #         for c in compound:
                    #             compoundpath = conditionpath[c]
                    #             dosage = list(compoundpath.keys())
                    #             for i in range(len(self.controlpath)):
                    #                 for d in dosage:
                    #                     grup = compoundpath[d]
                    #                     controlname = self.controlpath[i].split("/")
                    #                     controlmainpath = "/".join(controlname[:-1])
                    #                     if controlmainpath == conditionpath.name:
                    #                         self.hdf5normalization(grup, controlname[-1], self.controlpathdata[i])

    def get_combination(self, f, level):
        self.hdfnorm = Hdf5functions()
        if level > 0:
            try:
                keys = list(f.keys())
                for key in keys:
                    grp1 = f[key]
                    self.get_combination(grp1, level - 1)
            except AttributeError:
                pass
        else:
            self.hdfnorm.get_combination(f, self.compoundalone, self.possiblecombination,
                                         self.possiblecombinationsize)

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
            self.hdfnorm.medium_control(f, 1, self.mediumword, self.control, self.compoundalone,
                                        self.fieldsmedium, self.fieldsmediuminhibited)
            # compound = list(f.keys())
            # alone = []
            # mediumne = []
            # inhmedium = []
            # control = []
            # correcpositi = 0
            # correcpositiinhib = 0
            # if len(self.fieldsmedium) > 0:
            #     posit = self.fieldsmedium[-1][-1][-1]
            #     correcpositi = 1
            # if len(self.fieldsmediuminhibited) > 0:
            #     positinh = self.fieldsmediuminhibited[-1][-1][-1]
            #     correcpositiinhib = 1
            # # if f.name
            # for i in compound:
            #     if "Control" in i and "Condition" in f.name:
            #         cont = f[i]
            #         k = list(cont.keys())
            #         cont2 = cont[k[0]]
            #         control.append(cont2.name.split("/")[2:])
            #     elif self.mediumword in i and '_BR' in f.name:
            #         medium = f[i]
            #         k = list(medium.keys())
            #         a = medium[k[0]]
            #         cond = list(a.keys())
            #         makefields = a.name.split("/")
            #         listtemp = makefields[2:]
            #         temmed = []
            #         teminh = []
            #         for c in cond:
            #             lista2 = listtemp.copy()
            #             lista3 = listtemp.copy()
            #             if "_" in c and 'inhibited' in c:
            #                 lista3[3] += "_" + c.split("_")[-1]
            #                 if lista3[3] not in teminh:
            #                     if correcpositiinhib == 1:
            #                         posit += 1
            #                         lista3.append(positinh)
            #                     else:
            #                         lista3.append(len(inhmedium))
            #                     teminh.append(lista3[3])
            #                     inhmedium.append(lista3)
            #             if lista2[3] not in temmed:
            #                 temmed.append(lista2[3])
            #                 if correcpositi == 1:
            #                     posit += 1
            #                     lista2.append(posit)
            #                 else:
            #                     lista2.append(len(mediumne))
            #                 mediumne.append(lista2)
            #     elif i not in alone:
            #         alone.append(i)
            # if len(control) != 0:
            #     self.control.append(control)
            # self.compoundalone.append(alone)
            # if len(mediumne) > 0:
            #     self.fieldsmedium.append(mediumne)
            # if len(inhmedium) > 0:
            #     self.fieldsmediuminhibited.append(inhmedium)

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

    # ## ADD predicted and bliss score to the HDF5 file, method called by main

    def add_predictedbiscore(self, groupname, biscore, predicted):
        self.file = h5py.File(self.file_path, "r+")
        br = list(self.file.keys())
        for i in br:
            if "_BR" == i:
                self.brpath = self.file[i]

        self.hdfnorm.add_comboinformation(groupname, biscore, predicted, self.file,
                                          self.possiblecombination, self.conditionpath,
                                          self.possiblecombinationsize)
        self.file.close()

class Individualcombinationstructure(BRHdf5database):
    def __init__(self, count, file_name_br, experiment_type, branch, file_name,
                 header, elapse, date_info, date, data, std, medium):
        super().__init__(count, file_name_br, experiment_type, branch, file_name,
                         header, elapse, date_info, date, data, std, medium)
        # ## once it is none the type of experiment we can save the raw information on Data base

        cell_line_BR.append([])
        seeding_BR.append([])
        condition_BR.append([])
        compound_BR.append([])
        concentration.append([])
        branch_BR.append([])
        data_BR.append([])
        std_BR.append([])
        file_name_each.append([])

        cell_line_BR[count].append(self.header[0][0])
        seeding_BR[count].append(self.header[0][1])
        condition_BR[count].append(self.header[0][2])
        branch_BR[count].append(self.branch)
        data_BR[count].append(self.data)
        std_BR[count].append(self.std)
        file_name_each[count].append(self.experiment_name)
        for i in self.header:
            compound_BR[count].append(i[3])
            concentration[count].append(i[4])
            if i[3] not in compound_control_list:
                compound_control_list.append(i[3])

        count += 1
        self.coun = count
        self.open_file()
        self.structure_data()
        if len(self.header[0][-1]) > 1:
            self.computaverage()
        self.close_file()


class Individualgeneticperturbagen(BRHdf5database):
    def __init__(self, count, file_name_br, experiment_type, branch, file_name,
                 header, elapse, date_info, date, data, std, medium):
        super().__init__(count, file_name_br, experiment_type, branch, file_name,
                         header, elapse, date_info, date, data, std, medium)
        # ## once it is none the type of experiment we can save the raw information on Data base

        cell_line_BR.append([])
        seeding_BR.append([])
        condition_BR.append([])
        compound_BR.append([])
        concentration.append([])
        branch_BR.append([])
        data_BR.append([])
        std_BR.append([])
        file_name_each.append([])

        cell_line_BR[count].append(self.header[0][0])
        seeding_BR[count].append(self.header[0][1])
        condition_BR[count].append(self.header[0][2])
        branch_BR[count].append(self.branch)
        data_BR[count].append(self.data)
        std_BR[count].append(self.std)
        file_name_each[count].append(self.experiment_name)
        self.controlnamenew = []
        self.controlnamenewdata = []
        for i in self.header:
            compound_BR[count].append(i[3])
            concentration[count].append(i[4])
            if "Control" in i[3] or "control" in i[3]:
                name = i[3].split("_")
                i[3] = name[-1]
                self.controlnamenew.append(name[-1])
                self.controlnamenewdata.append(i[-1])
                compound_control_list.append(self.controlnamenew)
        self.controlnamenew = "_".join(self.controlnamenew)

        count += 1

        self.coun = count
        self.open_file()
        self.structure_dataarray()
        if len(self.header[0][-1]) > 1:
            self.computaverage()
        self.close_file()

    def structure_dataarray(self):
        # create main information (elapse, complete date and date infor(first and last date)
        self.main_grup = self.f.create_group(self.individual_file_name)
        self.create_main_information()

        # create groups
        self.create_first_arraysubgroup_array()

    def create_first_arraysubgroup_array(self):
        # build hdf5 according to the data structure, only 1 data contition
        if len(self.header[0][-1]) != 1:
            for j in self.header:
                try:
                    a = self.main_grup.create_group("/".join(j[:-1]))
                    self.eachpath.append(a)
                    for k in j[-1]:
                        a.create_dataset('data_' + str(k), data=self.data_ini[:, k])
                except ValueError:
                    existpath = "/".join(j[:-1])
                    existpath2 = self.main_grup[existpath]
                    for k in j[-1]:
                        existpath2.create_dataset('data_' + str(k), data=self.data_ini[:, k])
            controlmainpath = self.main_grup.create_group("/".join(j[:-3]) + "/Control" + "/" + self.controlnamenew)
            self.eachpath.append(controlmainpath)
            for i in self.controlnamenewdata:
                for h in i:
                    controlmainpath.create_dataset('data_' + str(h), data=self.data_ini[:, h])
        else:
            for i in self.header:
                a = self.file.create_group("/".join(i[:-1]))
                a.create_dataset('data', data=self.data_ini[:, i[-1]])
                a.create_dataset('std', data=self.std_ini[:, i[-1]])


class BRcompoundscreenseveralcontrol(BRHdf5database):
    def __init__(self, count, file_name_br, experiment_type, branch, file_name,
                 header, elapse, date_info, date, data, std, medium):
        super().__init__(count, file_name_br, experiment_type, branch, file_name,
                         header, elapse, date_info, date, data, std, medium)

        # self.control_name = []
        self.std_info = ["STD_BR"]
        self.control_data_single = None

        self.open_file()
        # compute mean and STD
        self.get_pathways()
        # ##
        # normalize only the mean
        self.normalizeseveralcontrol()
        self.hdfnorm.get_fieldsseveralcontrolmain(self.brpath, 5, 1, self.fields, self.data, self.inhibited,
                                                  self.normalized_perc, self.normalized, self.std_inh, self.std,
                                                  self.normalizedtranslation, self.datamedium, self.inhibitedmedium,
                                                  self.normalized_percmedium, self.normalizedmedium, self.std_inhmedium,
                                                  self.stdmedium, self.confinterval, self.normalizedtranslationmedium)
        self.get_medium_control(self.brpath, 3)
        self.close_file()


class BRcompoundscreenonecontrol(BRHdf5database):
    def __init__(self, count, file_name_br, experiment_type, branch, file_name,
                 header, elapse, date_info, date, data, std, medium):
        super().__init__(count, file_name_br, experiment_type, branch, file_name,
                         header, elapse, date_info, date, data, std, medium)

        # self.control_name = []
        self.std_info = ["STD_BR"]
        self.control_data_single = None

        self.open_file()
        # compute mean and STD
        self.get_pathways()
        # ##
        # normalize only the mean
        self.normalizeonecontrol()
        self.hdfnorm.get_fieldsonecontrolmain(self.brpath, 5, 1, self.fields, self.control, self.data,
                                              self.inhibited, self.normalized_perc, self.normalized,
                                              self.std_inh, self.std, self.normalizedtranslation,
                                              self.datamedium, self.inhibitedmedium, self.normalized_percmedium,
                                              self.normalizedmedium, self.std_inhmedium, self.stdmedium,
                                              self.confinterval, self.normalizedtranslationmedium)
        self.get_compoundalone(self.brpath, 3)

        self.close_file()


class BRcombinationstructure(BRHdf5database):
    def __init__(self, count, file_name_br, experiment_type, branch, file_name,
                 header, elapse, date_info, date, data, std, medium):
        super().__init__(count, file_name_br, experiment_type, branch, file_name,
                         header, elapse, date_info, date, data, std, medium)

        self.celline = []
        self.seeding = []
        self.condition = []
        # self.control_name = []
        self.std_info = ["STD_BR"]
        self.control_data_single = None

        self.open_file()
        # compute mean and STD
        self.get_pathways()
        # ##
        # normalize only the mean
        self.normalizeonecontrol()
        # self.get_fieldsonecontrol(self.brpath, 5)
        self.hdfnorm.get_fieldsonecontrolmain(self.brpath, 5, 1, self.fields, self.control, self.data, self.inhibited,
                                              self.normalized_perc, self.normalized, self.std_inh, self.std,
                                              self.normalizedtranslation, self.datamedium, self.inhibitedmedium,
                                              self.normalized_percmedium, self.normalizedmedium, self.std_inhmedium,
                                              self.stdmedium, self.confinterval, self.normalizedtranslationmedium)
        self.get_combination(self.brpath, 3)
        # self.get_medium_control(self.file, 3)

        self.close_file()


class BRgeneticchemicalperturbagem(BRHdf5database):
    def __init__(self, count, file_name_br, experiment_type, branch, file_name,
                 header, elapse, date_info, date, data, std, medium):
        super().__init__(count, file_name_br, experiment_type, branch, file_name,
                         header, elapse, date_info, date, data, std, medium)

        self.celline = []
        self.seeding = []
        self.branch = []
        self.condition = []
        # self.control_name = []
        self.std_info = ["STD_BR"]
        self.control_data_single = None

        self.open_file()
        # compute mean and STD
        self.get_pathways()
        # normalize only the mean
        self.normalizeseveralcontrol()
        self.hdfnorm.get_fieldsseveralcontrolmain(self.brpath, 5, 1, self.fields, self.data, self.inhibited,
                                                  self.normalized_perc, self.normalized, self.std_inh, self.std,
                                                  self.normalizedtranslation, self.datamedium, self.inhibitedmedium,
                                                  self.normalized_percmedium, self.normalizedmedium, self.std_inhmedium,
                                                  self.stdmedium, self.confinterval, self.normalizedtranslationmedium)

        self.hdfnorm.get_genetcombmain(self.brpath, 2, self.mediumword, self.possiblecombination,
                                       self.possiblecombinationsize)
        self.get_medium_control(self.brpath, 3)

        self.compoundalonearranged = self.hdfnorm.get_branch(1, self.brpath, self.branch,
                                                             self.compoundalone)
        self.compoundalone = self.compoundalonearranged
        self.close_file()

class BRgeneticperturbagem(BRHdf5database):
    def __init__(self, count, file_name_br, experiment_type, branch, file_name,
                 header, elapse, date_info, date, data, std, medium):
        super().__init__(count, file_name_br, experiment_type, branch, file_name,
                         header, elapse, date_info, date, data, std, medium)
        self.celline = []
        self.seeding = []
        self.condition = []
        # self.control_name = []
        self.std_info = ["STD_BR"]
        # this is the only experimental condition that accepts more than 1 control
        self.open_file()
        # compute mean and STD
        self.get_pathways()

        self.normalizeonecontrol()

        self.hdfnorm.get_fieldsonecontrolmain(self.brpath, 5, 1, self.fields, self.control, self.data, self.inhibited,
                                              self.normalized_perc, self.normalized, self.std_inh, self.std,
                                              self.normalizedtranslation, self.datamedium, self.inhibitedmedium,
                                              self.normalized_percmedium, self.normalizedmedium, self.std_inhmedium,
                                              self.stdmedium, self.confinterval, self.normalizedtranslationmedium)
        # self.get_fieldsonecontrol(self.brpath, 5)
        self.get_compoundalone(self.brpath, 3)

        self.close_file()
