import h5py
import os
import numpy as np
from scipy import stats
## common functions between save_hdf5file and save_hdf5brfiles

class Hdf5functions:

    def __init__(self):
        self.f = None
        self.celline = None
        self.seeding = None
        self.condition = None
        self.conditionpath = None
        self.controlpath = None
        self.controlpathdata = None
        self.mediumword = None
        self.medium = None
        self.data = None
        self.inhibited = None
        self.normalized_perc = None
        self.normalized = None
        self.std_inh = None
        self.std = None
        self.normalizedtranslation = None
        self.datamedium = None
        self.inhibitedmedium = None
        self.normalized_percmedium = None
        self.normalizedmedium = None
        self.std_inhmedium = None
        self.stdmedium = None
        self.normalizedtranslationmedium = None
        self.fields = None
        self.confinterval = None
        self.br = None

    # Normalizations
    def get_normalizeseveralcontrol(self, f, level, exp, k, celline,
                                    seeding, condition, controlpath, controlpathdata, medium, mediumword):
        self.f = f
        self.celline = celline
        self.seeding = seeding
        self.condition = condition
        self.controlpath = controlpath
        self.controlpathdata = controlpathdata
        self.mediumword =mediumword
        self.medium = medium

        level += 1
        if exp != 0:
            cellpath = self.f[exp][k]
        else:
            cellpath = self.f[k]
        self.celline.append(k)
        seeding = list(cellpath.keys())
        for j in seeding:
            self.seeding.append(j)
            level += 1
            seedingpath = cellpath[j]
            condition = list(seedingpath.keys())
            self.condition.append(condition)
            for g in condition:
                level += 1
                conditionpath = seedingpath[g]
                compound = list(conditionpath.keys())
                for c in compound:
                    compoundpath = conditionpath[c]
                    for i in range(len(self.controlpath)):
                        if compoundpath.name == self.controlpath[i]:
                            dosage = list(compoundpath.keys())
                            for d in dosage:
                                grup = compoundpath[d]
                                controlname = self.controlpath[i].split("/")
                                self.hdf5normalization(grup, controlname[-1], self.controlpathdata[i])
                        elif "Condition" not in compoundpath.name:
                            splitpath = self.controlpath[i].split("/")
                            pathmain = "/".join(splitpath[:-2])
                            pathcontrolcompound = splitpath[-1]
                            if pathmain in compoundpath.name and pathcontrolcompound in compoundpath.name:
                                dosage = list(compoundpath.keys())
                                for d in dosage:
                                    grup = compoundpath[d]
                                    controlname = self.controlpath[i].split("/")
                                    self.hdf5normalization(grup, controlname[-1], self.controlpathdata[i])
                            elif self.mediumword in c:
                                singledosage = list(compoundpath.keys())
                                grup = compoundpath[singledosage[0]]
                                new4 = self.controlpath[i].split("/")
                                new3 = "/".join(new4[:-2])
                                if new3 in conditionpath.name:
                                    self.hdf5normalization(grup, new4[-1], self.controlpathdata[i])
                        if c in self.medium:
                            singledosage = list(compoundpath.keys())
                            grup = compoundpath[singledosage[0]]
                            new = self.controlpath[i].split("/")
                            new2 = "/".join(new[:-1])
                            if new2 in conditionpath.name:
                                self.hdf5normalization(grup, new[-1], self.controlpathdata[i])
                        if "Control" in c:
                            singledosage = list(compoundpath.keys())
                            grup = compoundpath[singledosage[0]]
                            new = compoundpath.name.split("/")
                            new2 = new[-1].split("_")
                            if len(new2) > 3:
                                new3 = "_".join(new2[-2:])
                                new4 = "/".join(new[:-1])
                                complete = new4 + "/" + new3
                                if complete == self.controlpath[i]:
                                    self.hdf5normalization(grup, new[-1], self.controlpathdata[i])
                            else:
                                new3 = "_".join(new2[-1:])
                                new4 = "/".join(new[:-1])
                                complete = new4 + "/" + new3
                                if complete == self.controlpath[i]:
                                    self.hdf5normalization(grup, new[-1], self.controlpathdata[i])

    def get_normalizeonecontrol(self, f, exp, k, celline, seeding, condition, conditionpath, controlpath,
                                controlpathdata, medium, mediumword):

        self.f = f
        self.celline = celline
        self.seeding = seeding
        self.condition = condition
        self.conditionpath = conditionpath
        self.controlpath = controlpath
        self.controlpathdata = controlpathdata
        self.mediumword = mediumword
        self.medium = medium
        if exp != 0:
            cellpath = self.f[exp][k]
        else:
            cellpath = self.f[k]

        seeding = list(cellpath.keys())
        self.celline.append(k)
        for j in seeding:
            seedingpath = cellpath[j]
            condition = list(seedingpath.keys())
            self.seeding.append(j)
            for g in condition:
                conditionpath = seedingpath[g]
                compound = list(conditionpath.keys())
                self.condition.append(g)
                self.conditionpath.append(conditionpath.name)
                for c in compound:
                    compoundpath = conditionpath[c]
                    dosage = list(compoundpath.keys())
                    for i in range(len(self.controlpath)):
                        for d in dosage:
                            grup = compoundpath[d]
                            controlname = self.controlpath[i].split("/")
                            controlmainpath = "/".join(controlname[:-1])
                            if controlmainpath == conditionpath.name:
                                self.hdf5normalization(grup, controlname[-1], self.controlpathdata[i])

    @staticmethod
    def hdf5normalization(grp, control, controldata):

        if 'data' in list(grp):
            grp.create_dataset('normalizedperc' + "_" + control,
                               data=np.divide(grp['data'][:], controldata) * 100)
            grp.create_dataset('normalized' + "_" + control,
                               data=np.divide(grp['data'][:], controldata))
            grp.create_dataset('normalized_translation' + "_" + control,
                               data=np.divide(grp['data'][:], controldata) - 1)
            grp.create_dataset('inhibited' + "_" + control,
                               data=(1 - np.divide(grp['data'][:], controldata)))
            grp.create_dataset('stdinh' + "_" + control,
                               data=(np.divide(grp['std'][:], controldata)))
        else:
            grp.create_dataset('normalizedperc' + "_" + control,
                               data=np.divide(grp['mean'][:], controldata) * 100)
            grp.create_dataset('normalized' + "_" + control,
                               data=np.divide(grp['mean'][:], controldata))
            grp.create_dataset('normalized_translation' + "_" + control,
                               data=np.divide(grp['mean'][:], controldata) - 1)
            grp.create_dataset('inhibited' + "_" + control,
                               data=(1 - np.divide(grp['mean'][:], controldata)))
            grp.create_dataset('stdinh' + "_" + control,
                               data=(np.divide(grp['std_BR'][:], controldata)))

    # combinations common for drug combination and genetic-chemical perturbagem
    def get_combination(self, f, compoundalone, possiblecombination, possiblecombinationsize):
        self.compoundalone = compoundalone
        self.possiblecombination = possiblecombination
        self.possiblecombinationsize = possiblecombinationsize

        combination = []
        alone = []
        combinsize = []
        compound = list(f.keys())
        for i in compound:
            temp = []
            if "_" in i and "Control" not in i:
                each = i.split("_")
                temp.append(i)
                for j in each:
                    temp.append(j)
                if len(temp) > 0:
                    combination.append(temp)
                    combinsize.append([])
            elif "Control" not in i:
                alone.append(i)
        if len(alone) > 0:
            self.compoundalone.append(alone)
        if len(combination) > 0:
            self.possiblecombination.append(combination)
        for comb in range(len(combination)):
            for m in range(len(combination[comb])):
                for com in compound:
                    tes = f[com]
                    if com == combination[comb][m]:
                        conce = len(list(tes.keys()))
                        combinsize[comb].append(conce)
        if len(combinsize) > 0:
            self.possiblecombinationsize.append(combinsize)

    # add combination results to hdf5 file--> common for drugcombination and genetic-chemical perturbagem
    def add_comboinformation(self, groupname, biscore, predicted, path, possiblecombination, conditionpath,
                             possiblecombinationsize):

        for i in range(len(possiblecombination)):
            for j in range(len(possiblecombination[i])):
                grp = path.create_group(conditionpath[i] + "/" + "combination_"
                                        + possiblecombination[i][j][0])

                for u in range(possiblecombinationsize[i][j][0]):
                    grp2 = grp.create_group(groupname[i][j][u])
                    grp2 = grp[groupname[i][j][u]]
                    grp2.create_dataset('predicted_value', data=predicted[i][j][:, u])
                    grp2.create_dataset('BI_Score', data=biscore[i][j][:, u])

    # common to drug screens with more than 1 control and genetic-chemical perturbagem
    def get_fieldsseveralcontrolmain(self, f, level, br, fields, data, inhibited, normalized_perc, normalized,
                                     std_inh, std, normalizedtranslation, datamedium, inhibitedmedium, normalized_percmedium,
                                     normalizedmedium, std_inhmedium, stdmedium, confinterval,
                                     normalizedtranslationmedium):
        self.data = data
        self.inhibited = inhibited
        self.normalized_perc = normalized_perc
        self.normalized = normalized
        self.std_inh = std_inh
        self.std =std
        self.normalizedtranslation =normalizedtranslation
        self.datamedium = datamedium
        self.inhibitedmedium = inhibitedmedium
        self.normalized_percmedium = normalized_percmedium
        self.normalizedmedium =normalizedmedium
        self.std_inhmedium = std_inhmedium
        self.stdmedium = stdmedium
        self.normalizedtranslationmedium = normalizedtranslationmedium
        self.fields = fields
        self.confinterval = confinterval
        self.br = br

        self.get_fieldsseveralcontrol(f, level)

    def get_fieldsseveralcontrol(self, f, level):
        if level > 0:
            try:
                keys = list(f.keys())
                for key in keys:
                    grp1 = f[key]
                    self.get_fieldsseveralcontrol(grp1, level - 1)
            except AttributeError:
                pass
        else:
            makefields = f.name.split("/")
            if self.br == 1:
                listtemp = makefields[2:]

            else:
                listtemp = makefields[1:]
            if self.mediumword not in f.name:
                if len(listtemp) >= 5:
                    self.get_data(f, listtemp)
            else:
                if len(listtemp) >= 5:
                    processeddata = list(f.keys())
                    for k in processeddata:
                        if self.br == 1:
                            self.save_listbr(f, k, self.datamedium, self.inhibitedmedium,
                                             self.normalized_percmedium,
                                             self.normalizedmedium, self.std_inhmedium, self.stdmedium,
                                             self.confinterval, self.normalizedtranslationmedium)
                        else:
                            self.save_list(f, k, self.datamedium, self.inhibitedmedium,
                                           self.normalized_percmedium, self.normalizedmedium,
                                           self.std_inhmedium, self.stdmedium, self.normalizedtranslationmedium)


    def get_fieldsonecontrolmain(self, f, level, br, fields, control, data, inhibited, normalized_perc, normalized,
                                     std_inh, std, normalizedtranslation, datamedium, inhibitedmedium, normalized_percmedium,
                                     normalizedmedium, std_inhmedium, stdmedium, confinterval,
                                     normalizedtranslationmedium):
        self.data = data
        self.inhibited = inhibited
        self.normalized_perc = normalized_perc
        self.normalized = normalized
        self.std_inh = std_inh
        self.std =std
        self.normalizedtranslation =normalizedtranslation
        self.datamedium = datamedium
        self.inhibitedmedium = inhibitedmedium
        self.normalized_percmedium = normalized_percmedium
        self.normalizedmedium =normalizedmedium
        self.std_inhmedium = std_inhmedium
        self.stdmedium = stdmedium
        self.normalizedtranslationmedium = normalizedtranslationmedium
        self.fields = fields
        self.control = control
        self.confinterval = confinterval
        self.br = br

        self.get_fieldsonecontrol(f, level)

    def get_fieldsonecontrol(self, f, level):
        if level > 0:
            try:
                keys = list(f.keys())
                for key in keys:
                    grp1 = f[key]
                    self.get_fieldsonecontrol(grp1, level-1)
            except AttributeError:
                pass
        else:
            makefields = f.name.split("/")
            if self.br == 1:
                listtemp = makefields[2:]
            else:
                listtemp = makefields[1:]
            if len(listtemp) >= 5:
                self.get_data(f, listtemp)
                if "Control" in f.name:
                    self.control.append(listtemp)

    def get_data(self, f, listtemp):
        listtemp.append(len(self.fields))
        self.fields.append(listtemp)
        processeddata = list(f.keys())
        for k in processeddata:
            if self.br == 1:
                self.save_listbr(f, k, self.data, self.inhibited, self.normalized_perc, self.normalized,
                                 self.std_inh, self.std, self.confinterval, self.normalizedtranslation)
            else:
                self.save_list(f, k, self.data, self.inhibited, self.normalized_perc, self.normalized,
                               self.std_inh, self.std, self.normalizedtranslation)

    @staticmethod
    def save_list(f, k, datamedium, inhibitedmedium, normalized_percmedium, normalizedmedium,
                  std_inhmedium, stdmedium, normalizedtranslation):
        if 'data' == k:
            datamedium.append(f[k][:])
        elif 'inhibited' in k:
            inhibitedmedium.append((f[k][:]))
        elif 'normalizedperc' in k:
            normalized_percmedium.append((f[k][:]))
        elif 'normalized_translation' in k:
            normalizedtranslation.append((f[k][:]))
        elif 'normalized' in k:
            normalizedmedium.append((f[k][:]))
        elif 'stdinh' in k:
            std_inhmedium.append((f[k][:]))
        elif 'std' in k:
            stdmedium.append(f[k][:])

    @staticmethod
    def save_listbr(f, k, datamedium, inhibitedmedium, normalized_percmedium, normalizedmedium,
                  std_inhmedium, stdmedium, confinterval, normalizedtranslation):
        if 'mean' == k:
            datamedium.append(f[k][:])
        elif 'inhibited' in k:
            inhibitedmedium.append((f[k][:]))
        elif 'normalizedperc' in k:
            normalized_percmedium.append((f[k][:]))
        elif 'normalized_translation' in k:
            normalizedtranslation.append((f[k][:]))
        elif 'normalized' in k:
            normalizedmedium.append((f[k][:]))
        elif 'stdinh' in k:
            std_inhmedium.append((f[k][:]))
        elif 'std_BR' in k:
            stdmedium.append(f[k][:])
        elif '95%CI_BR' in k:
            confinterval.append(f[k][:])

    def medium_control(self, f, br, mediumword, control, compoundalone, fieldsmedium, fieldsmediuminhibited):
        self.br = br
        self.mediumword =mediumword
        self.control = control
        self.compoundalone = compoundalone
        self.fieldsmedium = fieldsmedium
        self.fieldsmediuminhibited = fieldsmediuminhibited

        compound = list(f.keys())
        alone = []
        mediumne = []
        inhmedium = []
        control = []
        correcpositi = 0
        correcpositiinhib = 0
        if len(self.fieldsmedium) > 0:
            posit = self.fieldsmedium[-1][-1][-1]
            correcpositi = 1
        if len(self.fieldsmediuminhibited) > 0:
            positinh = self.fieldsmediuminhibited[-1][-1][-1]
            correcpositiinhib = 1
        for i in compound:
            if "Control" in i and "Condition" in f.name:
                cont = f[i]
                k = list(cont.keys())
                cont2 = cont[k[0]]
                if self.br == 1:
                    control.append(cont2.name.split("/")[2:])
                else:
                    control.append(cont2.name.split("/")[1:])
            elif self.mediumword in i:
                medium = f[i]
                k = list(medium.keys())
                a = medium[k[0]]
                cond = list(a.keys())
                makefields = a.name.split("/")
                if self.br == 1:
                    listtemp = makefields[2:]
                else:
                    listtemp = makefields[1:]
                temmed = []
                teminh = []
                for c in cond:
                    lista2 = listtemp.copy()
                    lista3 = listtemp.copy()
                    if "_" in c and 'inhibited' in c:
                        lista3[3] += "_" + c.split("_")[-1]
                        if correcpositiinhib == 1:
                            posit += 1
                            lista3.append(positinh)
                        else:
                            lista3.append(len(inhmedium))
                        teminh.append(lista3[3])
                        inhmedium.append(lista3)
                    if lista2[3] not in temmed:
                        temmed.append(lista2[3])
                        if correcpositi == 1:
                            posit += 1
                            lista2.append(posit)
                        else:
                            lista2.append(len(mediumne))
                        mediumne.append(lista2)
            elif i not in alone:
                alone.append(i)
        if len(control) != 0:
            self.control.append(control)
        self.compoundalone.append(alone)
        if len(mediumne) > 0:
            self.fieldsmedium.append(mediumne)
        if len(inhmedium) > 0:
            self.fieldsmediuminhibited.append(inhmedium)

    # function only for genetic-chemical perturbagem
    def get_genetcombmain(self, f, level, mediumword, possiblecombination, possiblecombinationsize):
        self.f = f
        self.mediumword = mediumword
        self.possiblecombination = possiblecombination
        self.possiblecombinationsize = possiblecombinationsize
        level = level
        self.get_genetcomb(self.f, level)

    def get_genetcomb(self, f, level):
        if level > 0:
            try:
                keys = list(f.keys())
                for key in keys:
                    grp1 = f[key]
                    self.get_genetcomb(grp1, level - 1)
            except AttributeError:
                pass
        else:
            condition = list(f.keys())
            temp = []
            tempsize = []
            for i in condition:
                if i != 'Condition':
                    # this is the genetic perturbation in combination with the compound
                    combination = []
                    combinsize = []
                    comp = f[i]
                    compound = list(comp.keys())
                    for k in compound:
                        temp = []
                        temsize = []
                        if "_" in k and self.mediumword not in k:
                            each = k.split("_")
                            temp.append(k)
                            grp = comp[k]
                            temsize.append(len(list(grp.keys())))
                            temsize.append(len(list(grp.keys())))
                            temsize.append(1)
                            for j in each:
                                temp.append(j)
                            if len(temp) <= 3:
                                combination.append(temp)
                                combinsize.append(temsize)
                    self.possiblecombination.append(combination)
                    self.possiblecombinationsize.append(combinsize)

    def get_branch(self, br, file, branch, compoundalone):
        self.file = file
        self.branch = branch
        self.compoundalone = compoundalone
        self.compoundalonearranged = []
        self.br = br
        # cellpath, seedingpath, conditionpath, compound = self.hdf5walk(self.file, cell)
        # skip always 0Elapse; 000Date and 00Date_info
        if self.br == 1:
            lista1 = list(self.file.keys())
        else:
            lista1 = list(self.file.keys())[4:]
        for k in lista1:
            bran1 = [k]
            cellpath = self.file[k]
            for j in list(cellpath.keys()):
                seedingpath = cellpath[j]
                bran2 = bran1 + [j]
                for g in list(seedingpath.keys()):
                    bran3 = bran2 + [g]
                    com = [c for c in list(seedingpath[g].keys()) if "Control" not in c]
                    bran3.append(com)
                    self.branch.append(bran3)
        header_count = len(self.branch[::2])
        headers = [[] for i in range(0, header_count)]  # list_possible_combination
        newlist = headers
        co = 0
        count = 2
        for n in range(len(newlist)):
            newlist[n] = compoundalone[co:count]
            co = count
            count += 2
        self.compoundalonearranged = newlist
        return self.compoundalonearranged