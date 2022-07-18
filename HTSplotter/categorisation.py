import numpy as np
import sys
import os


class Categorisation:
    def __init__(self, header_info):

        # from header object in main
        self.headers = header_info.newheader
        self.medium = header_info.medium
        self.control = header_info.controlname
        self.stdinfo = header_info.stdinfo

        self.keylistmaingroups = ["Cell", "Seeding", "Condition", "Compound", "Concentration", "Unit", "Position"]
        self.diccompoundgroupkey = {}  # contains compound tested
        self.dicmediumgroupkey = {}  # contains all medium and cells only
        self.diccontrolgroupkey = {}  # contains all controls

        self.cellline = []
        self.seeding = []
        self.condition = []
        self.compound = []
        self.compoundmedium = []
        self.branch = []
        self.combination = []
        self.concentration = []

        self.group = []
        # get number of cell lines, number of different seeding

        # Experiment categorization :
        # compound screen=> Concentration >3; Condition = "Condition"; Combination = 0
        # compound combination => Condition = "Condition"; Combination != 0
        # genetic perturbagen => Concentration = 1; Condition = "Condition" ; Combination = 0
        # genetic-chemical perturbagen => Condition != "Condition"; Combination !=0
        # "Compound screen"; "Genetic perturbagen"; "Compound combination";"genetic-chemical perturbagen"
        self.experimentype = []  #
        self.compoundscreen = []

        #

        self.get_groups()
        self.get_groupsnumber()
        self.get_branch()
        self.get_experimentype()

    def get_groups(self):
        notmedium = 0
        notcontrol = 0
        for i in self.headers:
            if i[0] not in self.cellline and len(i[0]) != 0:
                self.cellline.append(i[0])
            if i[1] not in self.seeding and len(i[1]) != 0:
                self.seeding.append(i[1])
            if i[2] not in self.condition and len(i[2]) != 0:
                self.condition.append(i[2])
            if i[3] not in self.compound and len(i[3]) != 0:
                if i[3] not in self.medium and i[3] not in self.control:
                    self.compound.append(i[3])
            if i[3] not in self.compoundmedium:
                if i[3] not in self.control and len(i[3]) != 0:
                    self.compoundmedium.append(i[3])
                if "_" in i[3] and i[3] not in self.combination and i[3] not in self.control:
                    self.combination.append(i[3])

            if i[:4] not in self.group and len(i[:4]) != 0:
                self.group.append(i[:4])

    def get_groupsnumber(self):
        # make dictionalryies for each conditions that belongs to the main groups : compound, medium, control
        # group contains all situations
        self.diccompoundgroupkey = {comp: None for comp in range(len(self.compound))}
        self.dicmediumgroupkey = {med: None for med in range(len(self.medium))}
        self.diccontrolgroupkey = {contr: None for contr in range(len(self.control))}

        for cond in self.headers:
            for k in range(len(self.medium)):
                if self.medium[k] == cond[3]:
                    meddict = self.get_dictionary(cond)
                    self.dicmediumgroupkey[k] = meddict
            for j in range(len(self.control)):
                if self.control[j] == cond[3]:
                    contdict = self.get_dictionary(cond)
                    self.diccontrolgroupkey[j] = contdict
            for i in range(len(self.compound)):
                if self.compound[i] == cond[3]:
                    compdict = self.get_dictionary(cond)
                    if self.diccompoundgroupkey[i] == None:
                        self.diccompoundgroupkey[i] = compdict
                    else:
                        if compdict["Cell"] == self.diccompoundgroupkey[i]["Cell"]:
                            if compdict["Seeding"] == self.diccompoundgroupkey[i]["Seeding"]:
                                if compdict["Condition"] == self.diccompoundgroupkey[i]["Condition"]:
                                    if compdict["Compound"] == self.diccompoundgroupkey[i]["Compound"]:
                                        self.diccompoundgroupkey[i]["Concentration"].append(
                                            compdict["Concentration"][0])
                                        self.diccompoundgroupkey[i]["Unit"].append(
                                            compdict["Unit"][0])
                                        self.diccompoundgroupkey[i]["Position"].append(compdict["Position"][0])

    def get_branch(self):
        list = []
        for i in self.cellline:
            for j in self.seeding:
                for h in self.condition:
                    ctemp = []
                    for u in self.compoundmedium:
                        for k in self.group:
                            if j in k[1] and i in k[0] and h in k[2] and u == k[-1] and u not in ctemp:
                                ctemp.append(u)
                        li = [i, j, h, ctemp]
                        if li not in self.branch and len(ctemp) != 0:
                            self.branch.append(li)

    def get_concentration(self):
        for i in self.diccompoundgroupkey:
            self.concentration.append(len(self.diccompoundgroupkey[i]["Concentration"]))

    def get_experimentype(self):
        # Experiment categorization :
        # compound screen=> Concentration >3; Condition = "Condition"; Combination = 0
        # compound combination => Condition = "Condition"; Combination != 0
        # genetic perturbagen => Concentration = 1; Condition = "Condition" ; Combination = 0
        # genetic-chemical perturbagen => Condition != "Condition"; Combination !=0
        compoundscreen = 0  # not
        geneticpertur = 0  # not
        self.get_concentration()

        if len(self.condition) == 1 and len(self.combination) == 0:
            for i in self.diccompoundgroupkey:
                if len(self.diccompoundgroupkey[i]["Concentration"]) >= 2:
                    compoundscreen = 1
                if len(self.diccompoundgroupkey[i]["Concentration"]) == 1:
                    geneticpertur = 1
                    # self.experimentype = "Genetic_perturbagen"
            self.compoundscreen.append(len(self.diccompoundgroupkey[i]["Concentration"]))
            if compoundscreen == 1:
                self.experimentype = "drug_screen"
                print(self.experimentype, len(self.compound), len(self.control), self.stdinfo[0])
            if geneticpertur == 1 and compoundscreen == 0:
                self.experimentype = "Genetic_perturbagen"
                print(self.experimentype, len(self.control), len(self.compound))
        if len(self.combination) > 0 and len(self.condition) == 1 and len(self.condition) == len(self.cellline):
            self.experimentype = "drug_combination"
            print(self.experimentype, len(self.combination), len(self.control), self.stdinfo[0])
        if len(self.condition) > 1 and len(self.combination) > 0:
            self.experimentype = "genetic-chemical_perturbagen"
            print(self.experimentype, len(self.control), len(self.compound), self.stdinfo[0])
        if len(self.experimentype) == 0:
            self.experimentype = 'something went wrong'

    @staticmethod
    def get_dictionary(templist):
        tempdict = {}
        if len(templist[0]) == 0:
            tempdict["Cell"] = 'Unidentified'
        if len(templist[0]) != 0:
            tempdict["Cell"] = templist[0]
        if len(templist[2]) == 0:
            tempdict["Seeding"] = 'Unidentified'
        if len(templist[2]) != 0:
            tempdict["Seeding"] = templist[1]
        if len(templist[2]) == 0:
            tempdict["Condition"] = 'Unidentified'
        if len(templist[2]) != 0:
            tempdict["Condition"] = templist[2]
        if len(templist[3]) == 0:
            tempdict["Compound"] = 'Unidentified'
            tempdict["Concentration"] = ['unidentified']
            tempdict["Unit"] = ['unidentified']
            tempdict["Position"] = ['unidentified']
        if len(templist[3]) != 0:
            tempdict["Compound"] = templist[3]
            b = templist[4].split(" ")
            tempdict["Concentration"] = [b[0]]
            tempdict["Unit"] = [b[1]]
            tempdict["Position"] = [templist[5]]

        return tempdict
