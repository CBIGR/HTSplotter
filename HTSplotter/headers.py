import numpy as np
import sys
import os

class Headers:
    def __init__(self, headers):
        self.headersinitial = headers

        # #Get main information
        self.stdinfo = None
        self.headersinfo = []  # header once it is none if the data has not hot the std on headers
        self.stdheader = []
        # headers for no std
        self.headersgroup = []  # headers with data position if repetetive
        self.intermediaryheader = []  # # intermediary header, structure cell line and repetetive positions

        self.newheader = []  # final header, structure accepted by HTSplotter
        self.outputheader = []  # header for the output file
        # information to be saves on a txt file
        self.celline = []
        self.seeding = []
        self.conditions = []
        self.compound = []
        self.units = []

        self.controlname = []
        self.situations = []
        self.branch = []
        self.mysinglecompounddict = {"SingleCompound": [], "Concentration": [], "Units": []}
        self.mycombodict = {"Compoundcomb": [], "Concentration": [], "Units": []}
        # get sub informations
        self.controlheader = []
        self.cellcompinf = []
        self.medium = []
        self.mediumcondition = []
        self.compoundcompinfo = []
        self.controlinfo = []  # get the position of the control
        self.compoundsingle = []
        self.concencompound = []
        self.compondcomb = []  # get combination names
        self.numbercombinations = []  # get the number of combinations
        self.combinationdetail = []  # get detailed information about the combination
        self.experiment_type = None

        # error information
        self.errorheader = []

        # experiment definition
        self.experimenttype = None
        # Run main methods
        # self.split_compound_cell()
        self.get_stdinfo()

    def get_stdinfo(self):
        stdinfo = []
        for i in self.headersinitial:
            a = i.split(" ")
            count = 0
            for j in a:
                if '(S' in j or '(s' in j:
                    self.stdheader.append(i)
                    std = " ".join(a[count:]).replace("(", "").replace(")", "")
                    if std not in stdinfo:
                        stdinfo.append(std)
                count += 1
        if len(stdinfo) == 0:
            self.stdinfo = stdinfo
        else:
            self.stdinfo = stdinfo  # gives the type of STD and if it exists

    # ##methods to be used by main
    def headerswithstd(self):
        # split in half
        self.headersinfo = self.headersinitial[:(len(self.headersinitial) - len(self.stdheader))]
        if len(self.headersinfo) != len(self.stdheader):
            print("the number of STD headers is different from the experimental condition! "
                  "Please, check your inputfile")
            self.errorheader.append("the number of STD headers is different from the experimental condition! "
                                    "Please, check your inputfile")
        # get split each condition by compound and cell line
        self.getsplit_compound_cell()
        # get a intermediary header structure--> [cell line information and then compounds, position from data collumns]
        self.get_cell_compound_order()
        # gather repetitive information and create a intermediary header
        self.verify_repetitive_conditions()
        self.replacesybol()
        # crucial step to have the final headers structure!
        self.get_maininfo()
        # self.get_header_lend()
        # try:
        #     self.get_branch()
        # except:
        #     pass
        # self.get_experimenttype()

    def headersnostd(self):
        self.headersinfo = self.headersinitial
        # get split each condition by compound and cell line
        self.getsplit_compound_cell()
        # get a header structure--> [cell line information and then compounds, position from data collumns]
        self.get_cell_compound_order()
        # gather repetitive information and create a intermediary header
        self.verify_repetitive_conditions()
        self.replacesybol()
        # crucial step to have the final headers structure!
        self.get_maininfo()
        # self.get_header_lend()
        # try:
        #     self.get_branch()
        # except:
        #     pass
        #
        # self.get_experimenttype()

    # ##methods to be used by the headers
    def getsplit_compound_cell(self):
        count = 0
        for i in range(len(self.headersinfo)):
            # verify if the data has double space, to adjust to 1 space
            if "  " in self.headersinfo[i]:
                self.headersinfo[i] = " ".join(self.headersinfo[i].split())
            # txt file, each condition is separated by comma
            if "," in self.headersinfo[i]:
                a = self.headersinfo[i].split(",")
                b = " ".join(a[:-1])
                self.compoundcompinfo.append(b)
                self.cellcompinf.append(a[-1])
                count += 1
            # csv file, each condition is separated by ;
            if ";" in self.headersinfo[i]:
                a = self.headersinfo[i].split(";")
                self.compoundcompinfo.append(a[0])
                self.cellcompinf.append(a[-1])
                count += 1
            # txt file, each condition is separated by space
            if "," not in self.headersinfo[i]:
                # important consideration (the compound information is from a[(-1-3):-1]
                a = self.headersinfo[i].split(" ")
                if "(" in a[1]:
                    new = a[5:-1]
                    new = " ".join(new)
                    cel = " ".join(a[:5])
                else:
                    new = a[4:]
                    new = " ".join(new)
                    cel = " ".join(a[:4])
                self.compoundcompinfo.append(new)
                self.cellcompinf.append(cel)
                count += 1
        if count != len(self.headersinfo):
            # save to txterror file
            print('Not all headers have cell line information separated from compound information')
            b = 'Not all headers have cell line information separated from compound information'
            self.errorheader.append(b)
            # sys.exit()  # stops the program

    def get_cell_compound_order(self):
        # get a header experimental condition with the following structure:
        # cellline, seeding, compound, concentration, units, compoundx, concentration, unit, at the final position
        # the position of the data column
        count = 0
        for i in range(len(self.compoundcompinfo)):
            cell = []
            b = self.compoundcompinfo[i].split(" ")
            a = self.cellcompinf[i].split(" ")
            b.append(count)
            cell.append(a[0])
            # in case of having a different information between ()
            if "(" in a[1]:
                joined = "".join(a[2:])
                cell.append(joined)
            else:
                joined = "".join(a[1:])
                cell.append(joined)
            self.headersgroup.append(cell + b)
            count += 1

    def replacesybol(self):
        # replace symbol "/" by "per"
        for i in range(len(self.intermediaryheader)):
            for j in range(len(self.intermediaryheader[i])):
                try:
                    self.intermediaryheader[i][j] = self.intermediaryheader[i][j].replace("/", "per")
                except:
                    pass

    def verify_repetitive_conditions(self):
        for i in self.headersgroup:
            if i[:-1] not in self.intermediaryheader:
                self.intermediaryheader.append(i[:-1])
        for i in self.intermediaryheader:
            i.append([])
        for k in self.headersgroup:
            for j in self.intermediaryheader:
                if k[:-1] == j[:-1]:
                    j[-1].append(k[-1])
            # self.intermediaryheader.append(i[:-1])

    def get_maininfo(self):
        # the intermediary header has the followin structure:
        # first position always the cell line name
        # second position always the seeding
        # from the third to the penultimate position we have the compounds
        # the last position is the collum number from the input file
        control = 0  # indicates if the contidion is a control
        for i in self.intermediaryheader:
            newheader = []
            outheader= []
            # identify if the condition is a control. This information has to be on the cell line
            if "Control" in i[0] or "control" in i[0]:
                # remove "control" word from the cell line
                a = i[0].split("_")
                i[0] = "_".join(a[:-1])
                control = 1  # indicates that the header is a control
                # get the position of the control
                self.controlinfo.append(i[-1])

            self.get_celllineinfo(i[0], i[1])

            comp, conc_unit, condition, control, units, concentration = self.get_compoundinfo(i[2:-1], control)
            # if comp not in self.compound and "Control" not in comp:
            #     print(">>>>",comp)
            #     self.compound.append(comp)
            if condition not in self.conditions:
                self.conditions.append(condition)
            newheader.append(i[0])
            newheader.append(i[1])
            newheader.append(condition)
            newheader.append(comp)
            newheader.append(conc_unit)
            newheader.append(i[-1])
            outheader.append(i[0])
            outheader.append(i[1])
            outheader.append(condition)
            outheader.append(comp)
            outheader.append(concentration)
            outheader.append(units)
            outheader.append(i[-1])
            self.newheader.append(newheader)
            if 'Control' in comp:
                self.controlheader.append(outheader)
            self.outputheader.append(outheader)

    def get_celllineinfo(self, i, j):
        if i not in self.celline:
            self.celline.append(i)
        if j not in self.seeding:
            self.seeding.append(j)

    def get_compoundinfo(self, f, infocontrol):
        self.get_medium(f)
        comp = []
        concentration = []
        units = []
        conc_uni = []
        condition = []
        if len(f) % 2 == 0:
            if len(f)/2 != 3:
                new = f[-1][-1]
                f[-1] = f[-1].replace(new, '')
                if new.isdigit():
                    self.errorheader.append(f)
                else:
                    f.append(new)
        if len(f) % 3 == 0:
            comp, condition, concentration, units, conc_uni, infocontrol = self.name_len_divi_3(f, infocontrol)
        if infocontrol == 1:
            comp = "Control_" + comp
            if comp not in self.controlname:
                self.controlname.append(comp)
            infocontrol = 0
        if len(concentration) == 0:
            self.errorheader.append(concentration)
        if len(concentration) > 0:
            if not concentration.isdigit():
                self.errorheader.append(concentration)
        if len(condition) == 0:
            condition = "Condition"

        return comp, conc_uni, condition, infocontrol, units, concentration

    def get_comp_group(self, f, comp, condition, concentration, units, conc_uni):
        # if len(f) = 3 --> means 1: compound, concentration, units
        # if len(f)/3 != 1--> means more than 1 compound
        if len(f) == 3:
            if comp not in self.compoundsingle and comp != self.medium and "Control" not in comp:
                # in this condition we have 1 compound that is not medium or control
                self.compoundsingle.append(comp)
            self.mysinglecompounddict["SingleCompound"].append(comp)
            self.mysinglecompounddict["Concentration"].append(concentration)
            self.mysinglecompounddict["Units"].append(units)
        if len(f)/3 != 1:
            comb = "_".join(f[0::3])
            if comb not in self.compondcomb:
                self.numbercombinations.append(len(f[0::3]))
                self.compondcomb.append(comp)
            self.mycombodict["Compoundcomb"].append(comp)
            self.mycombodict["Concentration"].append(concentration)
            self.mycombodict["Units"].append(units)

        if comp not in self.compoundsingle and comp != self.medium and 'Control' not in comp and len(f) == 3:
            self.compound.append(comp)
            self.compoundsingle.append(comp)

    @staticmethod
    def name_len_divi_2(f, comp, concentration, units):
        condition = []
        for j in range(int(len(f) / 2)):
            comp.append(f[j * 2])
            concentration.append((f[j * 2 + 1]))
            units.append(f[j * 2 + 2])
        if len(comp) > 1:
            comp = "_".join(comp)
            concentration = "_".join(concentration)
            units = "_".join(units)
            conc_uni = concentration + " " + units
        elif len(comp) == 1:
            if "_" in comp[0]:
                g = comp[0].split("_")
                comp = g[0]
                condition = g[-1]
            else:
                comp = comp[0]
            concentration = concentration[0]
            units = units[0]
            conc_uni = concentration + " " + units

        return comp, condition, concentration, units, conc_uni

    @staticmethod
    def name_len_divi_3(f, infocontrol):
        condition = []
        comp = f[0::3]
        concentration = f[1::3]
        units = f[2::3]
        if len(comp) > 1:
            comp = ["_".join(comp)]
            concentration = "_".join(concentration)
            units = "_".join(units)
            conc_uni = concentration + " " + units
        elif len(comp) == 1:
            if "_" in comp[0]:
                g = comp[0].split("_")
                if len(g) > 1 and infocontrol != 1:
                    condition = g[-1]
                if len(g) >= 3 and infocontrol == 1:
                    condition = g[-1]
            concentration = concentration[0]
            units = units[0]
            conc_uni = concentration + " " + units
        return comp[0], condition, concentration, units, conc_uni, infocontrol

    # @staticmethod
    def get_medium(self, name):
        medium = ["Only_Cells", "OnlyCells", "only_cells", "Onlycells", "onlycells",
                  "onlyCells", "XXOnlyCells", "CellsOnly", "cellsonly", 'OnlyCells_', 'CellsOnly_']

        if name[0] in medium:
            # name = "zzzOnlyCells"
            # name[0] = "OnlyCells"
            print(name[0])
            self.medium.append(name[0])  #name
            # name[0] = "zzzOnlyCells"
        if '_' in name[0]:
            a = name[0].split("_")
            if a[0] in medium:
                name[0] = "OnlyCells" + "_" + a[1]
                self.medium.append(name[0])
