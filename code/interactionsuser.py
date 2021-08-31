import os
import time

class Outputfile:
    def __init__(self, file, biologicalreplicate, information, error, concentration, stdinfo, diccompoundgroup,
                 diccontrolgroup, dicmediumgroup, experimentype, control, compound, combination, medium,
                 elapsed, condition, outputheaders, errorheader):

        self.filename = information
        self.biologicalreplicate = biologicalreplicate
        self.errorfile = error
        self.file = file
        self.stdinfo = stdinfo
        self.outputheader = outputheaders
        self.errorheader = errorheader

        self.diccompoundgroup = diccompoundgroup
        self.diccontrolgroup = diccontrolgroup
        self.dicmediumgroup = dicmediumgroup
        self.experimentype = experimentype
        self.control = control
        self.condition = condition
        self.compound = compound
        self.combination = combination
        self.concentration = concentration
        self.medium = medium
        self.elapsed = elapsed
        self.controlnameerror = None
        self.error = []
        self.readout = None # indicates the readout for the axis
        # self.errorheader = []
        self.keylistmaingroups = ["Cell", "Seeding", "Condition", "Compound", "Concentration", "Unit", "Position"]

        self.get_check()
        self.setfile()
        self.set_text()
        self.closefile()

    def get_check(self):
        # self.error.append("check experiment conditions information")
        for i in range(len(self.outputheader)):
            for j in range(len(self.outputheader[i])):
                if len(self.outputheader[i][j]) == 0:
                    self.error = 1
                    self.outputheader[i][j] = "unidentified"
                if "_" == self.outputheader[i][j]:
                    self.error = 1
                    self.outputheader[i][j] = "unidentified"
            if self.error == 1:
                self.errorheader.append(self.outputheader[i])
        if len(self.control) == 0:
            self.error = 1
            self.errorheader.append("Control unidentified")
            self.controlnameerror = "Control unidentified"
        if self.error == 1:
            self.error = "check error file, please"
        # if len(self.errorheader) != 0:
        #     for k in self.errorheader:
        #         self.errorheader.append(k)

    def setfile(self):
        self.outputfile = open(self.filename, 'w')

    def closefile(self):
        self.outputfile.close()

    def set_text(self):
        # Experiment categorization :
        # compound screen=> Concentration >3; Condition = "Condition"; Combination = 0
        # compound combination => Condition = "Condition"; Combination != 0
        # genetic perturbagen => Concentration = 1; Condition = "Condition" ; Combination = 0
        # genetic-chemical perturbagen => Condition != "Condition"; Combination !=0
        # "Compound screen"; "Genetic perturbagen"; "Compound combination";"genetic-chemical perturbagen"
        self.outputfile.write('File name: ' + self.file + '\n')
        self.outputfile.write('\n')
        self.outputfile.write('HTSplotter Categorizes the experiments as follows: \n')
        self.outputfile.write('\t drug => Concentration >= 2; Condition = "Condition"; Combination = 0 \n')
        self.outputfile.write('\t drug combination => Condition = "Condition"; Combination != 0 \n')
        self.outputfile.write('\t genetic perturbagen => Concentration = 1; Condition = "Condition" ; '
                              'Combination = 0 \n')
        self.outputfile.write('\t genetic-chemical perturbagen => Condition != "Condition"; Combination !=0 \n')
        self.outputfile.write('\n')
        self.outputfile.write('Based on the input file the categorization is : ' + self.experimentype + '\n')
        self.outputfile.write('\n\t number of concentration for each compound: ')
        for conc in self.concentration:
            self.outputfile.write(str(conc) + ';')
        self.outputfile.write('\n')
        self.outputfile.write('\n\t combination: ')
        for j in self.combination:
            self.outputfile.write(j + ' ')
        self.outputfile.write('\n')
        self.outputfile.write('\n\t control: ')
        for k in self.control:
            self.outputfile.write(k + ' ')
        self.outputfile.write('\n')
        self.outputfile.write('\n\t condition: ')
        for i in self.condition:
            self.outputfile.write(i + ' ')
        self.outputfile.write('\n')
        if self.biologicalreplicate != 0:
            self.outputfile.write('please check carefully if the time points interval '
                                  'is the same across your biological replicates!' + '\n')
        self.outputfile.write('\n')
        self.outputfile.write('\t the number of time points is: ' + str(len(self.elapsed)) + '\t')
        for ela in self.elapsed:
            self.outputfile.write(str(ela) + '\t')

        self.outputfile.write('\n')
        self.outputfile.write('\n')
        self.outputfile.write('\n')
        self.outputfile.write('The following information was obtained from the file header \n')
        self.outputfile.write('\n')
        self.outputfile.write('Information about standard deviation: ' + self.stdinfo[0] + '\n')
        self.outputfile.write('\n')
        self.outputfile.write('Information about compounds: \n')
        self.recursive_dic(self.outputfile, self.diccompoundgroup)
        self.outputfile.write('\n')
        self.outputfile.write('Information about control: \n')
        self.recursive_dic(self.outputfile, self.diccontrolgroup)
        self.outputfile.write('\n')
        if len(self.dicmediumgroup) > 0:
            self.outputfile.write('Information about medium: \n')
            self.recursive_dic(self.outputfile, self.dicmediumgroup)
            self.outputfile.write('\n')

    def seterrorfile(self, info):
        error = open(self.errorfile, 'w')
        if info == 0:
            error.write("Error was not identified" + '\n' + '\n')
            error.write("cellline, seeding, condition, compound, concentration, units, position from the input file" +
                        '\n' + '\n')
            # error.write("'unidentified', means the information is missing" + '\n' + '\n')
            # for j in self.outputheader:
            #     j[-1] = str(j[-1])
            #     for i in range(len(j)):
            #         try:
            #             if len(j[i]) == 0:
            #                 j[i] = "unidentified"
            #         except TypeError:
            #             j[i] = str(j[i])
            #
            #     error.write(", ".join(j) + '\n')
            # error.write('\n')
        if info == 1:
            error.write("You have an error from your header, please check bellow " + '\n' + '\n')
            error.write("Information order: cellline, seeding, condition, compound, "
                        "concentration, units, position from the input file" +
                        '\n' + '\n')
            error.write("'unidentified', means that the information is missing." + '\n' +
                        "\t" + " Between square brackets is the column position from your input file" + '\n' + '\n')
            if self.controlnameerror != None:
                error.write(self.controlnameerror + ' : '+
                            'Please indicate the control adding "_Control" to the cell line name')
                error.write('\n')

            error.write('\n')
            error.write('Input file headers: ' + '\n')
            for j in self.outputheader:
                j[-1] = str(j[-1])
                for i in range(len(j)):
                    try:
                        if len(j[i]) == 0:
                            j[i] = "unidentified"
                    except TypeError:
                        j[i] = str(j[i])

                error.write(", ".join(j) + '\n')
            error.write('\n')
        error.close()

    @staticmethod
    def recursive_dic(file, name):
        for i in name:
            file.write('\n')
            try:
                file.write("\t Cell line: " + name[i]["Cell"] + "\n")
                file.write("\t Seeding: " + name[i]["Seeding"] + "\n")
                file.write("\t Condition: " + name[i]["Condition"] + "\n")
                file.write("\t Compound: " + name[i]["Compound"] + "\n")
                file.write("\t Number of Concentration tested: " + str(len(name[i]["Concentration"])) + '\n')
                file.write("\t Concentration range: ")
                for j in name[i]["Concentration"]:
                    file.write("\t" + j + ' ')
                file.write("\n")
                file.write("\t Number of Units tested: " + str(len(name[i]["Unit"])) + '\n')
            except TypeError:
                print(name[i])


class Inputfile:
    def __init__(self, file, experiment_type):
        self.fileinput = file
        self.inputfile = None
        self.information = []
        self.userinformation = None
        self.experiment_type = experiment_type

        self.get_fileinformation()
        self.close_file_remove()
        self.get_experimenttype()
        self.removeinputfile()

    def get_fileinformation(self):
        self.inputfile = open(self.fileinput, 'r')
        a = self.inputfile.readlines()
        count = 0
        for linha in a:
            b = linha.split('\n')
            self.information.append(b[0])
            count += 1

    def close_file_remove(self):
        self.inputfile.close()

    def get_experimenttype(self):
        readout = self.information[-2].split(' ')
        print(readout)
        if len(readout[-1]) == 0:
            self.readout = 'confluency'
        else:
            self.readout = readout[-1]
        a = self.information[-1].split(' ')
        self.userinformation = a[-1]
        self.experiment_type = self.userinformation

    def removeinputfile(self):

        os.remove(self.fileinput)
