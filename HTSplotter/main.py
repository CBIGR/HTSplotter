import os
import time
import psutil
import datetime
import sys
import numpy as np
# librarys made for this script
from HTSplotter.readfiles import Readfile
from HTSplotter.filenames import Filenames
from HTSplotter.headers import Headers
from HTSplotter.categorisation import Categorisation
from HTSplotter.interactionsuser import Outputfile, Inputfile
from HTSplotter.save_hdf5file import Compoundscreenonecontrol, Compoundscreen, \
    Compoundcombination, Geneticperturbagem, Geneticchemicalperturbagem
from HTSplotter.save_hdf5brfiles import BRHdf5database, Individualcombinationstructure, BRcompoundscreenseveralcontrol, \
    BRcombinationstructure, BRcompoundscreenonecontrol, BRgeneticperturbagem, BRgeneticchemicalperturbagem, \
    Individualgeneticperturbagen
# from plotting import Overtime
from HTSplotter.combination import ExperimentCombination
from HTSplotter.geneticperturbagem import ExperimentGeneticPerturbagem
from HTSplotter.compoundscreen import SingleCompound, SingleCompoundonecontrol
from HTSplotter.geneticchemicalperturbagem import GeneticChemicalPerturbagem

# from txtsavedata import Savetxt

class Analyser:

    def __init__(self):

        self.main_folder = None
        self.input_path = None
        self.information_extracted = None
        self.results_path = None

        self.biological_replicate = 0  # 0  # 0 if not, 1 if it is
        self.userinput = 0  # 0  # 1 if yes
        self.information_readout = None  # "confluency"  # default = confluency = 0; add effect name = 1
        self.readout_units = None  # "(%)"
        self.expected_effect = 0  # 0  # 0 = "inhibition"; 1 = enhanced
        self.file_name_br = '0'  # "BiologicalReplicateteste"
        self.synergy_method = 0  # 0-->Bliss; 1-->HSA; 2--> ZIP

        self.files_list = None

    def _set_main_folder_(self, value):
        if isinstance(value, str):
            self.main_folder = value
        else:
            raise AttributeError

    def _set_input_path_(self, value):
        if isinstance(value, str):
            self.input_path = value
        else:
            raise AttributeError

    def _set_information_extracted_(self, value):
        if isinstance(value, str):
            self.information_extracted = value
        else:
            raise AttributeError

    def _set_results_path_(self, value):
        if isinstance(value, str):
            self.results_path = value
        else:
            raise AttributeError

    def _set_biological_replicate_(self, value):
        if isinstance(value, int):
            self.biological_replicate = value
        else:
            raise AttributeError

    def _set_user_input_(self, value):
        if isinstance(value, int):
            self.user_input = value
        else:
            raise AttributeError

    def _set_information_readout_(self, value):
        if isinstance(value, str):
            self.information_readout = value
        else:
            raise AttributeError

    def _set_readout_units_(self, value):
        if isinstance(value, str):
            self.readout_units = value
        else:
            raise AttributeError

    def _set_expected_effect_(self, value):
        if isinstance(value, int):
            self.expected_effect = value
        else:
            raise AttributeError

    def _set_file_name_br_(self, value):
        if isinstance(value, str):
            self.file_name_br = value
        else:
            raise AttributeError

    def _set_files_list_(self, value):
        if isinstance(value, list):
            self.files_list = value
        else:
            raise AttributeError

    def _set_files_list_(self, value):
        if isinstance(value, list):
            self.synergy_method = value
        else:
            raise AttributeError

    def execute_from_file(self, input_file):

        ifile = open(input_file, 'r')
        count_lines = 0
        for _ in ifile:
            count_lines += 1
        ifile.seek(0, 0)

        count_global = 0
        while count_global < count_lines:

            self.files_list = []
            counter = 0
            for line in ifile:
                if counter == 0:
                    self.main_folder = line.split()[0]
                elif counter == 1:
                    aux = line.split()
                    if aux[0] == 'main_folder':
                        if aux[1] == '+':
                            self.input_path = self.main_folder + aux[2]
                        else:
                            self.input_path = self.main_folder
                    else:
                        self.input_path = aux[0]
                elif counter == 2:
                    aux = line.split()
                    if aux[0] == 'main_folder':
                        if aux[1] == '+':
                            self.information_extracted = self.main_folder + aux[2]
                        else:
                            self.information_extracted = self.main_folder
                    else:
                        self.information_extracted = aux[0]
                elif counter == 3:
                    aux = line.split()
                    if aux[0] == 'main_folder':
                        if aux[1] == '+':
                            self.results_path = self.main_folder + aux[2]
                        else:
                            self.results_path = self.main_folder
                    else:
                        self.results_path = aux[0]
                elif counter == 4:
                    self.biological_replicate = int(line.split()[0])
                elif counter == 5:
                    self.userinput = int(line.split()[0])
                elif counter == 6:
                    self.information_readout = line.split()[0]
                elif counter == 7:
                    self.readout_units = line.split()[0]
                elif counter == 8:
                    self.expected_effect = int(line.split()[0])
                elif counter == 9:
                    self.file_name_br = line.split()[0]
                else:
                    aux = line.split()
                    if aux[0] != 'end':
                        self.files_list.append(line.split()[0])
                    else:
                        break

                counter += 1

            ifile.readline()

            self.execute()
            count_global += counter + 2

    def execute(self):
    # if __name__ == '__main__':
        # start = time.time()
        # Number of logical CPUs in the system
        # p = psutil.Process()
        self.verify_inputs()

        start = time.time()

        p = psutil.Process()


        if len(self.information_readout) == 0:
            information_readout = "no information"
        if len(self.readout_units) == 0:
            readout_units = 'no information'

        if self.biological_replicate == 1:
            print("Processing Bioloical replicate")
            count = 0
            file_br_names = Filenames(self.file_name_br, self.input_path, self.results_path, self.information_extracted)
            # first go over each file before processing it
            # Check confirmation to the user for all files
            for eachfile in self.files_list:
                file_names = Filenames(eachfile, self.input_path, self.results_path, self.information_extracted)
                # read input file, one by one, extract:
                # date ,elapsed , header , data ,data_std ,std ,std_header , date_info
                file_info = Readfile(file_names.fileitxtnputpath)
                header_info = Headers(file_info.header)

                # Check if the header has STD or not
                if len(header_info.stdinfo) != 0:
                    header_info.headerswithstd()
                    file_info.get_split_data_std(len(header_info.headersinitial) - len(header_info.stdheader))
                else:
                    header_info.headersnostd()
                    header_info.stdinfo.append("STD computed by HTSplotter")
                    file_info.get_data_nostd()

                catego = Categorisation(header_info)

                out = Outputfile(eachfile, self.biological_replicate, file_names, catego, header_info, file_info.elapsed)

                # set error file
                if len(out.error) == 0:
                    print("HTSplotter did not identify any error from your input file")
                    # out.seterrorfile(0)
                if len(out.error) != 0:
                    print(out.error)
                    out.seterrorfile(1)
                    sys.exit()
                # check user confirmation
                if self.userinput == 1:
                    time_counter = 0
                    while not os.path.exists(self.input_path + "inputfile.txt"):
                        print("still running, BR")
                        time.sleep(1)
                        time_counter += 1
                    if os.path.exists(self.input_path + "inputfile.txt"):
                        userconfi = Inputfile(self.input_path + "inputfile.txt", catego.experimentype)
                        print(catego.experimentype, userconfi.experiment_type)
                        information_readout = userconfi.readout
                        print("readout!!!===>", userconfi.readout)
                        if catego.experimentype != userconfi.experiment_type:
                            out.seterrorfile(0)
                            print("=====>check error file", userconfi.experiment_type, catego.experimentype)
                            break
                        print("=====>ok", userconfi.experiment_type, catego.experimentype)
                        #out.seterrorfile(0)
            # Once the user confirmed that each file is okay we can continue
            for i in self.files_list:
                file_names = Filenames(i, self.input_path, self.results_path, self.information_extracted)
                # read input file, one by one, extract:
                # date ,elapsed , header , data ,data_std ,std ,std_header , date_info
                file_info = Readfile(file_names.fileitxtnputpath)
                header_info = Headers(file_info.header)

                # Check if the header has STD or not
                if len(header_info.stdinfo) != 0:
                    header_info.headerswithstd()
                    file_info.get_split_data_std(len(header_info.headersinitial) - len(header_info.stdheader))
                else:
                    header_info.headersnostd()
                    header_info.stdinfo.append("STD computed by HTSplotter")
                    file_info.get_data_nostd()

                catego = Categorisation(header_info)

                if catego.experimentype == "Genetic_perturbagen":
                    indcombdata = Individualgeneticperturbagen(count, file_br_names.filehdf5resultspath, i,
                                                               header_info, catego, file_info)
                else:
                    indcombdata = Individualcombinationstructure(count, file_br_names.filehdf5resultspath, i,
                                                                 header_info, catego, file_info)
                    print("check--> ", count)

                count += 1

            # # check if the input files have the same compounds and conditions
            if catego.experimentype == "drug_combination":
                print("Compound_combinationBR", catego.experimentype)
                brcombdata = BRcombinationstructure(count, file_br_names.filehdf5resultspath, self.file_name_br,
                                                    header_info, catego, file_info)
                print(len(brcombdata.fields), len(brcombdata.data), len(brcombdata.std), len(brcombdata.std_inh),
                      len(brcombdata.normalized_perc), len(brcombdata.normalized))
                for g in range(len(brcombdata.possiblecombination)):
                    print(brcombdata.possiblecombinationsize[g], brcombdata.possiblecombination[g])

                comb = ExperimentCombination(self.synergy_method, brcombdata, file_info, catego.branch,
                                             information_readout, readout_units, file_br_names,
                                             self.biological_replicate, self.file_name_br)
                # ###
                # add predicted and bliss score to the HDF5 file
                brcombdata.add_predictedbiscore(comb.comb_name_per_group, comb.effect_per_group,
                                                comb.synergy_score_per_group, self.synergy_method)
                # over time data:

                if len(brcombdata.elapse) > 1:
                    comb.inhibition(brcombdata.data, brcombdata.std)
                    if self.expected_effect == 0:
                        comb.inhibition(brcombdata.inhibited, brcombdata.std_inh, 2, 1)
                    elif self.expected_effect == 1:
                        comb.inhibition(brcombdata.normalizedtranslation, brcombdata.std_inh, 1, 1)
                    comb.doseresponse()
                else:
                    print("1 time point for combination biological replicates")
                    if self.expected_effect == 0:
                        comb.endpointinhibition(brcombdata.inhibited, brcombdata.std_inh, 2, 1)
                    if self.expected_effect == 1:
                        comb.endpointinhibition(brcombdata.normalizedtranslation, brcombdata.std_inh, 1, 1)
                    comb.doseresponse(1)
                comb.close_pdf()

            if catego.experimentype == "drug_screen":
                print("Compound_screen", catego.experimentype)
                if len(catego.compound) == len(catego.control):
                    # 1 control for each compound
                    print("Compound_screen several control BR", len(catego.compound), len(catego.control), catego.medium)
                    brcomscreendata = BRcompoundscreenseveralcontrol(count, file_br_names.filehdf5resultspath,
                                                                     self.file_name_br, header_info, catego, file_info)

                    print(len(brcomscreendata.fields), len(brcomscreendata.data), len(brcomscreendata.std),
                          len(brcomscreendata.std_inh), len(brcomscreendata.normalized_perc),
                          len(brcomscreendata.normalized))
                    for g in range(len(brcomscreendata.compoundalone)):
                        print(brcomscreendata.compoundalone[g])
                    print(brcomscreendata.medium)
                    singlecompond = SingleCompound(header_info.branch, file_info,
                                                   information_readout, readout_units,
                                                   self.biological_replicate, brcomscreendata,
                                                   file_br_names, self.file_name_br)

                    if len(brcomscreendata.elapse) > 1:
                        singlecompond.inhibitionovertime(brcomscreendata.data, brcomscreendata.std)
                        if self.expected_effect == 0:
                            singlecompond.inhibitionovertime(brcomscreendata.inhibited, brcomscreendata.std_inh, 2, 1)
                        elif self.expected_effect == 1:
                            singlecompond.inhibitionovertime(brcomscreendata.normalizedtranslation, brcomscreendata.std_inh,
                                                             1, 1)
                        singlecompond.doseresponse()
                    else:
                        singlecompond.doseresponse(1)

                    singlecompond.close_pdf()
                else:
                    print("drug_screen 1 control BR", len(catego.compound), len(catego.control), catego.medium)
                    brcomscreenonecontrol = BRcompoundscreenonecontrol(count, file_br_names.filehdf5resultspath,
                                                                       self.file_name_br, header_info, catego,
                                                                       file_info)
                    print("one control")
                    print(len(brcomscreenonecontrol.fields), len(brcomscreenonecontrol.data),
                          len(brcomscreenonecontrol.std), len(brcomscreenonecontrol.std_inh),
                          len(brcomscreenonecontrol.normalized_perc), len(brcomscreenonecontrol.normalized))
                    for g in range(len(brcomscreenonecontrol.compoundalone)):
                        print(brcomscreenonecontrol.compoundalone[g])
                    print(brcomscreenonecontrol.medium)

                    singlecomponecontrol = SingleCompoundonecontrol(header_info.branch, file_info,
                                                                    information_readout, readout_units,
                                                                    self.biological_replicate, brcomscreenonecontrol,
                                                                    file_br_names, self.file_name_br)

                    if len(brcomscreenonecontrol.elapse) > 1:
                        singlecomponecontrol.inhibitiononecontrolovertime(brcomscreenonecontrol.data,
                                                                          brcomscreenonecontrol.std)
                        if self.expected_effect == 0:
                            singlecomponecontrol.inhibitiononecontrolovertime(brcomscreenonecontrol.inhibited,
                                                                              brcomscreenonecontrol.std_inh, 2, 1)
                        elif self.expected_effect == 1:
                            singlecomponecontrol.inhibitiononecontrolovertime(brcomscreenonecontrol.normalizedtranslation,
                                                                              brcomscreenonecontrol.std_inh, 1, 1)
                        singlecomponecontrol.doseresponseonecontrol()
                    else:
                        singlecomponecontrol.doseresponseonecontrol(1)

                    singlecomponecontrol.close_pdf()

            if catego.experimentype == "Genetic_perturbagen":
                print("Genetic_perturbagen")
                # save information on HDF5 file

                BRhdf5genetic = BRgeneticperturbagem(count, file_br_names.filehdf5resultspath,
                                                     self.file_name_br, header_info, catego, file_info)
                print("did it?>>>")
                print(len(BRhdf5genetic.fields), len(BRhdf5genetic.data), len(BRhdf5genetic.std),
                      len(BRhdf5genetic.std_inh), len(BRhdf5genetic.normalized_perc), len(BRhdf5genetic.normalized))
                print(BRhdf5genetic.compoundalone)

                geneticpert = ExperimentGeneticPerturbagem(header_info, BRhdf5genetic, file_info, catego, file_br_names,
                                                           information_readout, readout_units, self.file_name_br,
                                                           self.biological_replicate)

                if len(BRhdf5genetic.elapse) > 1:
                    geneticpert.perturbagemovertime(BRhdf5genetic.data, BRhdf5genetic.std, 0, 0, 0)
                    if self.expected_effect == 0:
                        geneticpert.perturbagemovertime(BRhdf5genetic.inhibited, BRhdf5genetic.std_inh, 2, 1, 1)
                    elif self.expected_effect == 1:
                        geneticpert.perturbagemovertime(BRhdf5genetic.normalizedtranslation, BRhdf5genetic.std_inh, 1, 1, 1)

                else:
                    print("1 time point for genetic-chemical perturbagem")
                    if self.expected_effect == 0:
                        geneticpert.perturbagemendpoint(BRhdf5genetic.inhibited, BRhdf5genetic.std_inh, 2, 1)
                    elif self.expected_effect == 1:
                        geneticpert.perturbagemendpoint(BRhdf5genetic.normalizedtranslation, BRhdf5genetic.std_inh, 1, 1)

                geneticpert.close_pdf()

            if catego.experimentype == "genetic-chemical_perturbagen":
                print("genetic-chemical_perturbagen")
                # save information on HDF5 file
                brhdfgeneticchemical = BRgeneticchemicalperturbagem(count, file_br_names.filehdf5resultspath,
                                                                    self.file_name_br, header_info, catego, file_info)

                print(len(brhdfgeneticchemical.fields), len(brhdfgeneticchemical.data), len(brhdfgeneticchemical.std),
                      len(brhdfgeneticchemical.std_inh), len(brhdfgeneticchemical.normalized_perc),
                      len(brhdfgeneticchemical.normalized))
                print("medium information", len(brhdfgeneticchemical.fieldsmedium), len(brhdfgeneticchemical.datamedium),
                      len(brhdfgeneticchemical.stdmedium), len(brhdfgeneticchemical.std_inhmedium),
                      len(brhdfgeneticchemical.normalized_percmedium), len(brhdfgeneticchemical.normalizedtranslationmedium))
                print(brhdfgeneticchemical.possiblecombination)

                geneticchemical = GeneticChemicalPerturbagem(self.synergy_method, brhdfgeneticchemical, file_info,
                                                             file_br_names, information_readout, readout_units,
                                                             self.file_name_br, self.biological_replicate)

                if len(brhdfgeneticchemical.elapse) > 1:
                    geneticchemical.confluencyovertime(brhdfgeneticchemical.data, brhdfgeneticchemical.std, 0, 0)
                    if self.expected_effect == 0:
                        geneticchemical.confluencyovertime(brhdfgeneticchemical.inhibited, brhdfgeneticchemical.std_inh,
                                                           2, 1)
                    elif self.expected_effect == 1:
                        geneticchemical.confluencyovertime(brhdfgeneticchemical.normalizedtranslation,
                                                           brhdfgeneticchemical.std_inh, 1, 1)
                    geneticchemical.doseresponse()
                else:
                    print("are we here?")
                    if self.expected_effect == 0:
                        geneticchemical.endpointinhibition(brhdfgeneticchemical.inhibited, brhdfgeneticchemical.std_inh,
                                                           2, 1)
                    elif self.expected_effect == 1:
                        geneticchemical.endpointinhibition(brhdfgeneticchemical.normalizedtranslation,
                                                           brhdfgeneticchemical.std_inh, 1, 1)
                    geneticchemical.doseresponse(1)

                geneticchemical.close_pdf()

        else:

            for i in self.files_list:
                # creat different file names and paths:
                file_names = Filenames(i, self.input_path, self.results_path, self.information_extracted)
                print("file", i)
                # read input file, one by one, extract:
                # date ,elapsed , header , data ,data_std ,std ,std_header , date_info
                file_info = Readfile(file_names.fileitxtnputpath)
                header_info = Headers(file_info.header)

                # Check if the header has STD or not
                if len(header_info.stdinfo) != 0:
                    header_info.headerswithstd()
                    file_info.get_split_data_std(len(header_info.headersinitial) - len(header_info.stdheader))
                else:
                    header_info.headersnostd()
                    header_info.stdinfo.append("STD computed by HTSplotter")
                    file_info.get_data_nostd()

                catego = Categorisation(header_info)

                out = Outputfile(i, self.biological_replicate, file_names, catego, header_info, file_info.elapsed)

                # set error file
                if len(out.error) == 0:
                    print("HTSplotter did not identify any error from your input file")
                    # out.seterrorfile(0)
                if len(out.error) != 0:
                    print(out.error)
                    out.seterrorfile(1)
                    sys.exit()

                # check user confirmation
                if self.userinput == 1:
                    time_counter = 0
                    while not os.path.exists(self.input_path + "inputfile.txt"):
                        print("still running")
                        time.sleep(1)
                        time_counter += 1
                    if os.path.exists(self.input_path + "inputfile.txt"):
                        userconfi = Inputfile(self.input_path + "inputfile.txt", catego.experimentype)
                        print(catego.experimentype, userconfi.experiment_type)
                        if catego.experimentype != userconfi.experiment_type:
                            out.seterrorfile(0)
                            print("=====>check error file", userconfi.experiment_type, catego.experimentype)
                            break
                        print("=====>ok", userconfi.experiment_type, catego.experimentype)
                        # out.seterrorfile(0)

                # for each experiment type: save information on HDF5 file, process and visualization
                if catego.experimentype == "drug_screen":
                    if len(catego.compound) == len(catego.control):
                        # 1 control for each compound
                        print("drug_screen", len(catego.compound), len(catego.control), catego.medium)
                        # save information on HDF5 file
                        hdfcompound = Compoundscreen(file_names.filehdf5resultspath, header_info.newheader,
                                                     catego, file_info)

                        print(len(hdfcompound.fields), len(hdfcompound.data), len(hdfcompound.std),
                              len(hdfcompound.std_inh), len(hdfcompound.normalized_perc),
                              len(hdfcompound.normalized))

                        print("medium information", len(hdfcompound.datamedium),
                              len(hdfcompound.stdmedium), len(hdfcompound.std_inhmedium),
                              len(hdfcompound.normalized_percmedium), len(hdfcompound.normalizedmedium))

                        singlecompond = SingleCompound(header_info.branch, file_info,
                                                       information_readout, readout_units, self.biological_replicate,
                                                       hdfcompound, file_names, i)

                        if len(hdfcompound.elapse) > 1:
                            singlecompond.inhibitionovertime(hdfcompound.data, hdfcompound.std)
                            if self.expected_effect == 0:
                                singlecompond.inhibitionovertime(hdfcompound.inhibited, hdfcompound.std_inh, 2, 1)
                            elif self.expected_effect == 1:
                                singlecompond.inhibitionovertime(hdfcompound.normalizedtranslation, hdfcompound.std_inh,
                                                                 1, 1)
                            singlecompond.doseresponse()
                        else:
                            singlecompond.doseresponse(1)
                        singlecompond.close_pdf()

                    else:
                        # 1 control for all compounds
                        print("drug_screen", len(catego.compound), len(catego.control))
                        # save information on HDF5 file
                        hdfcompoundonecontrol = Compoundscreenonecontrol(file_names.filehdf5resultspath,
                                                                         header_info.newheader, catego, file_info)

                        print(len(hdfcompoundonecontrol.fields), len(hdfcompoundonecontrol.data),
                              len(hdfcompoundonecontrol.std), len(hdfcompoundonecontrol.std_inh),
                              len(hdfcompoundonecontrol.normalized_perc), len(hdfcompoundonecontrol.normalized))
                        singlecomponecontrol = SingleCompoundonecontrol(header_info.branch, file_info,
                                                                        information_readout, readout_units,
                                                                        self.biological_replicate,
                                                                        hdfcompoundonecontrol, file_names, i)

                        if len(hdfcompoundonecontrol.elapse) > 1:
                            singlecomponecontrol.inhibitiononecontrolovertime(hdfcompoundonecontrol.data,
                                                                              hdfcompoundonecontrol.std)
                            if self.expected_effect == 0:
                                singlecomponecontrol.inhibitiononecontrolovertime(hdfcompoundonecontrol.inhibited,
                                                                                  hdfcompoundonecontrol.std_inh, 2, 1)
                            elif self.expected_effect == 1:
                                singlecomponecontrol.inhibitiononecontrolovertime(
                                    hdfcompoundonecontrol.normalizedtranslation, hdfcompoundonecontrol.std_inh, 1, 1)
                            singlecomponecontrol.doseresponseonecontrol()
                        else:
                            singlecomponecontrol.doseresponseonecontrol(1)
                        singlecomponecontrol.close_pdf()

                if catego.experimentype == "drug_combination":
                    print("drug_combination")
                    print('0-->Bliss; 1-->HSA; 2-->ZIP, the synergy method is: ', self.synergy_method)
                    # save information on HDF5 file
                    hdfcompoundcomb = Compoundcombination(file_names.filehdf5resultspath, header_info.newheader,
                                                          catego, file_info)

                    print(len(hdfcompoundcomb.fields), len(hdfcompoundcomb.data), len(hdfcompoundcomb.std),
                          len(hdfcompoundcomb.std_inh), len(hdfcompoundcomb.normalized_perc),
                          len(hdfcompoundcomb.normalized))

                    comb = ExperimentCombination(self.synergy_method, hdfcompoundcomb, file_info, catego.branch,
                                                 information_readout, readout_units, file_names,
                                                 self.biological_replicate, i)

                    # add predicted and bliss score to the HDF5 file
                    hdfcompoundcomb.add_predictedbiscore(comb.comb_name_per_group,
                                                         comb.synergy_score_per_group,
                                                         comb.effect_per_group,
                                                         self.synergy_method)

                    # # over time data:
                    if len(hdfcompoundcomb.elapse) > 1:
                        comb.inhibition(hdfcompoundcomb.data, hdfcompoundcomb.std)
                        if self.expected_effect == 0:
                            comb.inhibition(hdfcompoundcomb.inhibited, hdfcompoundcomb.std_inh, 2, 1)
                        elif self.expected_effect == 1:
                            comb.inhibition(hdfcompoundcomb.normalizedtranslation, hdfcompoundcomb.std_inh, 1, 1)
                        comb.doseresponse(1)

                    else:
                        print("1 time point for combo")
                        if self.expected_effect == 0:
                            comb.endpointinhibition(hdfcompoundcomb.inhibited, hdfcompoundcomb.std_inh, 2, 1)
                        elif self.expected_effect == 1:
                            comb.endpointinhibition(hdfcompoundcomb.normalizedtranslation, hdfcompoundcomb.std_inh, 1, 1)
                        comb.doseresponse(1)

                    comb.close_pdf()

                if catego.experimentype == "Genetic_perturbagen":
                    print("Genetic_perturbagen")
                    # save information on HDF5 file

                    hdfgenetic = Geneticperturbagem(file_names.filehdf5resultspath, header_info.newheader,
                                                    catego, file_info)

                    print(len(hdfgenetic.fields), len(hdfgenetic.data), len(hdfgenetic.std),
                          len(hdfgenetic.std_inh), len(hdfgenetic.normalized_perc), len(hdfgenetic.normalized))
                    print(hdfgenetic.compoundalone)

                    geneticpert = ExperimentGeneticPerturbagem(header_info, hdfgenetic, file_info, catego, file_names,
                                                               information_readout, readout_units, i,
                                                               self.biological_replicate)

                    if len(hdfgenetic.elapse) > 1:
                        geneticpert.perturbagemovertime(hdfgenetic.data, hdfgenetic.std, 0, 0, 0)
                        if self.expected_effect == 0:
                            geneticpert.perturbagemovertime(hdfgenetic.inhibited, hdfgenetic.std_inh, 2, 1, 1)
                        elif self.expected_effect == 1:
                            geneticpert.perturbagemovertime(hdfgenetic.normalizedtranslation, hdfgenetic.std_inh, 1, 1, 1)

                    else:
                        print("1 time point for genetic perturbagem")
                        if self.expected_effect == 0:
                            geneticpert.perturbagemendpoint(hdfgenetic.inhibited, hdfgenetic.std_inh, 2, 1)
                        elif self.expected_effect == 1:
                            geneticpert.perturbagemendpoint(hdfgenetic.normalizedtranslation, hdfgenetic.std_inh, 1, 1)

                    geneticpert.close_pdf()

                if catego.experimentype == "genetic-chemical_perturbagen":
                    print("genetic-chemical_perturbagen")
                    # save information on HDF5 file
                    hdfgeneticchemical = Geneticchemicalperturbagem(file_names.filehdf5resultspath, header_info.newheader,
                                                                    catego, file_info)

                    print(len(hdfgeneticchemical.fields), len(hdfgeneticchemical.data), len(hdfgeneticchemical.std),
                          len(hdfgeneticchemical.std_inh), len(hdfgeneticchemical.normalized_perc),
                          len(hdfgeneticchemical.normalized))
                    print("medium information", len(hdfgeneticchemical.fieldsmedium), len(hdfgeneticchemical.datamedium),
                          len(hdfgeneticchemical.stdmedium), len(hdfgeneticchemical.std_inhmedium),
                          len(hdfgeneticchemical.normalized_percmedium), len(hdfgeneticchemical.normalizedmedium))
                    #
                    geneticchemical = GeneticChemicalPerturbagem(self.synergy_method, hdfgeneticchemical, file_info,
                                                                 file_names,
                                                                 information_readout, readout_units,
                                                                 i, self.biological_replicate)

                    if len(hdfgeneticchemical.elapse) > 1:
                        geneticchemical.confluencyovertime(hdfgeneticchemical.data, hdfgeneticchemical.std, 0, 0)
                        if self.expected_effect == 0:
                            geneticchemical.confluencyovertime(hdfgeneticchemical.inhibited, hdfgeneticchemical.std_inh,
                                                               2, 1)
                        elif self.expected_effect == 1:
                            geneticchemical.confluencyovertime(hdfgeneticchemical.normalizedtranslation,
                                                               hdfgeneticchemical.std_inh, 1, 1)
                        geneticchemical.doseresponse()

                    else:
                        print("1 time point for genetic-chemical perturbagem")
                        if self.expected_effect == 0:
                            geneticchemical.endpointinhibition(hdfgeneticchemical.inhibited, hdfgeneticchemical.std_inh,
                                                               2, 1)
                        elif self.expected_effect == 1:
                            geneticchemical.endpointinhibition(hdfgeneticchemical.normalizedtranslation,
                                                               hdfgeneticchemical.std_inh, 1, 1)
                        geneticchemical.doseresponse(1)

                    geneticchemical.close_pdf()

    def verify_inputs(self):

        if self.input_path is None:
            self.input_path = self.main_folder.copy()
        if self.information_extracted is None:
            self.information_extracted = self.main_folder.copy()
        if self.results_path is None:
            self.results_path = self.main_folder.copy()

if __name__ == '__main__':
    htsplotter = Analyser()

    htsplotter.main_folder = "HTSplotter/" \
                             "experiment_type/drug_combination/"  # add your folder path e.g. GitHub path

    htsplotter.input_path = htsplotter.main_folder + "inputfile/"  # + input path folder
    htsplotter.information_extracted = htsplotter.main_folder + 'information_extracted_files/'  # + information_extracted path folder
    htsplotter.results_path = htsplotter.main_folder + 'output_results/'  # + output_results path folder

    htsplotter.biological_replicate = 0  # 0 if not, 1 if it is
    htsplotter.userinput = 0  # 1 if yes
    htsplotter.information_readout = "confluency"  # e.g.: "confluency"
    htsplotter.readout_units = "(%)"  # e.g. # (%)
    htsplotter.expected_effect = 0  # 0 = "inhibition"; 1 = enhanced
    htsplotter.file_name_br = "BiologicalReplicateteste"  # file name
    htsplotter.synergy_method = 0  # 0-->Bliss; 1-->HSA; 2--> ZIP

    # GitHub experiment_type folder
    # drug_combination/
    htsplotter.files_list = ['drug_combination_several_time_points_repetitive_conditions',
                             'drug_combination_screen_1timepoint']  # all files to be analyzed
    # drug/
    # htsplotter.files_list = ['drugscreen_1timepoint',
    #                          'drugscreen_severaltimepoint_1control',
    #                          'drugscreen_severaltimepoint_severalcontrol']
    # genetic_pertubagen/
    # htsplotter.files_list = ['gene_perturbagen_1timepoint_1control',
    #                          'gene_perturbagen_severaltimepoints']
    # genetic-chemical_perturbagens/
    # htsplotter.files_list = ['genetic-chemical_perturbagen_1time_point',
    #                          'genetic-chemical_perturbagen_several-time_points']
    # all files to be analyzed
    htsplotter.execute()


