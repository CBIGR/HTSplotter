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
from HTSplotter.plotting import Overtime
from HTSplotter.combination import ExperimentCombination
from HTSplotter.geneticperturbagem import ExperimentGeneticPerturbagem
from HTSplotter.compoundscreen import SingleCompound, SingleCompoundonecontrol
from HTSplotter.geneticchemicalperturbagem import GeneticChemicalPerturbagem

from HTSplotter.txtsavedata import Savetxt


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

        self.verify_inputs()

        start = time.time()

        p = psutil.Process()

        if len(self.information_readout) == 0:
            self.information_readout = "no information"
        if len(self.readout_units) == 0:
            self.readout_units = 'no information'

        if self.biological_replicate == 1:

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

                catego = Categorisation(header_info.newheader, header_info.medium, header_info.controlname,
                                        header_info.stdinfo)

                out = Outputfile(eachfile, self.biological_replicate, file_names.information_extractedfile,
                                 file_names.errorfile, catego.concentration, header_info.stdinfo,
                                 catego.diccompoundgroupkey, catego.diccontrolgroupkey, catego.dicmediumgroupkey,
                                 catego.experimentype, catego.control, catego.compound, catego.combination,
                                 catego.medium,
                                 catego.compoundscreen, catego.condition,
                                 header_info.outputheader, header_info.errorheader)
                # set error file

                if len(out.error) != 0:

                    out.seterrorfile(1)
                    sys.exit()
                # check user confirmation
                if self.userinput == 1:
                    time_counter = 0
                    while not os.path.exists(self.input_path + "inputfile.txt"):

                        time.sleep(1)
                        time_counter += 1
                    if os.path.exists(self.input_path + "inputfile.txt"):
                        userconfi = Inputfile(self.input_path + "inputfile.txt", catego.experimentype)

                        self.information_readout = userconfi.readout

                        if catego.experimentype != userconfi.experiment_type:
                            out.seterrorfile(0)
                            break

                        # out.seterrorfile(0)
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

                catego = Categorisation(header_info.newheader, header_info.medium, header_info.controlname,
                                        header_info.stdinfo)

                if catego.experimentype == "Genetic_perturbagen":
                    indcombdata = Individualgeneticperturbagen(count, file_br_names.filehdf5resultspath,
                                                               catego.experimentype, catego.branch, i,
                                                               header_info.newheader, file_info.elapsed,
                                                               file_info.date_info, file_info.date, file_info.data,
                                                               file_info.std, header_info.medium)
                else:
                    indcombdata = Individualcombinationstructure(count, file_br_names.filehdf5resultspath,
                                                                 catego.experimentype, catego.branch, i,
                                                                 header_info.newheader, file_info.elapsed,
                                                                 file_info.date_info, file_info.date, file_info.data,
                                                                 file_info.std, header_info.medium)
                count += 1

            # # check if the input files have the same compounds and conditions
            if catego.experimentype == "drug_combination":

                brcombdata = BRcombinationstructure(count, file_br_names.filehdf5resultspath, catego.experimentype,
                                                    catego.branch, self.file_name_br, header_info.newheader,
                                                    file_info.elapsed,
                                                    file_info.date_info, file_info.date, file_info.data, file_info.std,
                                                    header_info.medium)

                comb = ExperimentCombination(brcombdata.fields, file_info.elapsed, catego.branch,
                                             brcombdata.control, brcombdata.compoundalone,
                                             brcombdata.possiblecombination,
                                             brcombdata.possiblecombinationsize, brcombdata.inhibited,
                                             brcombdata.std_inh, brcombdata.celline, brcombdata.seeding,
                                             brcombdata.condition, brcombdata.std_info, brcombdata.medium,
                                             self.information_readout, self.readout_units)
                # Save data as txt file
                Savetxt(file_br_names.fileioriginaldatapath, file_info.date_info, brcombdata.fields,
                        brcombdata.data, file_info.elapsed)
                Savetxt(file_br_names.fileinhibiteddatapath, file_info.date_info, brcombdata.fields,
                        brcombdata.inhibited, file_info.elapsed)
                Savetxt(file_br_names.filenormalizedatapath, file_info.date_info, brcombdata.fields,
                        brcombdata.normalized_perc, file_info.elapsed)
                comb.bi_score_save_txt(file_br_names.fileblisscorpath, comb.bi_score_per_group)
                comb.bi_score_save_txt(file_br_names.filepredictedblisscorpath, comb.predicted_per_group)
                # ###

                # add predicted and bliss score to the HDF5 file
                brcombdata.add_predictedbiscore(comb.comb_name_per_group, comb.predicted_per_group,
                                                comb.bi_score_per_group)
                # over time data:
                comb.get_ic_txt_path(file_br_names.fileictxtresultspath)
                comb.open_pdf(file_br_names.filepdfresultspath, self.file_name_br)
                if len(brcombdata.elapse) > 1:
                    comb.confluencyovertime(brcombdata.data, brcombdata.std)
                    if self.expected_effect == 0:
                        comb.inhibition(brcombdata.inhibited, brcombdata.std_inh, 2, 1)
                    elif self.expected_effect == 1:
                        comb.inhibition(brcombdata.normalizedtranslation, brcombdata.std_inh, 1, 1)

                else:

                    if self.expected_effect == 0:
                        comb.endpointinhibition(brcombdata.inhibited, brcombdata.std_inh, 2, 1)
                    if self.expected_effect == 1:
                        comb.endpointinhibition(brcombdata.normalizedtranslation, brcombdata.std_inh, 1, 1)
                comb.close_pdf()

            if catego.experimentype == "drug_screen":

                if len(catego.compound) == len(catego.control):
                    # 1 control for each compound

                    brcomscreendata = BRcompoundscreenseveralcontrol(count, file_br_names.filehdf5resultspath,
                                                                     catego.experimentype, catego.branch, self.file_name_br,
                                                                     header_info.newheader, file_info.elapsed,
                                                                     file_info.date_info, file_info.date,
                                                                     file_info.data,
                                                                     file_info.std, header_info.medium)

                    singlecompond = SingleCompound(header_info.branch, brcomscreendata.celline,
                                                   brcomscreendata.seeding, brcomscreendata.fields,
                                                   file_info.elapsed, brcomscreendata.control,
                                                   brcomscreendata.compoundalone, brcomscreendata.condition,
                                                   brcomscreendata.std_info, brcomscreendata.fieldsmedium,
                                                   brcomscreendata.datamedium, brcomscreendata.stdmedium,
                                                   brcomscreendata.fieldsmediuminhibited,
                                                   brcomscreendata.inhibitedmedium,
                                                   brcomscreendata.std_inhmedium,
                                                   brcomscreendata.normalizedtranslationmedium,
                                                   self.information_readout, self.readout_units)
                    singlecompond.get_ic_txt_path(file_br_names.fileictxtresultspath)
                    singlecompond.open_pdf(file_br_names.filepdfresultspath, self.file_name_br)
                    if len(brcomscreendata.elapse) > 1:
                        singlecompond.confluencyovertime(brcomscreendata.data, brcomscreendata.std)
                        if self.expected_effect == 0:
                            singlecompond.inhibitionovertime(brcomscreendata.inhibited, brcomscreendata.std_inh, 2, 1)
                        elif self.expected_effect == 1:
                            singlecompond.inhibitionovertime(brcomscreendata.normalizedtranslation,
                                                             brcomscreendata.std_inh,
                                                             1, 1)

                    singlecompond.doseresponse(brcomscreendata.normalized_perc, brcomscreendata.std)
                    singlecompond.close_pdf()
                else:

                    brcomscreenonecontrol = BRcompoundscreenonecontrol(count, file_br_names.filehdf5resultspath,
                                                                       catego.experimentype, catego.branch,
                                                                       self.file_name_br,
                                                                       header_info.newheader, file_info.elapsed,
                                                                       file_info.date_info, file_info.date,
                                                                       file_info.data,
                                                                       file_info.std, header_info.medium)

                    singlecomponecontrol = SingleCompoundonecontrol(header_info.branch, brcomscreenonecontrol.celline,
                                                                    brcomscreenonecontrol.seeding,
                                                                    brcomscreenonecontrol.fields,
                                                                    file_info.elapsed, brcomscreenonecontrol.control,
                                                                    brcomscreenonecontrol.compoundalone,
                                                                    brcomscreenonecontrol.condition,
                                                                    brcomscreenonecontrol.std_info,
                                                                    brcomscreenonecontrol.fieldsmedium,
                                                                    brcomscreenonecontrol.datamedium,
                                                                    brcomscreenonecontrol.stdmedium,
                                                                    brcomscreenonecontrol.fieldsmediuminhibited,
                                                                    brcomscreenonecontrol.inhibitedmedium,
                                                                    brcomscreenonecontrol.std_inhmedium,
                                                                    brcomscreenonecontrol.normalizedtranslationmedium,
                                                                    self.information_readout, self.readout_units)

                    singlecomponecontrol.get_ic_txt_path(file_br_names.fileictxtresultspath)
                    singlecomponecontrol.open_pdf(file_br_names.filepdfresultspath, self.file_name_br)
                    if len(brcomscreenonecontrol.elapse) > 1:
                        singlecomponecontrol.confluencyonecontrolovertime(brcomscreenonecontrol.data,
                                                                          brcomscreenonecontrol.std)
                        if self.expected_effect == 0:
                            singlecomponecontrol.inhibitiononecontrolovertime(brcomscreenonecontrol.inhibited,
                                                                              brcomscreenonecontrol.std_inh, 2, 1)
                        elif self.expected_effect == 1:
                            singlecomponecontrol.inhibitiononecontrolovertime(
                                brcomscreenonecontrol.normalizedtranslation,
                                brcomscreenonecontrol.std_inh, 1, 1)

                    singlecomponecontrol.doseresponseonecontrol(brcomscreenonecontrol.normalized_perc,
                                                                brcomscreenonecontrol.std)
                    singlecomponecontrol.close_pdf()

            if catego.experimentype == "Genetic_perturbagen":
                # save information on HDF5 file
                BRhdf5genetic = BRgeneticperturbagem(count, file_br_names.filehdf5resultspath,
                                                     catego.experimentype, catego.branch, self.file_name_br,
                                                     header_info.newheader, file_info.elapsed,
                                                     file_info.date_info, file_info.date, file_info.data,
                                                     file_info.std, header_info.medium)

                geneticpert = ExperimentGeneticPerturbagem(header_info.branch, BRhdf5genetic.celline,
                                                           BRhdf5genetic.seeding,
                                                           BRhdf5genetic.fields, file_info.elapsed,
                                                           BRhdf5genetic.control,
                                                           BRhdf5genetic.compoundalone, BRhdf5genetic.condition,
                                                           BRhdf5genetic.std_info, BRhdf5genetic.fieldsmedium,
                                                           BRhdf5genetic.datamedium,
                                                           BRhdf5genetic.stdmedium, BRhdf5genetic.fieldsmediuminhibited,
                                                           BRhdf5genetic.inhibitedmedium, BRhdf5genetic.std_inhmedium,
                                                           self.information_readout, self.readout_units)
                geneticpert.open_pdf(file_br_names.filepdfresultspath, self.file_name_br)
                if len(BRhdf5genetic.elapse) > 1:
                    geneticpert.perturbagemovertime(BRhdf5genetic.data, BRhdf5genetic.std, 0, 0, 0)
                    if self.expected_effect == 0:
                        geneticpert.perturbagemovertime(BRhdf5genetic.inhibited, BRhdf5genetic.std_inh, 2, 1, 1)
                    elif self.expected_effect == 1:
                        geneticpert.perturbagemovertime(BRhdf5genetic.normalizedtranslation, BRhdf5genetic.std_inh, 1,
                                                        1, 1)

                else:

                    if self.expected_effect == 0:
                        geneticpert.perturbagemendpoint(BRhdf5genetic.inhibited, BRhdf5genetic.std_inh, 2, 1, 1)
                    elif self.expected_effect == 1:
                        geneticpert.perturbagemendpoint(BRhdf5genetic.normalizedtranslation, BRhdf5genetic.std_inh, 1,
                                                        1, 1)

                geneticpert.close_pdf()

            if catego.experimentype == "genetic-chemical_perturbagen":

                # save information on HDF5 file
                brhdfgeneticchemical = BRgeneticchemicalperturbagem(count, file_br_names.filehdf5resultspath,
                                                                    catego.experimentype, catego.branch, self.file_name_br,
                                                                    header_info.newheader, file_info.elapsed,
                                                                    file_info.date_info, file_info.date, file_info.data,
                                                                    file_info.std, header_info.medium)


                geneticchemical = GeneticChemicalPerturbagem(brhdfgeneticchemical.fields, file_info.elapsed,
                                                             brhdfgeneticchemical.branch, brhdfgeneticchemical.control,
                                                             brhdfgeneticchemical.compoundalone,
                                                             brhdfgeneticchemical.possiblecombination,
                                                             brhdfgeneticchemical.possiblecombinationsize,
                                                             brhdfgeneticchemical.inhibited,
                                                             brhdfgeneticchemical.std_inh, brhdfgeneticchemical.celline,
                                                             brhdfgeneticchemical.seeding,
                                                             brhdfgeneticchemical.condition,
                                                             brhdfgeneticchemical.std_info,
                                                             brhdfgeneticchemical.fieldsmedium,
                                                             brhdfgeneticchemical.fieldsmediuminhibited,
                                                             brhdfgeneticchemical.datamedium,
                                                             brhdfgeneticchemical.normalizedtranslationmedium,
                                                             brhdfgeneticchemical.inhibitedmedium,
                                                             brhdfgeneticchemical.stdmedium,
                                                             brhdfgeneticchemical.std_inhmedium,
                                                             self.information_readout, self.readout_units)

                geneticchemical.get_ic_txt_path(file_br_names.fileictxtresultspath)
                geneticchemical.open_pdf(file_br_names.filepdfresultspath, self.file_name_br)
                if len(brhdfgeneticchemical.elapse) > 1:
                    geneticchemical.confluencyovertime(brhdfgeneticchemical.data, brhdfgeneticchemical.std, 0, 0)
                    if self.expected_effect == 0:
                        geneticchemical.confluencyovertime(brhdfgeneticchemical.inhibited, brhdfgeneticchemical.std_inh,
                                                           2, 1)
                    elif self.expected_effect == 1:
                        geneticchemical.confluencyovertime(brhdfgeneticchemical.normalizedtranslation,
                                                           brhdfgeneticchemical.std_inh, 1, 1)

                else:

                    if self.expected_effect == 0:
                        geneticchemical.endpointinhibition(brhdfgeneticchemical.inhibited, brhdfgeneticchemical.std_inh,
                                                           2, 1)
                    elif self.expected_effect == 1:
                        geneticchemical.endpointinhibition(brhdfgeneticchemical.normalizedtranslation,
                                                           brhdfgeneticchemical.std_inh, 1, 1)

                geneticchemical.doseresponse(brhdfgeneticchemical.normalized_perc, brhdfgeneticchemical.std_inh)
                geneticchemical.close_pdf()

        else:

            for i in self.files_list:
                # creat different file names and paths:
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

                catego = Categorisation(header_info.newheader, header_info.medium, header_info.controlname,
                                        header_info.stdinfo)

                out = Outputfile(i, self.biological_replicate, file_names.information_extractedfile, file_names.errorfile,
                                 catego.concentration, header_info.stdinfo, catego.diccompoundgroupkey,
                                 catego.diccontrolgroupkey, catego.dicmediumgroupkey, catego.experimentype,
                                 catego.control, catego.compound, catego.combination, catego.medium,
                                 file_info.elapsed, catego.condition, header_info.outputheader, header_info.errorheader)
                # set error file
                if len(out.error) != 0:
                    out.seterrorfile(1)
                    sys.exit()

                # check user confirmation
                if self.userinput == 1:
                    time_counter = 0
                    while not os.path.exists(self.input_path + "inputfile.txt"):

                        time.sleep(1)
                        time_counter += 1
                    if os.path.exists(self.input_path + "inputfile.txt"):
                        userconfi = Inputfile(self.input_path + "inputfile.txt", catego.experimentype)

                        if catego.experimentype != userconfi.experiment_type:
                            out.seterrorfile(0)
                            break
                        # out.seterrorfile(0)

                # for each experiment type: save information on HDF5 file, process and visualization
                if catego.experimentype == "drug_screen":
                    if len(catego.compound) == len(catego.control):
                        # 1 control for each compound
                        # save information on HDF5 file
                        hdfcompound = Compoundscreen(catego.experimentype, catego.branch,
                                                     file_names.filehdf5resultspath,
                                                     header_info.newheader, file_info.elapsed, file_info.date_info,
                                                     file_info.date, file_info.data, file_info.std, catego.stdinfo,
                                                     catego.medium, catego.compound)

                        singlecompond = SingleCompound(header_info.branch, hdfcompound.celline,
                                                       hdfcompound.seeding, hdfcompound.fields,
                                                       file_info.elapsed, hdfcompound.control,
                                                       hdfcompound.compoundalone, hdfcompound.condition,
                                                       catego.stdinfo, hdfcompound.fieldsmedium,
                                                       hdfcompound.datamedium, hdfcompound.stdmedium,
                                                       hdfcompound.fieldsmediuminhibited,
                                                       hdfcompound.inhibitedmedium, hdfcompound.std_inhmedium,
                                                       hdfcompound.normalizedtranslationmedium,
                                                       self.information_readout, self.readout_units)

                        singlecompond.get_ic_txt_path(file_names.fileictxtresultspath)
                        singlecompond.open_pdf(file_names.filepdfresultspath, i)
                        if len(hdfcompound.elapse) > 1:
                            singlecompond.confluencyovertime(hdfcompound.data, hdfcompound.std)
                            if self.expected_effect == 0:
                                singlecompond.inhibitionovertime(hdfcompound.inhibited, hdfcompound.std_inh, 2, 1)
                            elif self.expected_effect == 1:
                                singlecompond.inhibitionovertime(hdfcompound.normalizedtranslation, hdfcompound.std_inh,
                                                                 1, 1)

                        singlecompond.doseresponse(hdfcompound.normalized_perc, hdfcompound.std)

                        singlecompond.close_pdf()

                    else:
                        # 1 control for all compounds
                        # save information on HDF5 file
                        hdfcompoundonecontrol = Compoundscreenonecontrol(catego.experimentype, catego.branch,
                                                                         file_names.filehdf5resultspath,
                                                                         header_info.newheader, file_info.elapsed,
                                                                         file_info.date_info, file_info.date,
                                                                         file_info.data, file_info.std,
                                                                         catego.stdinfo, catego.medium, catego.compound)

                        singlecomponecontrol = SingleCompoundonecontrol(header_info.branch,
                                                                        hdfcompoundonecontrol.celline,
                                                                        hdfcompoundonecontrol.seeding,
                                                                        hdfcompoundonecontrol.fields,
                                                                        file_info.elapsed,
                                                                        hdfcompoundonecontrol.control,
                                                                        hdfcompoundonecontrol.compoundalone,
                                                                        hdfcompoundonecontrol.condition,
                                                                        catego.stdinfo,
                                                                        hdfcompoundonecontrol.fieldsmedium,
                                                                        hdfcompoundonecontrol.datamedium,
                                                                        hdfcompoundonecontrol.stdmedium,
                                                                        hdfcompoundonecontrol.fieldsmediuminhibited,
                                                                        hdfcompoundonecontrol.inhibitedmedium,
                                                                        hdfcompoundonecontrol.normalizedtranslationmedium,
                                                                        hdfcompoundonecontrol.std_inhmedium,
                                                                        self.information_readout, self.readout_units)

                        singlecomponecontrol.get_ic_txt_path(file_names.fileictxtresultspath)
                        singlecomponecontrol.open_pdf(file_names.filepdfresultspath, i)
                        if len(hdfcompoundonecontrol.elapse) > 1:
                            singlecomponecontrol.confluencyonecontrolovertime(hdfcompoundonecontrol.data,
                                                                              hdfcompoundonecontrol.std)
                            if self.expected_effect == 0:
                                singlecomponecontrol.inhibitiononecontrolovertime(hdfcompoundonecontrol.inhibited,
                                                                                  hdfcompoundonecontrol.std_inh, 2, 1)
                            elif self.expected_effect == 1:
                                singlecomponecontrol.inhibitiononecontrolovertime(
                                    hdfcompoundonecontrol.normalizedtranslation, hdfcompoundonecontrol.std_inh, 1, 1)

                        singlecomponecontrol.doseresponseonecontrol(hdfcompoundonecontrol.normalized_perc,
                                                                    hdfcompoundonecontrol.std)

                        singlecomponecontrol.close_pdf()

                if catego.experimentype == "drug_combination":
                    # save information on HDF5 file
                    hdfcompoundcomb = Compoundcombination(catego.experimentype, catego.branch,
                                                          file_names.filehdf5resultspath, header_info.newheader,
                                                          file_info.elapsed, file_info.date_info, file_info.date,
                                                          file_info.data, file_info.std, catego.stdinfo,
                                                          catego.medium, catego.compound)

                    comb = ExperimentCombination(hdfcompoundcomb.fields, file_info.elapsed, catego.branch,
                                                 hdfcompoundcomb.control, hdfcompoundcomb.compoundalone,
                                                 hdfcompoundcomb.possiblecombination,
                                                 hdfcompoundcomb.possiblecombinationsize, hdfcompoundcomb.inhibited,
                                                 hdfcompoundcomb.std_inh, hdfcompoundcomb.celline,
                                                 hdfcompoundcomb.seeding,
                                                 hdfcompoundcomb.condition, hdfcompoundcomb.std_info,
                                                 hdfcompoundcomb.medium, self.information_readout, self.readout_units)
                    # Save data as txt file
                    Savetxt(file_names.fileioriginaldatapath, file_info.date_info, hdfcompoundcomb.fields,
                            hdfcompoundcomb.data, file_info.elapsed)
                    Savetxt(file_names.fileinhibiteddatapath, file_info.date_info, hdfcompoundcomb.fields,
                            hdfcompoundcomb.inhibited, file_info.elapsed)
                    Savetxt(file_names.filenormalizedatapath, file_info.date_info, hdfcompoundcomb.fields,
                            hdfcompoundcomb.normalized_perc, file_info.elapsed)
                    comb.bi_score_save_txt(file_names.fileblisscorpath, comb.bi_score_per_group)
                    comb.bi_score_save_txt(file_names.filepredictedblisscorpath, comb.predicted_per_group)
                    # ###

                    # add predicted and bliss score to the HDF5 file
                    hdfcompoundcomb.add_predictedbiscore(comb.comb_name_per_group,
                                                         comb.bi_score_per_group, comb.predicted_per_group)
                    comb.get_ic_txt_path(file_names.fileictxtresultspath)
                    comb.open_pdf(file_names.filepdfresultspath, i)
                    # over time data:
                    if len(hdfcompoundcomb.elapse) > 1:
                        #
                        comb.confluencyovertime(hdfcompoundcomb.data, hdfcompoundcomb.std)
                        if self.expected_effect == 0:
                            comb.inhibition(hdfcompoundcomb.inhibited, hdfcompoundcomb.std_inh, 2, 1)
                        elif self.expected_effect == 1:
                            comb.inhibition(hdfcompoundcomb.normalizedtranslation, hdfcompoundcomb.std_inh, 1, 1)

                    else:
                        if self.expected_effect == 0:
                            comb.endpointinhibition(hdfcompoundcomb.inhibited, hdfcompoundcomb.std_inh, 2, 1)
                        elif self.expected_effect == 1:
                            comb.endpointinhibition(hdfcompoundcomb.normalizedtranslation, hdfcompoundcomb.std_inh, 1,
                                                    1)

                    comb.doseresponse(hdfcompoundcomb.normalized_perc, hdfcompoundcomb.std)
                    comb.close_pdf()

                if catego.experimentype == "Genetic_perturbagen":
                    # save information on HDF5 file
                    hdfgenetic = Geneticperturbagem(catego.experimentype, catego.branch, file_names.filehdf5resultspath,
                                                    header_info.newheader, file_info.elapsed, file_info.date_info,
                                                    file_info.date, file_info.data, file_info.std, catego.stdinfo,
                                                    catego.medium, catego.compound)

                    geneticpert = ExperimentGeneticPerturbagem(header_info.branch, hdfgenetic.celline,
                                                               hdfgenetic.seeding,
                                                               hdfgenetic.fields, file_info.elapsed, hdfgenetic.control,
                                                               hdfgenetic.compoundalone, hdfgenetic.condition,
                                                               catego.stdinfo, hdfgenetic.fieldsmedium,
                                                               hdfgenetic.datamedium, hdfgenetic.stdmedium,
                                                               hdfgenetic.fieldsmediuminhibited,
                                                               hdfgenetic.inhibitedmedium, hdfgenetic.std_inhmedium,
                                                               self.information_readout, self.readout_units)
                    Savetxt(file_names.fileioriginaldatapath, file_info.date_info, hdfgenetic.fields,
                            hdfgenetic.data, file_info.elapsed)

                    geneticpert.open_pdf(file_names.filepdfresultspath, i)
                    if len(hdfgenetic.elapse) > 1:
                        geneticpert.perturbagemovertime(hdfgenetic.data, hdfgenetic.std, 0, 0, 0)
                        if self.expected_effect == 0:
                            geneticpert.perturbagemovertime(hdfgenetic.inhibited, hdfgenetic.std_inh, 2, 1, 1)
                        elif self.expected_effect == 1:
                            geneticpert.perturbagemovertime(hdfgenetic.normalizedtranslation, hdfgenetic.std_inh, 1, 1,
                                                            1)

                    else:
                        if self.expected_effect == 0:
                            geneticpert.perturbagemendpoint(hdfgenetic.inhibited, hdfgenetic.std_inh, 2,
                                                            1, 1)
                        elif self.expected_effect == 1:
                            geneticpert.perturbagemendpoint(hdfgenetic.normalizedtranslation, hdfgenetic.std_inh, 1, 1,
                                                            1)

                    geneticpert.close_pdf()

                if catego.experimentype == "genetic-chemical_perturbagen":
                    # save information on HDF5 file
                    hdfgeneticchemical = Geneticchemicalperturbagem(catego.experimentype, catego.branch,
                                                                    file_names.filehdf5resultspath,
                                                                    header_info.newheader,
                                                                    file_info.elapsed, file_info.date_info,
                                                                    file_info.date,
                                                                    file_info.data, file_info.std, catego.stdinfo,
                                                                    catego.medium, catego.compound)
                    #
                    geneticchemical = GeneticChemicalPerturbagem(hdfgeneticchemical.fields, file_info.elapsed,
                                                                 hdfgeneticchemical.branch, hdfgeneticchemical.control,
                                                                 hdfgeneticchemical.compoundalone,
                                                                 hdfgeneticchemical.possiblecombination,
                                                                 hdfgeneticchemical.possiblecombinationsize,
                                                                 hdfgeneticchemical.inhibited,
                                                                 hdfgeneticchemical.std_inh, hdfgeneticchemical.celline,
                                                                 hdfgeneticchemical.seeding,
                                                                 hdfgeneticchemical.condition,
                                                                 hdfgeneticchemical.std_info,
                                                                 hdfgeneticchemical.fieldsmedium,
                                                                 hdfgeneticchemical.fieldsmediuminhibited,
                                                                 hdfgeneticchemical.datamedium,
                                                                 hdfgeneticchemical.normalizedtranslationmedium,
                                                                 hdfgeneticchemical.inhibitedmedium,
                                                                 hdfgeneticchemical.stdmedium,
                                                                 hdfgeneticchemical.std_inhmedium,
                                                                 self.information_readout, self.readout_units)

                    geneticchemical.get_ic_txt_path(file_names.fileictxtresultspath)
                    Savetxt(file_names.fileioriginaldatapath, file_info.date_info, hdfgeneticchemical.fields,
                            hdfgeneticchemical.data, file_info.elapsed)
                    Savetxt(file_names.fileinhibiteddatapath, file_info.date_info, hdfgeneticchemical.fields,
                            hdfgeneticchemical.inhibited, file_info.elapsed)
                    Savetxt(file_names.filenormalizedatapath, file_info.date_info, hdfgeneticchemical.fields,
                            hdfgeneticchemical.normalizedtranslation, file_info.elapsed)

                    geneticchemical.bi_score_save_txt(file_names.fileblisscorpath, geneticchemical.bi_score_per_group)
                    geneticchemical.bi_score_save_txt(file_names.filepredictedblisscorpath,
                                                      geneticchemical.predicted_per_group)
                    geneticchemical.open_pdf(file_names.filepdfresultspath, i)
                    if len(hdfgeneticchemical.elapse) > 1:
                        geneticchemical.confluencyovertime(hdfgeneticchemical.data, hdfgeneticchemical.std, 0, 0)
                        if self.expected_effect == 0:
                            geneticchemical.confluencyovertime(hdfgeneticchemical.inhibited, hdfgeneticchemical.std_inh,
                                                               2, 1)
                        elif self.expected_effect == 1:
                            geneticchemical.confluencyovertime(hdfgeneticchemical.normalizedtranslation,
                                                               hdfgeneticchemical.std_inh, 1, 1)

                    else:
                        if self.expected_effect == 0:
                            geneticchemical.endpointinhibition(hdfgeneticchemical.inhibited, hdfgeneticchemical.std_inh,
                                                               2, 1)
                        elif self.expected_effect == 1:
                            geneticchemical.endpointinhibition(hdfgeneticchemical.normalizedtranslation,
                                                               hdfgeneticchemical.std_inh, 1, 1)

                    geneticchemical.doseresponse(hdfgeneticchemical.normalized_perc, hdfgeneticchemical.std_inh)

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

    htsplotter.main_folder = "HTSplotter/"\
                             "experiment_type/drug_combination/"# add your folder path e.g. GitHub path

    htsplotter.input_path = htsplotter.main_folder + "inputfile/"# + input path folder
    htsplotter.information_extracted = htsplotter.main_folder + 'information_extracted_files/'# + information_extracted path folder
    htsplotter.results_path = htsplotter.main_folder + 'output_results/'# + output_results path folder

    htsplotter.biological_replicate = 0  # 0 if not, 1 if it is
    htsplotter.userinput = 0  # 1 if yes
    htsplotter.information_readout = "confluency"  # e.g.: "confluency"
    htsplotter.readout_units = "(%)" # e.g. # (%)
    htsplotter.expected_effect = 0  # 0 = "inhibition"; 1 = enhanced
    htsplotter.file_name_br = "BiologicalReplicateteste" # file name

    #GitHub experiment_type folder
    #drug_combination/
    htsplotter.files_list = ['drug_combination_several_time_points_repetitive_conditions',
                             'drug_combination_screen_1timepoint'] # all files to be analyzed
    #drug/
    # htsplotter.files_list = ['drugscreen_1timepoint',
    #                          'drugscreen_severaltimepoint_1control',
    #                          'drugscreen_severaltimepoint_severalcontrol']
    #genetic_pertubagen/
    # htsplotter.files_list = ['gene_perturbagen_1timepoint_1control',
    #                          'gene_perturbagen_severaltimepoints']
    # genetic-chemical_perturbagens/
    # htsplotter.files_list = ['genetic-chemical_perturbagen_1time_point',
    #                          'genetic-chemical_perturbagen_several-time_points']
    # all files to be analyzed
    htsplotter.execute()

