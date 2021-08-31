import os
import time
import psutil
import datetime
import sys
import numpy as np
# librarys made for this script
from readfiles import Readfile
from filenames import Filenames
from headers import Headers
from categorisation import Categorisation
from interactionsuser import Outputfile, Inputfile
from save_hdf5file import Compoundscreenonecontrol, Compoundscreen, \
    Compoundcombination, Geneticperturbagem, Geneticchemicalperturbagem
from save_hdf5brfiles import BRHdf5database, Individualcombinationstructure, BRcompoundscreenseveralcontrol, \
    BRcombinationstructure, BRcompoundscreenonecontrol, BRgeneticperturbagem, BRgeneticchemicalperturbagem, \
    Individualgeneticperturbagen
from plotting import Overtime
from combination import ExperimentCombination
from geneticperturbagem import ExperimentGeneticPerturbagem
from compoundscreen import SingleCompound, SingleCompoundonecontrol
from geneticchemicalperturbagem import GeneticChemicalPerturbagem

from txtsavedata import Savetxt

if __name__ == '__main__':
    start = time.time()
    # Number of logical CPUs in the system
    p = psutil.Process()
    # print("psutil.cpu_count() = {0}".format(psutil.cpu_count()))
    # main folder where csv files are
    main_folder = "C:/Users/cdcarval/Dropbox (speleman lab)/Personal Lab/HTSplotter/experiment_type/drug_combination/"
    input_path = main_folder + "Inputfile/"
    information_extracted = main_folder + "Information_extracted_files/"
    results_path = main_folder + "Output_results/"

    biological_replicate = 0  # 0 if not, 1 if it is
    userinput = 0  # 1 if yes
    information_readout = "confluency"  # default = confluency = 0; add effect name = 1
    readout_units = "(%)"
    expected_effect = 0  # 0 = "inhibition"; 1 = enhanced
    file_name_br = "BiologicalReplicateteste"

    files_list = ['drug_combination_several_time_points']

    print("main folder=>", main_folder)
    print("files path", input_path)
    print("results path", results_path)
    print("files list", files_list)
    if len(information_readout) == 0:
        information_readout = "no information"
    if len(readout_units) == 0:
        readout_units = 'no information'

    if biological_replicate == 1:
        print("Processing Bioloical replicate")
        count = 0
        file_br_names = Filenames(file_name_br, input_path, results_path, information_extracted)
        # first go over each file before processing it
        # Check confirmation to the user for all files
        for eachfile in files_list:
            file_names = Filenames(eachfile, input_path, results_path, information_extracted)
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

            out = Outputfile(eachfile, biological_replicate, file_names.information_extractedfile,
                             file_names.errorfile, catego.concentration, header_info.stdinfo,
                             catego.diccompoundgroupkey, catego.diccontrolgroupkey, catego.dicmediumgroupkey,
                             catego.experimentype, catego.control, catego.compound, catego.combination, catego.medium,
                             catego.compoundscreen, catego.condition,
                             header_info.outputheader, header_info.errorheader)
            # set error file
            if len(out.error) == 0:
                print("HTSplotter did not identify any error from your input file")
                #out.seterrorfile(0)
            if len(out.error) != 0:
                print(out.error)
                out.seterrorfile(1)
                sys.exit()
            # check user confirmation
            if userinput == 1:
                time_counter = 0
                while not os.path.exists(input_path + "inputfile.txt"):
                    print("still running, BR")
                    time.sleep(1)
                    time_counter += 1
                if os.path.exists(input_path + "inputfile.txt"):
                    userconfi = Inputfile(input_path + "inputfile.txt", catego.experimentype)
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
        for i in files_list:
            file_names = Filenames(i, input_path, results_path, information_extracted)
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
            print("Compound_combinationBR", catego.experimentype)
            brcombdata = BRcombinationstructure(count, file_br_names.filehdf5resultspath, catego.experimentype,
                                                catego.branch, file_name_br, header_info.newheader, file_info.elapsed,
                                                file_info.date_info, file_info.date, file_info.data, file_info.std,
                                                header_info.medium)

            print(len(brcombdata.fields), len(brcombdata.data), len(brcombdata.std), len(brcombdata.std_inh),
                  len(brcombdata.normalized_perc), len(brcombdata.normalized))
            for g in range(len(brcombdata.possiblecombination)):
                print(brcombdata.possiblecombinationsize[g], brcombdata.possiblecombination[g])

            comb = ExperimentCombination(brcombdata.fields, file_info.elapsed, catego.branch,
                                         brcombdata.control, brcombdata.compoundalone,
                                         brcombdata.possiblecombination,
                                         brcombdata.possiblecombinationsize, brcombdata.inhibited,
                                         brcombdata.std_inh, brcombdata.celline, brcombdata.seeding,
                                         brcombdata.condition, brcombdata.std_info, brcombdata.medium,
                                         information_readout, readout_units)
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
            brcombdata.add_predictedbiscore(comb.comb_name_per_group, comb.predicted_per_group, comb.bi_score_per_group)
            # over time data:
            comb.get_ic_txt_path(file_br_names.fileictxtresultspath)
            comb.open_pdf(file_br_names.filepdfresultspath, file_name_br)
            if len(brcombdata.elapse) > 1:
                comb.confluencyovertime(brcombdata.data, brcombdata.std)
                if expected_effect == 0:
                    comb.inhibition(brcombdata.inhibited, brcombdata.std_inh, 2, 1)
                elif expected_effect == 1:
                    comb.inhibition(brcombdata.normalizedtranslation, brcombdata.std_inh, 1, 1)

            else:
                print("1 time point for combination biological replicates")
                if expected_effect == 0:
                    comb.endpointinhibition(brcombdata.inhibited, brcombdata.std_inh, 2, 1)
                if expected_effect == 1:
                    comb.endpointinhibition(brcombdata.normalizedtranslation, brcombdata.std_inh, 1, 1)
            comb.close_pdf()

        if catego.experimentype == "drug_screen":
            print("Compound_screen", catego.experimentype)
            if len(catego.compound) == len(catego.control):
                # 1 control for each compound
                print("Compound_screen several control BR", len(catego.compound), len(catego.control), catego.medium)
                brcomscreendata = BRcompoundscreenseveralcontrol(count, file_br_names.filehdf5resultspath,
                                                                 catego.experimentype, catego.branch, file_name_br,
                                                                 header_info.newheader, file_info.elapsed,
                                                                 file_info.date_info, file_info.date, file_info.data,
                                                                 file_info.std, header_info.medium)
                print(len(brcomscreendata.fields), len(brcomscreendata.data), len(brcomscreendata.std),
                      len(brcomscreendata.std_inh), len(brcomscreendata.normalized_perc),
                      len(brcomscreendata.normalized))
                for g in range(len(brcomscreendata.compoundalone)):
                    print(brcomscreendata.compoundalone[g])
                print(brcomscreendata.medium)
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
                                               information_readout, readout_units)
                singlecompond.get_ic_txt_path(file_br_names.fileictxtresultspath)
                singlecompond.open_pdf(file_br_names.filepdfresultspath, file_name_br)
                if len(brcomscreendata.elapse) > 1:
                    singlecompond.confluencyovertime(brcomscreendata.data, brcomscreendata.std)
                    if expected_effect == 0:
                        singlecompond.inhibitionovertime(brcomscreendata.inhibited, brcomscreendata.std_inh, 2, 1)
                    elif expected_effect == 1:
                        singlecompond.inhibitionovertime(brcomscreendata.normalizedtranslation, brcomscreendata.std_inh,
                                                         1, 1)

                singlecompond.doseresponse(brcomscreendata.normalized_perc, brcomscreendata.std)
                singlecompond.close_pdf()
            else:
                print("drug_screen 1 control BR", len(catego.compound), len(catego.control), catego.medium)
                brcomscreenonecontrol = BRcompoundscreenonecontrol(count, file_br_names.filehdf5resultspath,
                                                                   catego.experimentype, catego.branch, file_name_br,
                                                                   header_info.newheader, file_info.elapsed,
                                                                   file_info.date_info, file_info.date, file_info.data,
                                                                   file_info.std, header_info.medium)
                print("one control")
                print(len(brcomscreenonecontrol.fields), len(brcomscreenonecontrol.data),
                      len(brcomscreenonecontrol.std), len(brcomscreenonecontrol.std_inh),
                      len(brcomscreenonecontrol.normalized_perc), len(brcomscreenonecontrol.normalized))
                for g in range(len(brcomscreenonecontrol.compoundalone)):
                    print(brcomscreenonecontrol.compoundalone[g])
                print(brcomscreenonecontrol.medium)

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
                                                                information_readout, readout_units)

                singlecomponecontrol.get_ic_txt_path(file_br_names.fileictxtresultspath)
                singlecomponecontrol.open_pdf(file_br_names.filepdfresultspath, file_name_br)
                if len(brcomscreenonecontrol.elapse) > 1:
                    singlecomponecontrol.confluencyonecontrolovertime(brcomscreenonecontrol.data,
                                                                      brcomscreenonecontrol.std)
                    if expected_effect == 0:
                        singlecomponecontrol.inhibitiononecontrolovertime(brcomscreenonecontrol.inhibited,
                                                                          brcomscreenonecontrol.std_inh, 2, 1)
                    elif expected_effect == 1:
                        singlecomponecontrol.inhibitiononecontrolovertime(brcomscreenonecontrol.normalizedtranslation,
                                                                          brcomscreenonecontrol.std_inh, 1, 1)

                singlecomponecontrol.doseresponseonecontrol(brcomscreenonecontrol.normalized_perc,
                                                            brcomscreenonecontrol.std)
                singlecomponecontrol.close_pdf()

        if catego.experimentype == "Genetic_perturbagen":
            print("Genetic_perturbagen")
            # save information on HDF5 file
            BRhdf5genetic = BRgeneticperturbagem(count, file_br_names.filehdf5resultspath,
                                                 catego.experimentype, catego.branch, file_name_br,
                                                 header_info.newheader, file_info.elapsed,
                                                 file_info.date_info, file_info.date, file_info.data,
                                                 file_info.std, header_info.medium)
            print("did it?>>>")
            print(len(BRhdf5genetic.fields), len(BRhdf5genetic.data), len(BRhdf5genetic.std),
                  len(BRhdf5genetic.std_inh), len(BRhdf5genetic.normalized_perc), len(BRhdf5genetic.normalized))
            print(BRhdf5genetic.compoundalone)
            geneticpert = ExperimentGeneticPerturbagem(header_info.branch, BRhdf5genetic.celline, BRhdf5genetic.seeding,
                                                       BRhdf5genetic.fields, file_info.elapsed, BRhdf5genetic.control,
                                                       BRhdf5genetic.compoundalone, BRhdf5genetic.condition,
                                                       BRhdf5genetic.std_info, BRhdf5genetic.fieldsmedium,
                                                       BRhdf5genetic.datamedium,
                                                       BRhdf5genetic.stdmedium, BRhdf5genetic.fieldsmediuminhibited,
                                                       BRhdf5genetic.inhibitedmedium, BRhdf5genetic.std_inhmedium,
                                                       information_readout, readout_units)
            geneticpert.open_pdf(file_br_names.filepdfresultspath, file_name_br)
            if len(BRhdf5genetic.elapse) > 1:
                geneticpert.perturbagemovertime(BRhdf5genetic.data, BRhdf5genetic.std, 0, 0, 0)
                if expected_effect == 0:
                    geneticpert.perturbagemovertime(BRhdf5genetic.inhibited, BRhdf5genetic.std_inh, 2, 1, 1)
                elif expected_effect == 1:
                    geneticpert.perturbagemovertime(BRhdf5genetic.normalizedtranslation, BRhdf5genetic.std_inh, 1, 1, 1)

            else:
                print("1 time point for genetic-chemical perturbagem")
                if expected_effect == 0:
                    geneticpert.perturbagemendpoint(BRhdf5genetic.inhibited, BRhdf5genetic.std_inh, 2, 1, 1)
                elif expected_effect == 1:
                    geneticpert.perturbagemendpoint(BRhdf5genetic.normalizedtranslation, BRhdf5genetic.std_inh, 1, 1, 1)


            geneticpert.close_pdf()

        if catego.experimentype == "genetic-chemical_perturbagen":
            print("genetic-chemical_perturbagen")
            # save information on HDF5 file
            brhdfgeneticchemical = BRgeneticchemicalperturbagem(count, file_br_names.filehdf5resultspath,
                                                                catego.experimentype, catego.branch, file_name_br,
                                                                header_info.newheader, file_info.elapsed,
                                                                file_info.date_info, file_info.date, file_info.data,
                                                                file_info.std, header_info.medium)

            print(len(brhdfgeneticchemical.fields), len(brhdfgeneticchemical.data), len(brhdfgeneticchemical.std),
                  len(brhdfgeneticchemical.std_inh), len(brhdfgeneticchemical.normalized_perc),
                  len(brhdfgeneticchemical.normalized))
            print("medium information", len(brhdfgeneticchemical.fieldsmedium), len(brhdfgeneticchemical.datamedium),
                  len(brhdfgeneticchemical.stdmedium), len(brhdfgeneticchemical.std_inhmedium),
                  len(brhdfgeneticchemical.normalized_percmedium), len(brhdfgeneticchemical.normalizedtranslationmedium))
            print(brhdfgeneticchemical.possiblecombination)
            # print(brhdfgeneticchemical.branch)
            geneticchemical = GeneticChemicalPerturbagem(brhdfgeneticchemical.fields, file_info.elapsed,
                                                         brhdfgeneticchemical.branch, brhdfgeneticchemical.control,
                                                         brhdfgeneticchemical.compoundalone,
                                                         brhdfgeneticchemical.possiblecombination,
                                                         brhdfgeneticchemical.possiblecombinationsize,
                                                         brhdfgeneticchemical.inhibited,
                                                         brhdfgeneticchemical.std_inh, brhdfgeneticchemical.celline,
                                                         brhdfgeneticchemical.seeding, brhdfgeneticchemical.condition,
                                                         brhdfgeneticchemical.std_info,
                                                         brhdfgeneticchemical.fieldsmedium,
                                                         brhdfgeneticchemical.fieldsmediuminhibited,
                                                         brhdfgeneticchemical.datamedium,
                                                         brhdfgeneticchemical.normalizedtranslationmedium,
                                                         brhdfgeneticchemical.inhibitedmedium,
                                                         brhdfgeneticchemical.stdmedium,
                                                         brhdfgeneticchemical.std_inhmedium,
                                                         information_readout, readout_units)

            geneticchemical.get_ic_txt_path(file_br_names.fileictxtresultspath)
            geneticchemical.open_pdf(file_br_names.filepdfresultspath, file_name_br)
            if len(brhdfgeneticchemical.elapse) > 1:
                geneticchemical.confluencyovertime(brhdfgeneticchemical.data, brhdfgeneticchemical.std, 0, 0)
                if expected_effect == 0:
                    geneticchemical.confluencyovertime(brhdfgeneticchemical.inhibited, brhdfgeneticchemical.std_inh,
                                                       2, 1)
                elif expected_effect == 1:
                    geneticchemical.confluencyovertime(brhdfgeneticchemical.normalizedtranslation,
                                                       brhdfgeneticchemical.std_inh, 1, 1)

            else:
                print("are we here?")
                if expected_effect == 0:
                    geneticchemical.endpointinhibition(brhdfgeneticchemical.inhibited, brhdfgeneticchemical.std_inh,
                                                       2, 1)
                elif expected_effect == 1:
                    geneticchemical.endpointinhibition(brhdfgeneticchemical.normalizedtranslation,
                                                       brhdfgeneticchemical.std_inh, 1, 1)

            geneticchemical.doseresponse(brhdfgeneticchemical.normalized_perc, brhdfgeneticchemical.std_inh)
            geneticchemical.close_pdf()

    else:

        for i in files_list:
            # creat different file names and paths:
            file_names = Filenames(i, input_path, results_path, information_extracted)

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

            out = Outputfile(i, biological_replicate, file_names.information_extractedfile, file_names.errorfile,
                             catego.concentration, header_info.stdinfo, catego.diccompoundgroupkey,
                             catego.diccontrolgroupkey, catego.dicmediumgroupkey, catego.experimentype,
                             catego.control, catego.compound, catego.combination, catego.medium,
                             file_info.elapsed, catego.condition, header_info.outputheader, header_info.errorheader)
            # set error file
            if len(out.error) == 0:
                print("HTSplotter did not identify any error from your input file")
                # out.seterrorfile(0)
            if len(out.error) != 0:
                print(out.error)
                out.seterrorfile(1)
                sys.exit()

            # check user confirmation
            if userinput == 1:
                time_counter = 0
                while not os.path.exists(input_path + "inputfile.txt"):
                    print("still running")
                    time.sleep(1)
                    time_counter += 1
                if os.path.exists(input_path + "inputfile.txt"):
                    userconfi = Inputfile(input_path + "inputfile.txt", catego.experimentype)
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
                    hdfcompound = Compoundscreen(catego.experimentype, catego.branch, file_names.filehdf5resultspath,
                                                 header_info.newheader, file_info.elapsed, file_info.date_info,
                                                 file_info.date, file_info.data, file_info.std, catego.stdinfo,
                                                 catego.medium, catego.compound)
                    print(len(hdfcompound.fields), len(hdfcompound.data), len(hdfcompound.std),
                          len(hdfcompound.std_inh), len(hdfcompound.normalized_perc),
                          len(hdfcompound.normalized))

                    print("medium information", len(hdfcompound.datamedium),
                          len(hdfcompound.stdmedium), len(hdfcompound.std_inhmedium),
                          len(hdfcompound.normalized_percmedium), len(hdfcompound.normalizedmedium))

                    singlecompond = SingleCompound(header_info.branch, hdfcompound.celline,
                                                   hdfcompound.seeding, hdfcompound.fields,
                                                   file_info.elapsed, hdfcompound.control,
                                                   hdfcompound.compoundalone, hdfcompound.condition,
                                                   catego.stdinfo, hdfcompound.fieldsmedium,
                                                   hdfcompound.datamedium, hdfcompound.stdmedium,
                                                   hdfcompound.fieldsmediuminhibited,
                                                   hdfcompound.inhibitedmedium, hdfcompound.std_inhmedium,
                                                   hdfcompound.normalizedtranslationmedium,
                                                   information_readout, readout_units)

                    singlecompond.get_ic_txt_path(file_names.fileictxtresultspath)
                    singlecompond.open_pdf(file_names.filepdfresultspath, i)
                    if len(hdfcompound.elapse) > 1:
                        singlecompond.confluencyovertime(hdfcompound.data, hdfcompound.std)
                        if expected_effect == 0:
                            singlecompond.inhibitionovertime(hdfcompound.inhibited, hdfcompound.std_inh, 2, 1)
                        elif expected_effect == 1:
                            singlecompond.inhibitionovertime(hdfcompound.normalizedtranslation, hdfcompound.std_inh,
                                                             1, 1)

                    singlecompond.doseresponse(hdfcompound.normalized_perc, hdfcompound.std)

                    singlecompond.close_pdf()

                else:
                    # 1 control for all compounds
                    print("drug_screen", len(catego.compound), len(catego.control))
                    # save information on HDF5 file
                    hdfcompoundonecontrol = Compoundscreenonecontrol(catego.experimentype, catego.branch,
                                                                     file_names.filehdf5resultspath,
                                                                     header_info.newheader, file_info.elapsed,
                                                                     file_info.date_info, file_info.date,
                                                                     file_info.data, file_info.std,
                                                                     catego.stdinfo, catego.medium, catego.compound)
                    print(len(hdfcompoundonecontrol.fields), len(hdfcompoundonecontrol.data),
                          len(hdfcompoundonecontrol.std), len(hdfcompoundonecontrol.std_inh),
                          len(hdfcompoundonecontrol.normalized_perc), len(hdfcompoundonecontrol.normalized))
                    singlecomponecontrol = SingleCompoundonecontrol(header_info.branch, hdfcompoundonecontrol.celline,
                                                                    hdfcompoundonecontrol.seeding,
                                                                    hdfcompoundonecontrol.fields,
                                                                    file_info.elapsed, hdfcompoundonecontrol.control,
                                                                    hdfcompoundonecontrol.compoundalone,
                                                                    hdfcompoundonecontrol.condition,
                                                                    catego.stdinfo, hdfcompoundonecontrol.fieldsmedium,
                                                                    hdfcompoundonecontrol.datamedium,
                                                                    hdfcompoundonecontrol.stdmedium,
                                                                    hdfcompoundonecontrol.fieldsmediuminhibited,
                                                                    hdfcompoundonecontrol.inhibitedmedium,
                                                                    hdfcompoundonecontrol.normalizedtranslationmedium,
                                                                    hdfcompoundonecontrol.std_inhmedium,
                                                                    information_readout, readout_units)

                    singlecomponecontrol.get_ic_txt_path(file_names.fileictxtresultspath)
                    singlecomponecontrol.open_pdf(file_names.filepdfresultspath, i)
                    if len(hdfcompoundonecontrol.elapse) > 1:
                        singlecomponecontrol.confluencyonecontrolovertime(hdfcompoundonecontrol.data,
                                                                          hdfcompoundonecontrol.std)
                        if expected_effect == 0:
                            singlecomponecontrol.inhibitiononecontrolovertime(hdfcompoundonecontrol.inhibited,
                                                                              hdfcompoundonecontrol.std_inh, 2, 1)
                        elif expected_effect == 1:
                            singlecomponecontrol.inhibitiononecontrolovertime(
                                hdfcompoundonecontrol.normalizedtranslation, hdfcompoundonecontrol.std_inh, 1, 1)

                    singlecomponecontrol.doseresponseonecontrol(hdfcompoundonecontrol.normalized_perc,
                                                                hdfcompoundonecontrol.std)

                    singlecomponecontrol.close_pdf()

            if catego.experimentype == "drug_combination":
                print("drug_combination")
                # save information on HDF5 file
                hdfcompoundcomb = Compoundcombination(catego.experimentype, catego.branch,
                                                      file_names.filehdf5resultspath, header_info.newheader,
                                                      file_info.elapsed, file_info.date_info, file_info.date,
                                                      file_info.data, file_info.std, catego.stdinfo,
                                                      catego.medium, catego.compound)

                print(len(hdfcompoundcomb.fields), len(hdfcompoundcomb.data), len(hdfcompoundcomb.std),
                      len(hdfcompoundcomb.std_inh), len(hdfcompoundcomb.normalized_perc), len(hdfcompoundcomb.normalized))

                comb = ExperimentCombination(hdfcompoundcomb.fields, file_info.elapsed, catego.branch,
                                             hdfcompoundcomb.control, hdfcompoundcomb.compoundalone,
                                             hdfcompoundcomb.possiblecombination,
                                             hdfcompoundcomb.possiblecombinationsize, hdfcompoundcomb.inhibited,
                                             hdfcompoundcomb.std_inh, hdfcompoundcomb.celline, hdfcompoundcomb.seeding,
                                             hdfcompoundcomb.condition, hdfcompoundcomb.std_info,
                                             hdfcompoundcomb.medium, information_readout, readout_units)
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
                    if expected_effect == 0:
                        comb.inhibition(hdfcompoundcomb.inhibited, hdfcompoundcomb.std_inh, 2, 1)
                    elif expected_effect == 1:
                        comb.inhibition(hdfcompoundcomb.normalizedtranslation, hdfcompoundcomb.std_inh, 1, 1)

                else:
                    print("1 time point for combo")
                    if expected_effect == 0:
                        comb.endpointinhibition(hdfcompoundcomb.inhibited, hdfcompoundcomb.std_inh, 2, 1)
                    elif expected_effect == 1:
                        comb.endpointinhibition(hdfcompoundcomb.normalizedtranslation, hdfcompoundcomb.std_inh, 1, 1)

                comb.doseresponse(hdfcompoundcomb.normalized_perc, hdfcompoundcomb.std)
                comb.close_pdf()

            if catego.experimentype == "Genetic_perturbagen":
                print("Genetic_perturbagen")
                # save information on HDF5 file
                hdfgenetic = Geneticperturbagem(catego.experimentype, catego.branch, file_names.filehdf5resultspath,
                                                header_info.newheader, file_info.elapsed, file_info.date_info,
                                                file_info.date, file_info.data, file_info.std, catego.stdinfo,
                                                catego.medium, catego.compound)

                print(len(hdfgenetic.fields), len(hdfgenetic.data), len(hdfgenetic.std),
                      len(hdfgenetic.std_inh), len(hdfgenetic.normalized_perc), len(hdfgenetic.normalized))
                print(hdfgenetic.compoundalone)
                geneticpert = ExperimentGeneticPerturbagem(header_info.branch, hdfgenetic.celline, hdfgenetic.seeding,
                                                           hdfgenetic.fields, file_info.elapsed, hdfgenetic.control,
                                                           hdfgenetic.compoundalone, hdfgenetic.condition,
                                                           catego.stdinfo, hdfgenetic.fieldsmedium,
                                                           hdfgenetic.datamedium, hdfgenetic.stdmedium,
                                                           hdfgenetic.fieldsmediuminhibited,
                                                           hdfgenetic.inhibitedmedium, hdfgenetic.std_inhmedium,
                                                           information_readout, readout_units)
                Savetxt(file_names.fileioriginaldatapath, file_info.date_info, hdfgenetic.fields,
                        hdfgenetic.data, file_info.elapsed)

                geneticpert.open_pdf(file_names.filepdfresultspath, i)
                if len(hdfgenetic.elapse) > 1:
                    geneticpert.perturbagemovertime(hdfgenetic.data, hdfgenetic.std, 0, 0, 0)
                    if expected_effect == 0:
                        geneticpert.perturbagemovertime(hdfgenetic.inhibited, hdfgenetic.std_inh, 2, 1, 1)
                    elif expected_effect == 1:
                            geneticpert.perturbagemovertime(hdfgenetic.normalizedtranslation, hdfgenetic.std_inh, 1, 1, 1)

                else:
                    print("1 time point for genetic perturbagem")
                    if expected_effect == 0:
                        geneticpert.perturbagemendpoint(hdfgenetic.inhibited, hdfgenetic.std_inh, 2,
                                                        1, 1)
                    elif expected_effect == 1:
                        geneticpert.perturbagemendpoint(hdfgenetic.normalizedtranslation, hdfgenetic.std_inh, 1, 1, 1)

                geneticpert.close_pdf()

            if catego.experimentype == "genetic-chemical_perturbagen":
                print("genetic-chemical_perturbagen")
                # save information on HDF5 file
                hdfgeneticchemical = Geneticchemicalperturbagem(catego.experimentype, catego.branch,
                                                                file_names.filehdf5resultspath, header_info.newheader,
                                                                file_info.elapsed, file_info.date_info, file_info.date,
                                                                file_info.data, file_info.std, catego.stdinfo,
                                                                catego.medium, catego.compound)

                print(len(hdfgeneticchemical.fields), len(hdfgeneticchemical.data), len(hdfgeneticchemical.std),
                      len(hdfgeneticchemical.std_inh), len(hdfgeneticchemical.normalized_perc),
                      len(hdfgeneticchemical.normalized))
                print("medium information", len(hdfgeneticchemical.fieldsmedium), len(hdfgeneticchemical.datamedium),
                      len(hdfgeneticchemical.stdmedium), len(hdfgeneticchemical.std_inhmedium),
                      len(hdfgeneticchemical.normalized_percmedium), len(hdfgeneticchemical.normalizedmedium))
                #
                geneticchemical = GeneticChemicalPerturbagem(hdfgeneticchemical.fields, file_info.elapsed,
                                                             hdfgeneticchemical.branch, hdfgeneticchemical.control,
                                                             hdfgeneticchemical.compoundalone,
                                                             hdfgeneticchemical.possiblecombination,
                                                             hdfgeneticchemical.possiblecombinationsize,
                                                             hdfgeneticchemical.inhibited,
                                                             hdfgeneticchemical.std_inh, hdfgeneticchemical.celline,
                                                             hdfgeneticchemical.seeding, hdfgeneticchemical.condition,
                                                             hdfgeneticchemical.std_info,
                                                             hdfgeneticchemical.fieldsmedium,
                                                             hdfgeneticchemical.fieldsmediuminhibited,
                                                             hdfgeneticchemical.datamedium,
                                                             hdfgeneticchemical.normalizedtranslationmedium,
                                                             hdfgeneticchemical.inhibitedmedium,
                                                             hdfgeneticchemical.stdmedium,
                                                             hdfgeneticchemical.std_inhmedium,
                                                             information_readout, readout_units)

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
                    if expected_effect == 0:
                        geneticchemical.confluencyovertime(hdfgeneticchemical.inhibited, hdfgeneticchemical.std_inh,
                                                           2, 1)
                    elif expected_effect == 1:
                        geneticchemical.confluencyovertime(hdfgeneticchemical.normalizedtranslation,
                                                           hdfgeneticchemical.std_inh, 1, 1)

                else:
                    print("1 time point for genetic-chemical perturbagem")
                    if expected_effect == 0:
                        geneticchemical.endpointinhibition(hdfgeneticchemical.inhibited, hdfgeneticchemical.std_inh,
                                                           2, 1)
                    elif expected_effect == 1:
                        geneticchemical.endpointinhibition(hdfgeneticchemical.normalizedtranslation,
                                                           hdfgeneticchemical.std_inh, 1, 1)

                geneticchemical.doseresponse(hdfgeneticchemical.normalized_perc, hdfgeneticchemical.std_inh)

                geneticchemical.close_pdf()

    with p.oneshot():
        print(i)
        report = open("informaiton_CPU_time.txt", "a")
        report.write("\n")
        report.write("file name: " + i + '\n')
        report.write(datetime.datetime.fromtimestamp(p.create_time()).strftime("%Y-%m-%d %H:%M:%S") + '\n')
        report.write("cpu_times, user= " + '\t' + str(p.cpu_times()[0]) + '\n')
        report.write("cpu_times, system = " + '\t'+str(p.cpu_times()[1]) + '\n')
        report.write('cpu_percent= ' + '\t'+ str(p.cpu_percent()) + '\n')
        report.write("ppid= " + '\t'+ str(p.ppid()) + '\n')
        report.write("status =" + '\t'+ str(p.status()) + '\n')
        report.write("memory_info, rss= " +  '\t'+ str(p.memory_info()[0]) + '\n')
        report.write("memory_info, vms = " +  '\t'+ str(p.memory_info()[1]) + '\n')
        report.write("io_counters, read_count = " + '\t'+ str(p.io_counters()[0]) + '\n')
        report.write("io_counters, write_count = " + '\t'+ str(p.io_counters()[1]) + '\n')
        report.write("io_counters, read_bytes = " + '\t' + str(p.io_counters()[2]) + '\n')
        report.write("io_counters, write_bytes = " + '\t'+ str(p.io_counters()[3]) + '\n')
        report.close()
        # print("name= ", p.name()) # execute internal routine once collecting multiple info
        print("cpu times= ", p.cpu_times())  # return cached value seconds
        print('cpu_percent= ', p.cpu_percent())  # return cached value
        # print("creat_time= ", p.create_time())  # return cached value
        print(datetime.datetime.fromtimestamp(p.create_time()).strftime("%Y-%m-%d %H:%M:%S"))
        print("ppid= ", p.ppid())  # return cached value
        print("status =", p.status())  # return cached value
        # print("threads= ", p.threads())
        print("memory_info = ", p.memory_info())
        # print("memory maps =", p.memory_maps())
        print("io_counters = ", p.io_counters())
        print("check====", psutil.cpu_freq())
    # Number of logical CPUs in the system
    print("end=psutil.cpu_count() = {0}".format(psutil.cpu_count()))
