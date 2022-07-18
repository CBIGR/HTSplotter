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
    # main_folder = "HTSplotter_master_ugent/" \
    #               "experiment_type/drug_combination/"
    main_folder = "C:/Users/cdcarval/Dropbox (speleman lab)/Personal Lab/HTSplotter_master_ugent/" \
                  "experiment_type/"
    # main_folder = "C:/Users/cdcarval/Dropbox (speleman lab)/Personal Lab/Carol personal Lab/wet lab/siRNA/siTPX2/" \
    #               "Incucyte-rawdata/"
    # main_folder = "C:/Users/cdcarval/Dropbox (speleman lab)/Personal Lab/Carol personal Lab/wet lab/Drugging/Drug/" \
    #               "CR-31-B/CR-31-B_IC50/HTSplotter_data/"
    input_path = main_folder + 'drug_combination/'   # + "Incucyte_RawData"
    information_extracted = main_folder #+ "HTSplotter_results/" #+ 'HTSplotter_results/' #+ 'resultsHTSplotter/' # + "Information_extracted_files/"
    results_path = main_folder #+ "HTSplotter_results/"#+ 'HTSplotter_results/'  #+ 'resultsHTSplotter/'  # + "Output_results/"

    biological_replicate = 0  # 0 if not, 1 if it is
    userinput = 0  # 1 if yes
    information_readout = "confluency"  # default = confluency = 0; add effect name = 1
    readout_units = "(%)"
    expected_effect = 0  # 0 = "inhibition"; 1 = enhanced

    # files_list = ['drug_combination_several_time_pointsZIP']
    # #### files to test drug screen analysis
    # files_list = ['drugscreen_1timepoint']
    # files_list = ['drugscreen_severaltimepoint_severalcontrol']
    # files_list = ['Rep3 (1)'] #, 'Rep2 (1)', 'Rep3 (1)', 'Rep1 (1)/
    # files_list = ['drugscreen_severaltimepoint_1control']
    # files_list = ['20191125_siTPX2_CLBGA']
    # files_list = ['20191125_siTPX2_IMR32']

    # files_list = ['20191125_siTPX2_CLBGA_after24hseeding', '20191125_siTPX2_IMR32_after24hseeding']
    # files_list = ['20211203_CLBGA CR31B dosage range_BR2', '20211203_CLBGA CR31B dosage range_BR3',
    #               '20211203_CLBGA CR31B dosage range_BR4', '20211203_SKNBE2c CR31B dosage range_BR2',
    #               '20211203_SKNBE2c CR31B dosage range_BR3', '20211203_SKNBE2c CR31B dosage range_BR4']
    # files_list = ['20211203_CLBGA CR31B dosage range_BR2', '20211203_CLBGA CR31B dosage range_BR3',
    #               '20211203_CLBGA CR31B dosage range_BR4']
    files_list = ['20220707_SRA737_CLBGA_IC', '20220707_SRA737_IMR32_IC']
    # files_list = ['20211203_SKNBE2c CR31B dosage range_BR2',
    #               '20211203_SKNBE2c CR31B dosage range_BR3', '20211203_SKNBE2c CR31B dosage range_BR4']
    # files_list = ['gene_perturbagen_severaltimepoints']
    # files_list = [ '20191125_siTPX2_IMR32', '20200228_CLBGA_siTPX2',
                  # '20200228_IMR32_siTPX2', '20200228_SHSY5Y_siTPX2', '20200228_SKNBE2c_siTPX2',
    #               '20200306_CLBGA_siTPX2', '20200306_IMR32_siTPX2', '20200306_SHSY5Y_siTPX2',
    #               '20200306_SKNBE2c_siTPX2']
    # ### files to teste genetic-chemical perturbagem
    # files_list = ['202202_IMR32_siNTC_vs_si47']
    # files_list = ['genetic-chemical_perturbagen_several-time_points']
    # SOX11 data
    # files_list = ['SHEP_SOX11_OE_G2_Celst14-02-2019']
    # file_name_br = "BR_"
    synergy_method = 0  # 0-->Bliss; 1-->HSA; 2--> ZIP
    file_name_br = 'CLBGA CR31B dosage range_BR'
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

            catego = Categorisation(header_info)

            out = Outputfile(eachfile, biological_replicate, file_names, catego, header_info, file_info.elapsed)

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
            brcombdata = BRcombinationstructure(count, file_br_names.filehdf5resultspath, file_name_br,
                                                header_info, catego, file_info)
            print(len(brcombdata.fields), len(brcombdata.data), len(brcombdata.std), len(brcombdata.std_inh),
                  len(brcombdata.normalized_perc), len(brcombdata.normalized))
            for g in range(len(brcombdata.possiblecombination)):
                print(brcombdata.possiblecombinationsize[g], brcombdata.possiblecombination[g])

            comb = ExperimentCombination(synergy_method, brcombdata, file_info, catego.branch,
                                         information_readout, readout_units, file_br_names, biological_replicate,
                                         file_name_br)
            # ###
            # add predicted and bliss score to the HDF5 file
            brcombdata.add_predictedbiscore(comb.comb_name_per_group, comb.effect_per_group,
                                            comb.synergy_score_per_group, synergy_method)
            # over time data:

            if len(brcombdata.elapse) > 1:
                comb.inhibition(brcombdata.data, brcombdata.std)
                if expected_effect == 0:
                    comb.inhibition(brcombdata.inhibited, brcombdata.std_inh, 2, 1)
                elif expected_effect == 1:
                    comb.inhibition(brcombdata.normalizedtranslation, brcombdata.std_inh, 1, 1)
                comb.doseresponse()
            else:
                print("1 time point for combination biological replicates")
                if expected_effect == 0:
                    comb.endpointinhibition(brcombdata.inhibited, brcombdata.std_inh, 2, 1)
                if expected_effect == 1:
                    comb.endpointinhibition(brcombdata.normalizedtranslation, brcombdata.std_inh, 1, 1)
                comb.doseresponse(1)
            comb.close_pdf()

        if catego.experimentype == "drug_screen":
            print("Compound_screen", catego.experimentype)
            if len(catego.compound) == len(catego.control):
                # 1 control for each compound
                print("Compound_screen several control BR", len(catego.compound), len(catego.control), catego.medium)
                brcomscreendata = BRcompoundscreenseveralcontrol(count, file_br_names.filehdf5resultspath,
                                                                 file_name_br, header_info, catego, file_info)

                print(len(brcomscreendata.fields), len(brcomscreendata.data), len(brcomscreendata.std),
                      len(brcomscreendata.std_inh), len(brcomscreendata.normalized_perc),
                      len(brcomscreendata.normalized))
                for g in range(len(brcomscreendata.compoundalone)):
                    print(brcomscreendata.compoundalone[g])
                print(brcomscreendata.medium)
                singlecompond = SingleCompound(header_info.branch, file_info,
                                               information_readout, readout_units,
                                               biological_replicate, brcomscreendata,
                                               file_br_names, file_name_br)

                if len(brcomscreendata.elapse) > 1:
                    singlecompond.inhibitionovertime(brcomscreendata.data, brcomscreendata.std)
                    if expected_effect == 0:
                        singlecompond.inhibitionovertime(brcomscreendata.inhibited, brcomscreendata.std_inh, 2, 1)
                    elif expected_effect == 1:
                        singlecompond.inhibitionovertime(brcomscreendata.normalizedtranslation, brcomscreendata.std_inh,
                                                         1, 1)
                    singlecompond.doseresponse()
                else:
                    singlecompond.doseresponse(1)

                singlecompond.close_pdf()
            else:
                print("drug_screen 1 control BR", len(catego.compound), len(catego.control), catego.medium)
                brcomscreenonecontrol = BRcompoundscreenonecontrol(count, file_br_names.filehdf5resultspath,
                                                                   file_name_br, header_info, catego, file_info)
                print("one control")
                print(len(brcomscreenonecontrol.fields), len(brcomscreenonecontrol.data),
                      len(brcomscreenonecontrol.std), len(brcomscreenonecontrol.std_inh),
                      len(brcomscreenonecontrol.normalized_perc), len(brcomscreenonecontrol.normalized))
                for g in range(len(brcomscreenonecontrol.compoundalone)):
                    print(brcomscreenonecontrol.compoundalone[g])
                print(brcomscreenonecontrol.medium)

                singlecomponecontrol = SingleCompoundonecontrol(header_info.branch, file_info,
                                                                information_readout, readout_units,
                                                                biological_replicate, brcomscreenonecontrol,
                                                                file_br_names, file_name_br)

                if len(brcomscreenonecontrol.elapse) > 1:
                    singlecomponecontrol.inhibitiononecontrolovertime(brcomscreenonecontrol.data,
                                                                      brcomscreenonecontrol.std)
                    if expected_effect == 0:
                        singlecomponecontrol.inhibitiononecontrolovertime(brcomscreenonecontrol.inhibited,
                                                                          brcomscreenonecontrol.std_inh, 2, 1)
                    elif expected_effect == 1:
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
                                                 file_name_br, header_info, catego, file_info)
            print("did it?>>>")
            print(len(BRhdf5genetic.fields), len(BRhdf5genetic.data), len(BRhdf5genetic.std),
                  len(BRhdf5genetic.std_inh), len(BRhdf5genetic.normalized_perc), len(BRhdf5genetic.normalized))
            print(BRhdf5genetic.compoundalone)

            geneticpert = ExperimentGeneticPerturbagem(header_info, BRhdf5genetic, file_info, catego, file_br_names,
                                                       information_readout, readout_units, file_name_br,
                                                       biological_replicate)

            if len(BRhdf5genetic.elapse) > 1:
                geneticpert.perturbagemovertime(BRhdf5genetic.data, BRhdf5genetic.std, 0, 0, 0)
                if expected_effect == 0:
                    geneticpert.perturbagemovertime(BRhdf5genetic.inhibited, BRhdf5genetic.std_inh, 2, 1, 1)
                elif expected_effect == 1:
                    geneticpert.perturbagemovertime(BRhdf5genetic.normalizedtranslation, BRhdf5genetic.std_inh, 1, 1, 1)

            else:
                print("1 time point for genetic-chemical perturbagem")
                if expected_effect == 0:
                    geneticpert.perturbagemendpoint(BRhdf5genetic.inhibited, BRhdf5genetic.std_inh, 2, 1)
                elif expected_effect == 1:
                    geneticpert.perturbagemendpoint(BRhdf5genetic.normalizedtranslation, BRhdf5genetic.std_inh, 1, 1)

            geneticpert.close_pdf()

        if catego.experimentype == "genetic-chemical_perturbagen":
            print("genetic-chemical_perturbagen")
            # save information on HDF5 file
            brhdfgeneticchemical = BRgeneticchemicalperturbagem(count, file_br_names.filehdf5resultspath,
                                                                file_name_br, header_info, catego, file_info)

            print(len(brhdfgeneticchemical.fields), len(brhdfgeneticchemical.data), len(brhdfgeneticchemical.std),
                  len(brhdfgeneticchemical.std_inh), len(brhdfgeneticchemical.normalized_perc),
                  len(brhdfgeneticchemical.normalized))
            print("medium information", len(brhdfgeneticchemical.fieldsmedium), len(brhdfgeneticchemical.datamedium),
                  len(brhdfgeneticchemical.stdmedium), len(brhdfgeneticchemical.std_inhmedium),
                  len(brhdfgeneticchemical.normalized_percmedium), len(brhdfgeneticchemical.normalizedtranslationmedium))
            print(brhdfgeneticchemical.possiblecombination)

            geneticchemical = GeneticChemicalPerturbagem(synergy_method, brhdfgeneticchemical, file_info, file_br_names,
                                                         information_readout, readout_units, file_name_br,
                                                         biological_replicate)

            if len(brhdfgeneticchemical.elapse) > 1:
                geneticchemical.confluencyovertime(brhdfgeneticchemical.data, brhdfgeneticchemical.std, 0, 0)
                if expected_effect == 0:
                    geneticchemical.confluencyovertime(brhdfgeneticchemical.inhibited, brhdfgeneticchemical.std_inh,
                                                       2, 1)
                elif expected_effect == 1:
                    geneticchemical.confluencyovertime(brhdfgeneticchemical.normalizedtranslation,
                                                       brhdfgeneticchemical.std_inh, 1, 1)
                geneticchemical.doseresponse()
            else:
                print("are we here?")
                if expected_effect == 0:
                    geneticchemical.endpointinhibition(brhdfgeneticchemical.inhibited, brhdfgeneticchemical.std_inh,
                                                       2, 1)
                elif expected_effect == 1:
                    geneticchemical.endpointinhibition(brhdfgeneticchemical.normalizedtranslation,
                                                       brhdfgeneticchemical.std_inh, 1, 1)
                geneticchemical.doseresponse(1)

            geneticchemical.close_pdf()

    else:

        for i in files_list:
            # creat different file names and paths:
            file_names = Filenames(i, input_path, results_path, information_extracted)
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

            out = Outputfile(i, biological_replicate, file_names, catego, header_info, file_info.elapsed)

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
                    hdfcompound = Compoundscreen(file_names.filehdf5resultspath, header_info.newheader,
                                                 catego, file_info)

                    print(len(hdfcompound.fields), len(hdfcompound.data), len(hdfcompound.std),
                          len(hdfcompound.std_inh), len(hdfcompound.normalized_perc),
                          len(hdfcompound.normalized))

                    print("medium information", len(hdfcompound.datamedium),
                          len(hdfcompound.stdmedium), len(hdfcompound.std_inhmedium),
                          len(hdfcompound.normalized_percmedium), len(hdfcompound.normalizedmedium))

                    singlecompond = SingleCompound(header_info.branch, file_info,
                                                   information_readout, readout_units, biological_replicate,
                                                   hdfcompound, file_names, i)

                    if len(hdfcompound.elapse) > 1:
                        singlecompond.inhibitionovertime(hdfcompound.data, hdfcompound.std)
                        if expected_effect == 0:
                            singlecompond.inhibitionovertime(hdfcompound.inhibited, hdfcompound.std_inh, 2, 1)
                        elif expected_effect == 1:
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
                                                                    biological_replicate, hdfcompoundonecontrol,
                                                                    file_names, i)

                    if len(hdfcompoundonecontrol.elapse) > 1:
                        singlecomponecontrol.inhibitiononecontrolovertime(hdfcompoundonecontrol.data,
                                                                          hdfcompoundonecontrol.std)
                        if expected_effect == 0:
                            singlecomponecontrol.inhibitiononecontrolovertime(hdfcompoundonecontrol.inhibited,
                                                                              hdfcompoundonecontrol.std_inh, 2, 1)
                        elif expected_effect == 1:
                            singlecomponecontrol.inhibitiononecontrolovertime(
                                hdfcompoundonecontrol.normalizedtranslation, hdfcompoundonecontrol.std_inh, 1, 1)
                        singlecomponecontrol.doseresponseonecontrol()
                    else:
                        singlecomponecontrol.doseresponseonecontrol(1)
                    singlecomponecontrol.close_pdf()

            if catego.experimentype == "drug_combination":
                print("drug_combination")
                print('0-->Bliss; 1-->HSA; 2-->ZIP, the synergy method is: ', synergy_method)
                # save information on HDF5 file
                hdfcompoundcomb = Compoundcombination(file_names.filehdf5resultspath, header_info.newheader,
                                                      catego, file_info)

                print(len(hdfcompoundcomb.fields), len(hdfcompoundcomb.data), len(hdfcompoundcomb.std),
                      len(hdfcompoundcomb.std_inh), len(hdfcompoundcomb.normalized_perc),
                      len(hdfcompoundcomb.normalized))

                comb = ExperimentCombination(synergy_method, hdfcompoundcomb, file_info, catego.branch,
                                             information_readout, readout_units, file_names,
                                             biological_replicate, i)

                # add predicted and bliss score to the HDF5 file
                hdfcompoundcomb.add_predictedbiscore(comb.comb_name_per_group,
                                                     comb.synergy_score_per_group,
                                                     comb.effect_per_group,
                                                     synergy_method)

                # # over time data:
                if len(hdfcompoundcomb.elapse) > 1:
                    comb.inhibition(hdfcompoundcomb.data, hdfcompoundcomb.std)
                    if expected_effect == 0:
                        comb.inhibition(hdfcompoundcomb.inhibited, hdfcompoundcomb.std_inh, 2, 1)
                    elif expected_effect == 1:
                        comb.inhibition(hdfcompoundcomb.normalizedtranslation, hdfcompoundcomb.std_inh, 1, 1)
                    comb.doseresponse()

                else:
                    print("1 time point for combo")
                    if expected_effect == 0:
                        comb.endpointinhibition(hdfcompoundcomb.inhibited, hdfcompoundcomb.std_inh, 2, 1)
                    elif expected_effect == 1:
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
                                                           information_readout, readout_units, i, biological_replicate)

                if len(hdfgenetic.elapse) > 1:
                    geneticpert.perturbagemovertime(hdfgenetic.data, hdfgenetic.std, 0, 0, 0)
                    if expected_effect == 0:
                        geneticpert.perturbagemovertime(hdfgenetic.inhibited, hdfgenetic.std_inh, 2, 1, 1)
                    elif expected_effect == 1:
                        geneticpert.perturbagemovertime(hdfgenetic.normalizedtranslation, hdfgenetic.std_inh, 1, 1, 1)

                else:
                    print("1 time point for genetic perturbagem")
                    if expected_effect == 0:
                        geneticpert.perturbagemendpoint(hdfgenetic.inhibited, hdfgenetic.std_inh, 2, 1)
                    elif expected_effect == 1:
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
                geneticchemical = GeneticChemicalPerturbagem(synergy_method, hdfgeneticchemical, file_info, file_names,
                                                             information_readout, readout_units,
                                                             i, biological_replicate)

                if len(hdfgeneticchemical.elapse) > 1:
                    geneticchemical.confluencyovertime(hdfgeneticchemical.data, hdfgeneticchemical.std, 0, 0)
                    if expected_effect == 0:
                        geneticchemical.confluencyovertime(hdfgeneticchemical.inhibited, hdfgeneticchemical.std_inh,
                                                           2, 1)
                    elif expected_effect == 1:
                        geneticchemical.confluencyovertime(hdfgeneticchemical.normalizedtranslation,
                                                           hdfgeneticchemical.std_inh, 1, 1)
                    geneticchemical.doseresponse()

                else:
                    print("1 time point for genetic-chemical perturbagem")
                    if expected_effect == 0:
                        geneticchemical.endpointinhibition(hdfgeneticchemical.inhibited, hdfgeneticchemical.std_inh,
                                                           2, 1)
                    elif expected_effect == 1:
                        geneticchemical.endpointinhibition(hdfgeneticchemical.normalizedtranslation,
                                                           hdfgeneticchemical.std_inh, 1, 1)
                    geneticchemical.doseresponse(1)

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
