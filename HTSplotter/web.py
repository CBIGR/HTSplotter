#!/usr/bin/python
import os
import time
import psutil
import datetime
import sys
import numpy as np
import glob
import argparse
import shutil
from minio import Minio
from minio.error import S3Error

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

parser = argparse.ArgumentParser(description='Submit HTSplotter job from website')
parser.add_argument('-f', nargs=1, required=True, help='The directory containing file(s) that need to be analysed',
                    metavar='file_dir')
parser.add_argument('-b', nargs=1, required=True, choices=['yes', 'no'], help='Biological replicate analysis?',
                    metavar='bio_rep')
parser.add_argument('-n', nargs=1, required=True, help='Biological replicate desired filename', metavar='biorep_fn')
parser.add_argument('-e', nargs=1, required=True, choices=['inhibition', 'enhancement'], help='Expected effect',
                    metavar='exp_effect')
parser.add_argument('-i', nargs=1, required=True, help='Information readout', metavar='info_readout')
parser.add_argument('-r', nargs=1, required=True, help='Readout unit', metavar='readout_unit')
parser.add_argument('-sy', nargs=1, required=True, help='Synergism/antagonism caculation', metavar='syn_ant')
parser.add_argument('-u', nargs=1, required=True, choices=['yes','no'], help='User input corrected?', metavar='user_input')
parser.add_argument('-t', nargs='*', required=False, default='none', help='The experiment type(s)', metavar='exp_types')
parser.add_argument('-sh', nargs=1, required=True, help='S3 host', metavar='s3_host')
parser.add_argument('-sb', nargs=1, required=True, help='S3 bucket', metavar='s3_bucket')
parser.add_argument('-sa', nargs=1, required=True, help='S3 access', metavar='s3_access')
parser.add_argument('-ss', nargs=1, required=True, help='S3 key', metavar='s3_key')
# parser.add_argument('-m', nargs=1, required=True, help='The submitters email address (used for reporting)', metavar="user_email")

# # Example usage
# python3 /var/www/html/scripts/web.py -f /var/www/html/user_files/User_KSkyuLno_jasper_anckaert/ -b 'no' -n 'IMR32_BR' -e 'inhibition' -i 'confluency' -r '%' -u 'no' -sh 'play.min.io' -sb 'htsplotter' -sa 'Q3AM3UQ867SPQQA43P2F' -ss 'zuf+tfteSlswRu7BJ86wekitnifILbZam1KYY3TG'
# python3 /var/www/html/scripts/web.py -f /var/www/html/user_files/User_KSkyuLno_jasper_anckaert/ -b 'no' -n 'IMR32_BR' -e 'inhibition' -i 'confluency' -r '%' -u 'yes' -t 'drug_combination' -sh 'play.min.io' -sb 'htsplotter' -sa 'Q3AM3UQ867SPQQA43P2F' -ss 'zuf+tfteSlswRu7BJ86wekitnifILbZam1KYY3TG'
# python3 /var/www/html/scripts/web.py -f /var/www/html/user_files/User_3P5SX5YB_jasper_anckaert/ -b 'yes' -n 'IMR32_BR' -e 'inhibition' -i 'confluency' -r '%' -u 'no' -sh 'play.min.io' -sb 'htsplotter' -sa 'Q3AM3UQ867SPQQA43P2F' -ss 'zuf+tfteSlswRu7BJ86wekitnifILbZam1KYY3TG'
# python3 /var/www/html/scripts/web.py -f /var/www/html/user_files/User_3P5SX5YB_jasper_anckaert/ -b 'yes' -n 'IMR32_BR' -e 'inhibition' -i 'confluency' -r '%' -u 'yes' -t 'drug_combination' 'drug_combination' -sh 'play.min.io' -sb 'htsplotter' -sa 'Q3AM3UQ867SPQQA43P2F' -ss 'zuf+tfteSlswRu7BJ86wekitnifILbZam1KYY3TG'


if __name__ == '__main__':
	start = time.time()
	# Number of logical CPUs in the system
	p = psutil.Process()

	# Parse arguments
	args = parser.parse_args()
	# print(args)

	input_path = args.f[0] + "Inputfile/"
	results_path = args.f[0] + "Output_html/"
	# information_extracted = input_path + "Information_extracted_files/"
	information_extracted = args.f[0]

	biological_replicate = 1 if args.b[0] == 'yes' else 0 # 0 if not, 1 if it is
	# experiment_type_arrayed = 0  # 0 = no, 1= yes
	# each_file_each_plate = 1  # 0 no, 1 = yes
	userinput = 1 if args.u[0] == 'yes' else 0 # 1 if yes
	information_readout = args.i[0] #"confluency"  # default = confluency = 0; add effect name = 1
	readout_units = args.r[0] #"(%)"
	expected_effect = 0  if args.e[0] == 'inhibition' else 1 # 0 = "inhibition"; 1 = enhanced
	file_name_br = args.n[0]
	files_list = [f for f in os.listdir(input_path) if f.endswith('.txt')]
	exp_types = args.t
	synergy_method = 0 if args.sy[0] == 'Bliss' else 1 if args.sy[0] == 'HSA' else 2 if args.sy[0] == 'ZIP' else 2# 0-->Bliss; 1-->HSA; 2--> ZIP

	os.makedirs(results_path, exist_ok=True)

	if biological_replicate == 1:
		print("Processing Bioloical replicate")
		count = 0
		file_br_names = Filenames(file_name_br, input_path, results_path, information_extracted)

		if userinput == 0:

		# first go over each file before processing it
			for eachfile in files_list:
				file_names = Filenames(eachfile, input_path, results_path, information_extracted)

				# get potential file names
				file_hdf5_name = os.path.basename(file_names.filehdf5resultspath)
				file_pdf_name = os.path.basename(file_names.filepdfresultspath)
				file_txt_name = os.path.basename(file_names.fileitxtnputpath)
				file_png_name = os.path.basename(file_names.filepngresultspath)
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

				# set error file
				print("Experiment type"+str(count)+":",catego.experimentype)
				print("File"+str(count)+":",eachfile)
				# if len(out.error) == 0:
				# 	print("HTSplotter did not identify any error from your input file")
				# 	out.seterrorfile(0)
				if len(out.error) != 0:
					print("check error file, please")
					out.seterrorfile(1)
					sys.exit()
				count += 1
		# check user confirmation
		elif userinput == 1:
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
				

				# print("Output PDF: ",file_name_br)
				if catego.experimentype != exp_types[count]:
					print("=====>check error file", exp_types[count], catego.experimentype)
					break
				print("=====>ok", exp_types[count], catego.experimentype)

				count += 1
				# out.seterrorfile(0)

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

			# zip all files
			print('Zip all output files')
			output_filename = args.f[0]+results_path.replace("/var/www/html/user_files/", "").replace("/Output_html/", "")
			shutil.make_archive(output_filename, 'zip', input_path, 'Output_html')

			# Minio file upload
			upload_fn = output_filename+'.zip'
			upload_alias = os.path.basename(upload_fn)

			client = Minio(
				args.sh[0],#"play.min.io",
				access_key=args.sa[0],#"Q3AM3UQ867SPQQA43P2F",
				secret_key=args.ss[0],#"zuf+tfteSlswRu7BJ86wekitnifILbZam1KYY3TG",
			)
			# found = client.bucket_exists(args.sb[0])
			# if not found:
			# 	client.make_bucket(args.sb[0])
			# else:
			# 	print("Bucket",args.sb[0],"already exists")

			result = client.fput_object(args.sb[0], upload_alias, upload_fn)
			shutil.rmtree(input_path)

	else:
		for i in files_list:
			# create different file names and paths:
			file_names = Filenames(i, input_path, results_path, information_extracted)
			# get potential file names
			file_hdf5_name = os.path.basename(file_names.filehdf5resultspath)
			file_pdf_name = os.path.basename(file_names.filepdfresultspath)
			file_txt_name = os.path.basename(file_names.fileitxtnputpath)
			file_png_name = os.path.basename(file_names.filepngresultspath)

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
			if userinput == 0:
				print("Experiment type: ",catego.experimentype)
				print("File:",i)
				# if len(out.error) == 0:
				# 	print("HTSplotter did not identify any error from your input file")
				# 	out.seterrorfile(0)
				if len(out.error) != 0:
					print("check error file, please")
					out.seterrorfile(1)
					sys.exit()

			# check user confirmation
			elif userinput == 1:
				print("Output file: ",file_pdf_name[:-4]) # remove .pdf extension
				if catego.experimentype != exp_types[0]:
					print("=====>check error file", exp_types[0], catego.experimentype)
					break
				print("=====>ok", exp_types[0], catego.experimentype)
				out.seterrorfile(0)
				
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
						comb.inhibition(hdfcompoundcomb.data, hdfcompoundcomb.std, 1)
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

				# zip all files
				print('Zip all output files')
				output_filename = args.f[0]+results_path.replace("/var/www/html/user_files/", "").replace("/Output_html/", "")
				shutil.make_archive(output_filename, 'zip', results_path)

				# Minio file upload
				upload_fn = output_filename+'.zip'
				upload_alias = os.path.basename(upload_fn)

				client = Minio(
					args.sh[0],#"play.min.io",
					access_key=args.sa[0],#"Q3AM3UQ867SPQQA43P2F",
					secret_key=args.ss[0],#"zuf+tfteSlswRu7BJ86wekitnifILbZam1KYY3TG",
				)
				# found = client.bucket_exists(args.sb[0])
				# if not found:
				# 	client.make_bucket(args.sb[0])
				# else:
				# 	print("Bucket",args.sb[0],"already exists")

				result = client.fput_object(args.sb[0], upload_alias, upload_fn)
				# shutil.rmtree(input_path)
	with p.oneshot():
		print(i)
		report = open(information_extracted+"information_CPU_time.txt", "a")
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
