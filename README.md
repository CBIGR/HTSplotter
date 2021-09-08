# HTSplotter
An end-to-end data processing, analysis and visualisation tool for chemical and genetic in vitro perturbation screens

website [manual](HTSplotter/web/images/HTSplotterManual.pdf).

## Install

```
pip install HTSplotter==0.10
```

## Uninstall

```
pip uninstall HTSplotter==0.10
```

## Import HTSplotter library

```
import HTSplotter as HTSP
```

## Running Analysis

The HTSplotter analysis can be made by setting by introducing the input data from the shell or it can read the inputs from an input file. This allows the integration of HTSplotter into other scripts and the grouping of various analysis into asingle run, even in the case of different inputs.
 In both cases, the Analyser object must be initialized:

```
hts = HTSP.Analyser()
```

### Setting up inputs from shell

To setup the inputs, the attributes of the Analyser object must be changed to the proper ones:

```
# set main folder path
hts.main_folder = 'path/to/main/folder/'

# set readout information
hts.information_readout = 'confluency'

# set readout units
hts.readout_units = '(%)'

# in case of biological replicate analysis give a file name. If it is not a biological replicate analysis set 0
hts.file_name_br = 0

# list of file names
hts.files_list = ['drugscreen_1timepoint.txt', 'drug_combination_screen_1timepoint.txt']

# =0 not biological replicate (default); =1 biological replicate
hts.biological_replicate = 0

# =0 no user validation (default); =1 user validation
hts.userinput = 0

# =0 inhibition (defult); =1 enhancement
hts.expected_effect = 0
```

By default, all input files, log files (information extracted) and output files are set to the main folder path. However, these can be changed to diferent paths to ease data organization: 

```
# set input files path
hts.input_path = 'path/to/input/files/'

# set path for extracted information from the headers
hts.information_extracted = 'path/to/information/extracted/files/'

# set path for extracted information from the headers
hts.results_path = 'path/to/output/files/'
```

To run the analysis, after all inputs have been set, you must use the execute() routine:

```
# execute analysis
htsplotter.execute()
```

### Setting up inputs from input file

The input data can be introduced into an input file, where several analysis blocks can be set with different characteristics. Setting up HTSplotter in this way allows the sequential execution of several execute() methods, for each input data block.
Each analysis block must terminate with #.

```
path/to/to/main/directory/                             # add your folder path e.g. GitHub path
path/to/input/files/directory/                         # input path folder
path/to/information/extracted/files/directory          # information_extracted path folder
path/to/output/files/directory/                        # output_results path folder
0                                                      # biological_replicate 0 if not, 1 if it is
0                                                      # userinput 1 if yes
confluency                                             # information_readout e.g.: "confluency"
(%)                                                    # readout_units e.g. (%)
0                                                      # expected_effect 0 = "inhibition"; 1 = enhanced
0                                                      # file_name_br file name
drug_combination_screen_1timepoint                     # file 1
end                                                    # file 2 (if end, stop reading)
#
path/to/to/main/directory/                             # add your folder path e.g. GitHub path
path/to/input/files/directory/                         # input path folder
path/to/information/extracted/files/directory          # information_extracted path folder
path/to/output/files/directory/                        # output_results path folder
0                                                      # biological_replicate 0 if not, 1 if it is
0                                                      # userinput 1 if yes
confluency                                             # information_readout e.g.: "confluency"
(area)                                                 # readout_units e.g. (%)
0                                                      # expected_effect 0 = "inhibition"; 1 = enhanced
0                                                      # file_name_br file name
drugscreen_1timepoint				       # file 1
end                                                    # file 2 (if end, stop reading)
#
``` 

```
hts.execute_from_file('path/to/input/file/filename.txt')
```

 