# HTSplotter
An end-to-end data processing, analysis and visualisation tool for chemical and genetic in vitro perturbation screens

website [manual](HTSplotter/web/images/HTSplotterManual.pdf).

## Using as python library:

Install

```
pip install HTSplotter==0.10

```

Uninstall

```
pip uninstall HTSplotter==0.10

```

#### HTSplotter tutorial- running from python

##### Import HTSplotter library

```
import HTSplotter as HTSP

```


##### Running Analysis

Both start with initializing HTSplotter object

```
hts = HTSP.Analyser()
```

###### Setting up inputs from shell

Setting up inputs

```
# set main folder path
hts.main_folder = 'path/to/analyses_folder/'

# set input files path
hts.input_path = 'path/to/inputfiles/'

# set path for extracted information from the headers
hts.information_extracted = 'path/to/information_extracted_files/'

# set path for extracted information from the headers
hts.results_path = 'path/to/outputfiles/'

# set 0 if it is not a biological replicate analysis or 1 if it is
hts.biological_replicate = 0

# set 0 if user input is aimed, set 1 if not
hts.userinput = 0

# set readout information
hts.information_readout = 'confluency'

# set readout units
hts.readout_units = '(%)'

# expected effect, if inhibition, 0, if enhancement, 1
hts.expected_effect = 0

# in case of biological replicate analysis give a file name. If it is not a biological replicate analysis set 0
hts.file_name_br = 0

# list of file names
hts.files_list = ['drugscreen_1timepoint.txt', 'drug_combination_screen_1timepoint.txt']

```

Executing analysis

```
# execute analysis
htsplotter.execute()

```

###### Setting up inputs from input file

Executing analysis directly from input file

```
hts.execute_from_file('path/to/input/file/filename.txt')
```

#### prepare input_file
HTSplotter allows an input file, as txt, with all atributes and files. In this way, all files to be analysed with the same atributes can be grouped in one block
check [input file](input_file.txt) example
```
import HTSplotter as HTSP

# input file with different atributes and different files
hts.execute_from_file('PATH_from_your_input_file/input_file.txt')

```

