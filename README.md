# CLBTope: Prediction of linear and conformational B-cell epitopes.
## Introduction
CLBTope is developed to predict, scan, and, design the both types (linear/conformational) of B-cell epitopes using sequence information only. In the standalone version, Random Forest based model is implemented along with the BLAST search, named it as hybrid approach.
CLBTope is also available as web-server at https://webs.iiitd.edu.in/raghava/clbtope. Please read/cite the content about the clbtope for complete information including algorithm behind the approach.
____________________________________________________________________
## Model.zip (optional)
Compressed file containing the CLBTope model. 
The same can be downloaded from the standalone page of the website. 
Extract the contents of the compressed file to the same location.
After extraction, you should see the extracted folder in the standalone folder alongside the "Model.zip" file.

Note: The standalone version will automatically fetch the model, eliminating the need for repeated downloads. However, a stable internet connection is required for this process.
____________________________________________________________________

## PIP Installation
PIP version is also available for easy installation and usage of this tool. The following command is required to install the package 
```
pip install clbtope
```
To know about the available option for the pip package, type the following command:
```
clbtope -h
```

## Standalone
The Standalone version of transfacpred is written in python3 and following libraries are necessary for the successful run:
- scikit-learn
- Pandas
- Numpy
- blastp

## Minimum USAGE
To know about the available option for the stanadlone, type the following command:
```
python clbtope.py -h
```
To run the example, type the following command:
```
python3 clbtope.py -i example_input.fa
```
This will predict if the submitted sequences can B-cell epitope or not. It will use other parameters by default. It will save the output in "outfile.csv" in CSV (comma seperated variables).

## Full Usage
```
usage: python3 clbtope.py [-h]
                       [-i INPUT
                       [-o OUTPUT]
                       [-j {1,2,3}]
                       [-t THRESHOLD]
                       [-w {8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30}]
                       [-d {1,2}]
```
```
Please provide following arguments for successful run

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Input: protein or peptide sequence(s) in FASTA format
                        or single sequence per line in single letter code
  -o OUTPUT, --output OUTPUT
                        Output: File for saving results by default outfile.csv
  -j {1,2,3}, --job {1,2,3}
                        Job Type: 1:Predict, 2: Design, 3:Scan, by default 1
  -t THRESHOLD, --threshold THRESHOLD
                        Threshold: Value between 0 to 1 by default 0.16
  -w {8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30}, --winleng {8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30}
                        Window Length: 8 to 30 (scan mode only), by default 9
  -d {1,2}, --display {1,2}
                        Display: 1:B-cell epitope only, 2: All peptides, by default 1
```

**Input File:** It allow users to provide input in the FASTA format.

**Output File:** Program will save the results in the CSV format, in case user do not provide output file name, it will be stored in "outfile.csv".

**Threshold:** User should provide threshold between 0 and 1, by default its 0.53.

**Job:** User is allowed to choose between three different modules, such as, 1 for prediction, 2 for Designing and 3 for scanning, by default its 1.

**Window length**: User can choose any pattern length between 8 and 30 in long sequences. This option is available for only scanning module.

**Display type:** This option allow users to fetch either only B-cell epitope by choosing option 1 or prediction against all peptides by choosing option 2.

CLBTope Package Files
=======================
It contains the following files, brief descript of these files given below

INSTALLATION                    : Installations instructions

LICENSE                         : License information

README.md                       : This file provide information about this package

model.zip                       : This zipped file contains the compressed version of model

envfile                         : This file compeises of paths for the database and blastp executable

clbtope.py                      : Main python program

example_input.fa                : Example file contain peptide sequenaces in FASTA format

example_predict_output.csv      : Example output file for predict module

example_scan_output.csv         : Example output file for scan module

example_design_output.csv       : Example output file for design module
