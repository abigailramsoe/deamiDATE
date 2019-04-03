# deamiDATE 1.0

deamiDATE 1.0 is a site-specific deamidation tool for palaeproteomics.

Given MaxQuant result output (in particular, the evidence.txt and peptides.txt files), it produces two graphs, one showing bulk deamidation, one showing site-specific deamidation.
Deamidation is calculated per experiment, and per protein in experiment. If any filtering of experiments or proteins is required, it is recommended this is done before running deamiDATE.

## Regular use:

### Unix:
1. Download files from GitHub
2. Unzip
3. Navigate to deamiDATE 1.0 directory in the "Packaged Version" folder
4. Run program, with include an argument for the folder with MQ result files,
e.g. `./run ../Test Data/Tapir` must be a directory containing the evidence.txt and peptides.txt

Results will appear in the directory given in step 5.

## Advanced use:

You can also run it as a python script, there are some dependencies you will need to install, namely:
* matplotlib
* numpy
* pandas
* csv

This allows you to edit options that are hard-coded in, for example, change `show` to `True` to have plots pop up, instead of just save.

Run as above, with an MQ output folder as an argument, e.g. `./run.py ../Test Data/ Tapir`

## Test data

One folder of MaxQuant results is supplied as a test folder. It is a tapir sample from

Welker F, Collins MJ, Thomas JA, Wadsley M, Brace S, Cappellini E, Turvey ST, Reguero M, Gelfo JN, Kramarz A, Burger J, Thomas-Oates J, Ashford DA, Ashton PD, Rowsell K, Porter DM, Kessler B, Fischer R, Baessmann C, Kaspar S, Olsen JV, Kiley P, Elliott JA, Kelstrup CD, Mullin V, Hofreiter M, Willerslev E, Hublin JJ, Orlando L, Barnes I, MacPhee RD, Ancient proteins resolve the evolutionary history of Darwin's South American ungulates. Nature, 522(7554):81-4(2015) 

Raw files available [here](http://proteomecentral.proteomexchange.org/cgi/GetDataset?ID=PXD001411)

