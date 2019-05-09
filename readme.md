# deamiDATE 1.0

deamiDATE 1.0 is a site-specific deamidation tool for palaeproteomics.

Given MaxQuant result output (in particular, the `evidence.txt` and `peptides.txt` files), it produces two graphs, one showing bulk deamidation, one showing site-specific deamidation.
Deamidation is calculated per experiment, and per protein in experiment. If any filtering of experiments or proteins is required, it is recommended this is done before running deamiDATE.

* Note: only tested on UNIX * 

## Use:

### Normal Use:
1. Download files from GitHub
2. Unzip
3. Navigate to deamiDATE directory
4. Navigate to Packaged Version directory
5. Make `run` executable, e.g. `chmod +x run`
6. Run program e.g. `./run` with the following args:
	1. A directory containing at least `evidence.txt` and `peptides.txt`, e.g. `./run ../Test Data`
	2. A text file containing a list (one per line) of relevant proteins to "filter in", e.g. `./run ../Test Data ../protein_list.txt`

Results will appear in the directory of the MQ output files


### Advanced Use

*Note: only use if you want to mess with stuff*

This is just running the python script without the package details, will require a tonne of dependencies, you can work out which. 
Do this if you'd like to mess with default options etc.

Steps 1 - 3 as per Unix instructions, then

1. Make `run.py` executable, e.g. `chmod +x run`
2. Run program e.g. `./run` with the following args:
	1. A directory containing at least `evidence.txt` and `peptides.txt`, e.g. `./run ../Test Data`
	2. A text file containing a list (one per line) of relevant proteins to "filter in", e.g. `./run ../Test Data ../protein_list.txt`
 
Results will appear in the directory of the MQ output files


## Test Data

Test data is output of MaxQuant run on tapir files from:
Welker, F. et al. Ancient proteins resolve the evolutionary history of Darwin’s South American ungulates. Nature 522, 81–84 (2015).
