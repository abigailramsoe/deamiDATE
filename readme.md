# deamiDATE 1.0

deamiDATE 1.0 is a site-specific deamidation tool for palaeproteomics.

Given MaxQuant result output (in particular, the `evidence.txt` and `peptides.txt` files), it produces two graphs, one showing bulk deamidation, one showing site-specific deamidation.
Deamidation is calculated per experiment, and per protein in experiment. If any filtering of experiments or proteins is required, it is recommended this is done before running deamiDATE.

*Note: only tested on UNIX*

## Use:

* You will need Python 2.7
* You will need to install at least `numpy` and `pandas` (normally done with e.g. `pip install numpy`)

1. Make `run.py` executable, e.g. `chmod +x run`
2. Run program e.g. `./run.py` with the following args:
	1. A directory containing at least `evidence.txt` and `peptides.txt`, e.g. `./run Test_Data`
	2. [OPTIONAL] A text file containing a list (one per line) of relevant proteins to "filter in", e.g. `./run Test_Data Test_Data/protein_list.txt`

Results will appear in the directory in argument 1.

There are two result plots
* One site specific
* One "bulk"

There are three result csv files
* One with site-specific results
* One with bulk results, averaged per protein and per sample
* One with bulk results, but with each deamidation measurement reported (no averaging)

*Note: reported deamidation values are ALWAYS relative remaining N AND Q, NOT relative deamidation.*

## Common problems & solutions

Problems with evidence and peptide file path
* Change lines 22 and 23

Problems with finding `Info/asn.csv` or `Info/gln.csv`
* Change line 113 and 114

Filter out contaminants (default True)
* Change line 561 `filter_con = False`

Filter out reverse hits (default True)
* Change line 80 `if intensity != "" and "REV" not in protein:` to `if intensity != "":`

Show plots when created (default False)
* Change lines 563 and 564 `show = True`

Print csv of results (default True)
* Change lines 563 and 564 `to_print = False`

Change colormap
* In lines 384 - 390 change the colormap to another of your choice)

## Test Data

Test data is output of MaxQuant run on tapir files from:
Welker, F. et al. Ancient proteins resolve the evolutionary history of Darwin’s South American ungulates. Nature 522, 81–84 (2015).
