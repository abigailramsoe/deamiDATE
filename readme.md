# deamiDATE 1.0

deamiDATE 1.0 is a site-specific deamidation tool for palaeproteomics.

Given MaxQuant result output (in particular, the evidence.txt and peptides.txt files), it produces two graphs, one showing bulk deamidation, one showing site-specific deamidation.
Deamidation is calculated per experiment, and per protein in experiment. If any filtering of experiments or proteins is required, it is recommended this is done before running deamiDATE.

## Use:

### Unix:
1. Download files from GitHub
2. Unzip
3. Navigate to deamiDATE directory
4. Run program, with include an argument for the folder with MQ result files,
e.g. `./run ../Test Data/Tapir` must be a directory containing the evidence.txt and peptides.txt

Results will appear in the directory given in step 5.


### Windows

*Note: untested*

Steps 1 - 2 as per Unix instructions, then

1. Open command prompt (no idea)
2. `run C:\Users\Username\Test Data\Tapir` where the second argument is the directory containing the MQ results

Results will *hopefully* appear in the given MQ results directory.

## Frequent problems

*To be added after testing!*
