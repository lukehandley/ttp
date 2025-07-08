#!/bin/bash

# Activate the conda environment
source ~/miniforge3/bin/activate astroq

# Run packed example
python bin/runTTP.py -f examples/packed/inputs.txt

# Run packed-half-1 example
python bin/runTTP.py -f examples/packed-half-1/inputs.txt

# Run packed-half-2 example
python bin/runTTP.py -f examples/packed-half-2/inputs.txt 