# BayesianOpt4SurfaceMatching


This code optimize the geometry of inorganic interface crystals using automated DFT calculations via Bayesian Optimization approach.

## Requirements ##

1. Python 3.6+
2. NumPy
3. Pandas
4. ASE (https://wiki.fysik.dtu.dk/ase/)
5. pymatgen (https://pymatgen.org/)
6. bayesian-optimization https://github.com/fmfn/BayesianOptimization
7. Vienna Ab initio Simulation Package (VASP) https://www.vasp.at/

## Set up the input file (input.json) and the POSCAR files before running the code 

The input files include:
- The structure of interface, substrate, and film slabs.
- input.json: Includes general flags required in the VASP calculation.
The flags can be added or removed. More flag keys can be found in the ASE VASP calculator.

## Installation
* `python ./setup.py develop`

## Usage
Before running, change the environment variables VASP_RUN_COMMAND, OUTFILENAME, and VASP_PP_PATH.

* `python ./main.py`

