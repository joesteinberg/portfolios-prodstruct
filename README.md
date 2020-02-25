# portfolios-prodstruct
Supplementary materials for "International Portfolio Diversification and the Structure of Global Production" (RED, 2018)
Joseph B. Steinberg
University of Toronto
RED Manuscript RED-16-74R2

This document describes the data and programs needed to replicate my empirical and quantitative
analysis. The document is organized into three sections:

0. Datasets
1. Python scripts (data processing and figure creation)
2. Matlab programs (quantitative analysis and simple example)

The python scripts (section 1) should be run before the Matlab programs (section 2) in order to 
replicate the quantitative analysis. All of my analyses were performed on a Dell Workstation 
running Ubuntu Linux 12.04. All scripts and programs should work on any platform without
modification as long as the user has the requisite software packages (all of which are standard
in the economics commmunity). There are no random numbers generated in any of the parts of the 
analyses, so no seeds are needed.

///////////////////////////////////////////////////////////////////////////////////////////////

0. Datasets ("data_sources")

The two datasets below are contained in the supplementary materials. They are widely used
in the literature so I will not discuss them here.

0.1 Lane and Milesi-Feretti ("data_sources/ewn19702011.csv")
0.2 Penn World Tables 8.0 ("data_sources/pwt80.dta")

The World Input Output Database (WIOD) data that I use to calibrate the model require a lot of
storage space (20 gigabytes or so). I have stored them in HDF format. I have included the pythons
scripts that process this data in the suppplementary materials, but I have not included the data
themselves. They are available upon request, and I have included all intermediate files created by
these processing scripts so that the user can replicate all parts of my analysis.

///////////////////////////////////////////////////////////////////////////////////////////////

1. Python scripts ("programs/python")

Note: to run all python scripts, the user must have Numpy, Pandas, and Matplotlib installed.
These come standard with the Anaconda distribution which is common in scientific computing.
All latex output is written to "programs/python/latex." None of these scripts take more than
a minute or so to run on a standard computer (with the exception of 1.1 and 1.2, which the 
user does not need to run to replicate my work).

1.0 The script "regions.py" contains a few helper objects that are used to aggregate data
according to the scheme in Table 1. All other python scripts use this one, but the user
does not need to call it directly.

1.1 The script "prepare_wiod_data.py" loads the raw data from the World Input Output Database,
aggregates it according to the scheme laid out in table 1, and saves some intermediate files.
The raw data is over 20 gigabytes so I have not included it in the supplementary materials.
I have included the intermediate files (in the directory "programs/python/wiod_data") so you
do not need to run this script.

1.2 The script "proddata_regions.py" loads the raw data from the WIOD and computes some statistics
about the global production structure, including the intermediate share of trade used in Figure 1.
This script writes an intermediate file, "programs/python/results/proddata_regions.txt" which is
used by a later python script to create the figure. Again, I have not included the raw data
but I have included the intermediate file.

1.3 The script "iomats.py" loads the intermediate files created in 1.1 and constructs the
input-output tables used in the analysis.
- First it loads the intermediate data from 1.1 and creates the benchmark 1995 and 2011 input-output tables
- At this point, I manually loaded these into the Excel file "excel/4country_alt_scenarios.xlsm" and used a
  VBA macro to construct the counterfactual input-output tables. The online appendix describes these steps,
  and the supplemenbtary materials contains all the results of these steps.
- The user should now re-run iomats.py which will load the results of the Excel work, write latex tables 
  with the benchmark and counterfactual IO tables (Tables 3 and 4), and write text fies used by the matlab
  script "calibrate.m" below.

1.4 The script "iomats_stats.py" uses the input-output data from 1.2 and creates Table 5.

1.5 The script "portfolios.py" loads the Lane and Milesi-Feretti and Penn World Tables data and 
calculates country-level internatoinal portfolio diversification in each year between 1995 and 2011.
See the online appendix for details.

1.6 The script "plots_regions.py" uses the intermediate files from 1.2 and 1.5 to create Figure 1. It also
writes a text file with the region-level diversification data used by the matlab program below to calibrate
the wedges and write the latex tables.


///////////////////////////////////////////////////////////////////////////////////////////////

2. Matlab programs ("programs/matlab/")

Note: To run the Matlab code used in this paper you must have the symbolic toolbox and the optimization
toolbox. The quantitative model takes approximiately 10 minutes to run on my Dell Workstation with 
dual 6-core Intel Xeon CPUs. I have not used any parallelizations, though, so this runtime should not
vary significantly for other setups.

2.0.0 Uribe and Schmidt-Grohe Toolkit ("programs/matlab/usg_toolkit")

This directory contains scripts written by Martin Uribe abd Stephanie Schmidt-Grohe that help
set up the symbolic equilibrium system and take derivatives. The user does not need to interface
with these at all.

2.0.1 Devereux and Sutherland (2011) example code ("programs/matlab/devereux_sutherland_example")

This directory contains the code to solve one of the examples from Devereux and Sutherland (2011).
It contains some helper scripts that are used in solving for the linearized decision rules. The
user does not need to interface with these at all.

2.1 Quantitative model ("programs/matlab/quantitative_model")

The script "quant_model.m" is the only script the user needs to run to reproduce the quantitative
analysis. This script calls a number of other scripts directly:
- "analytical_params_vars_eqns.m" sets up the model variables and equations as symbolics
- "analytical_derivs.m" taes symbolic derivatives
- "calibration.m" calibrates model parameters to input-output data from step 1
- "zero_order_portfolios.m" calibrates portfolio wedges and solves for steady-state portfolios
- "calibration_latex1.m" and "calibration_latex2.m" write parameters to latex tables (tables 2 and 6)
- "results_latex.m" writes results to latex tables (tables 7-10)
Some of these scripts call others as well. Notable ones are:
- "fileio.m" loads the calibration data
- "store_ss.m"
- "ls2ds.m" and "ls_solution_js.m" construct the state-space representation of the equilibrium system
  and solve for linearized decision rule coefficients (see Devereux and Sutherland, 2011)...
  these scripts call some of the helper files from 2.0.1 above
- "zero_order_portfolios2.m" and "share_func.m" are helper functions for solving for portfolios

2.2 Lucas Tree model ("programs/matlab/lucas_tree_model")

Tihs directory contains scripts to solve for steady state portfolios in a simple Lucas tree model
with an arbitrary number of countries. This example demonstrates how to use my extension of the DS
method without all the extra stuff needed in the quantitative model (input-output structure, etc.)
The user just needs to run the script "lucas_trees.m." The other scripts are similar to those in
2.1.
