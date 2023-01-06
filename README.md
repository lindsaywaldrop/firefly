# Modeling odor capture by fireflies

This project models odor capture by beetle antennae in the family Lampyridae. 

### Instructions for Use 

The modeling consists of three phases: 

 1. Production of 2D antenna models based on morphometrics on 15 species of lampyrid beetle antennae. These individual models are produced in src/Rscripts/Species_models. Source each species script to generate the appropriate files to run morphology model in IBAMR and MATLAB.
 
 2. Simulation of steady-state fluid flow at low Reynolds numbers of the 2D models using the Immersed Boundary Method with Adaptive Mesh Refinement (IBAMR, https://ibamr.github.io/). This code is located in src/ibamr. Edit lines 6 and 7 of the Makefile to include the paths to your build of IBAMR. In shell, run the compile_ibamr_sniff.sh file in src/bash using: <br />
 `$ sh compile_ibamr_sniff.sh $WD $Species` <br />
 where WD is your working directory (top level) and Species is the species you wish to run. Next, create input2d files by running the src/set_input2d.sh file with: <br />
 `$ sh set_input2d.sh $WD $Species` <br />
 You may then run IBAMR in src/bin using:  <br />
 `$ ./main2d input2dn` <br />
 where n is the number of the simulation for that species you would like to complete. Scripts src/bash/runkeckibamrjobs.sh and runibamr.job provides code to automate running jobs on SLURM-based clusters. 
 
 3. Simulation of odor convection and capture in MATLAB using the IBAMR fluid fields based on code in the publication Waldrop et al. 2018 (https://doi.org/10.1007/s10886-018-1017-2). This code is located in src/matlab.
 You may run single simulations in the MATLAB console using:<br />
 `> entsniff(WD,Species,i,'water',n)`<br />
 where i is the simulation number(s) and n is the number of runs to use in parallel. (If you do not have parallel computing set up on MATLAB, use 1. If using parallel, the number of runs in i should be equal to n.)
 4. Additional analysis of values can be found in the R Markdown notebook in doc/Odor-ResultsWorkbook.Rmd. This workbook can be used to run analysis on results of each species using the code in src/Rscripts. 
 
### Software Requirements

This code base is compatible with IBAMR v 0.1.11, R version 4.1.1, and MATLAB release 2019b. 

In R, install required R packages with: <br />
`> source("./src/Rscripts/install-R-packages.R")`

### Contact

For comments or questions, please contact waldrop@chapman.edu.




