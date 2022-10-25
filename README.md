## Git repository for ALSPAC LifeCycle green space and child
## cardiovascular outcomes analysis code (B3930)

This repository contains analysis scripts supporting the LifeCycle
project examining access to green space, and its interaction with
socioeconomic position, on child cardiovascular outcomes using a 
structured life-course modelling approach.

A description of each file is as follows:
 - LifeCycle_GreenSpace_SimulationScript.r: R script demonstrating how
 to apply these life course methods, with interactions between multiple
 explosure, in a simple simulated dataset
 - LifeCycle_GreenSpace_SimulationScript_Stata.do: Stata script showing
 how to apply these methods in a simple simulated dataset, but in Stata
 - LifeCycle_GreenSpace_SimulationStudyScript.r: R script to test and
 calibrate the methods of the formal simulation
 - Sims_HPC: Folder containing R and slurm shell scripts to run the formal
 simulation study, using the University of Bristol's High Performance
 Computing suite (https://www.bristol.ac.uk/acrc/high-performance-computing/)
 - LifeCycle_GreenSpace_CombiningHPCSimulations.r: R script to bring together
 and analyse the formal simulation results
 - LifeCycle_GreenSpace_ALSPAC_AnalysisScript.r: R script to run the applied
 example in a UK birth cohort (ALSPAC)

Note that ALSPAC data access is through a system of managed open access. 
Information about access to ALSPAC data is given on the ALSPAC website 
(http://www.bristol.ac.uk/alspac/researchers/access/). The datasets 
used in these scripts are linked to ALSPAC project number B3930; if you 
are interested in accessing these datasets, please quote this number 
during your application.
