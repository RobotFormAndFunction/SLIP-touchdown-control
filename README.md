# SLIP-touchdown-control

This repository should contain all files needed to create data sets of velocity  magnitudes as a function of touchdown angle and velocity angle for SLIP like runners.  

Please reach out if files are missing or anything is broken, I may have missed a file.

All files will need to be added to the search path for the program to run

4 Main parts to this repository:

**1) Data set creation**

Run Centered_data_creation to creat the data set.  You can set the resolution and bounds of the touchdown and velocity angles.  As configured this calls a mex file to create the data set as it reduces runtime signifcantly.  The compiled mex file does the same as the included data_set_centered.

**2) Realtive stiffness relationship**

This creates the mapping and data sets for the determining the slope and intercept of the froude number and relative stiffness.  Run froude_vary_krel

**3) Perturbation recovery**

Run td_pertabation_calc (I didn't know how to spell perturbation when I made the file name)

This allows you to use the lookup tables with varying resolution, create the optimal solution, and apply noise to the optimal solution

**4) Rough terrain**

To use touchdown control on rough terrain run: rough_terrain_LU_control to get results

This allows you to use the lookup tables with varying resolution, create the optimal solution, and apply noise to the optimal solution on rough terrrain
