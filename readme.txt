Run the following command in MATLAB to 

addpath(genpath(cat(2, pwd, foldername)))

Replace foldername with the name of the folder containing the code, and pwd with the directory
containing the folder.

Eg. addpath(genpath(cat(2, 'C:\Users\User\Downloads', '\swarm_code')))


Short explanation of what the files do below

C01.m, C02.m etc: Objective functions used for the numerical experiments
gen_transrot.m : used to generate the rotation matrices and/or translation vectors where relevant. They should be stored in the "data" folder. Some examples generated are provided.
mPSO.m, hmPSO.m: Code of the respective modified PSO algorithms. Similarly for CSO and BAT.
hmSwarmRunner.m: Code to setup and run the algorithms for comparison, including the saving of output / results.

Some driver files previously used are provided, eg. pso_c04.m. You can also see STAGING.m for another example.
The driver files can be used to provide and keep track of some parameters for the experimental setup. (Useful if running multiple experiments with slightly different parameters)