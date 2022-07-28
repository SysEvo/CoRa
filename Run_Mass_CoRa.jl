## Running in julia terminal
	cd("C:\\Users\\ese_1\\OneDrive\\Documentos\\CoRa")
       using Pkg;
       using CSV;
       using DelimitedFiles;
       Pkg.activate(".");		# Activate local environment (requiere '.toml' files)
	iARG = (mm = "BNFv1",  # Label for motif file
       ex = "1250Set1",      # Label for parameters file
       pp = :mY,         # Label for perturbation type
       ax = :mY);    # Label for condition/environment
       pars = CSV.File("InputFiles\\ARGS_BNFv1_Mass_Par_1250Set1.csv") # Core parameters
       # print = "short"      # Flag for the output file including/excluding the steady states (Options: "short", "all")
	include("Mass_CoRa_Main_BNF.jl");