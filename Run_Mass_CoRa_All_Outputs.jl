## Running in julia terminal
	cd("C:\\Users\\ese_1\\OneDrive\\Documentos\\CoRa")
       using Pkg
       using CSV
       using DelimitedFiles
       Pkg.activate(".\\CoRa");		# Activate local environment (requiere '.toml' files)
	iARG = (mm = "ATFv2",  # Label for motif file
       ex = "1250Set2",      # Label for parameters file
       pp = :mY,         # Label for perturbation type
       ax = :mY)    # Label for condition/environment
       pars = CSV.File("InputFiles\\ARGS_ATFv2_Mass_Par_1250Set2.csv") # Core parameters
       superstrict = false
	include("Mass_CoRa_Main_ATF_All_Outputs.jl")