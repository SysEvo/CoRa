## Running in julia terminal
	cd("C:\\Users\\ese_1\\OneDrive\\Documentos\\CoRa")
       using Pkg;
       using CSV;
       using DelimitedFiles;
       Pkg.activate(".\\CoRa");		# Activate local environment (requiere '.toml' files)
	iARG = (mm = "ATFv2",  # Label for motif file
       ex = "1250Set2",      # Label for parameters file
       pp = :mY,         # Label for perturbation type
       ax = :mY);    # Label for condition/environment
       pars = CSV.File("InputFiles\\ARGS_ATFv2_Mass_Par_1250Set2.csv"); # Core parameters
       strict = true  # Should a steady state be found, and then refound in the next iteration, and this difference not be in accordance to the rtol given, said SS will be obtained regardless
       gap_tol = 10 # This is how wide we accept gaps to be in lines
	include("Mass_CoRa_Main_ATF.jl");
