## Running in julia terminal
	cd("C:\\Users\\ese_1\\OneDrive\\Escritorio\\CoRa")
       using Pkg;
       using CSV;
       using DelimitedFiles;
       Pkg.activate(".\\CoRa");		# Activate local environment (requiere '.toml' files)
	iARG = (mm = "FADv2",  # Label for motif file
       ex = "1250Set2",      # Label for parameters file
       pp = :mY,         # Label for perturbation type
       ax = :mY);    # Label for condition/environment
       pars = CSV.File("InputFiles\\ARGS_FADv2_Mass_Par_1250Set2.csv"); # Core parameters
       key_names = (:g, :mY, :gY, :mU, :gU, :mW, :gW, :e0, :eP, :eM, :mUs); #The names of the parameters of the systems
       strict = true;  # Should a steady state be found, and then refound in the next iteration, and this difference not be in accordance to the rtol given, said SS will be obtained regardless
       gap_size = 20; # This is the size of a region to be examined for NaN gaps
       gap_tol = 0.5; # How much of the examined region can be NaNs before the program kills it
       rtol = 1e-10; # This is the error tolerance for the differential equations solver to know it has reached a satisfactory steady state
	include("Mass_CoRa_Main.jl");
