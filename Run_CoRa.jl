## Running in julia terminal
       cd("C:\\Users\\ese_1\\OneDrive\\Documentos\\CoRa")
       #using Pkg; 
       #using BenchmarkTools
       #Pkg.activate(".");		# Activate local environment (requiere '.toml' files)
	iARG = (mm = "ATFv2",  # Label for motif file
       ex = "Fig2B",      # Label for parameters file
       pp = :mY,         # Label for perturbation type
       ax = :mY,         # Label for condition/environment
       an = "ExDyn");    # Chose analysis type (Options: ExSSs, ExDyn, CoRams, OptDY)
	include("CoRa_Main_v2.5.jl");
