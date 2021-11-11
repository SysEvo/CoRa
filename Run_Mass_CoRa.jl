## Running in julia terminal
	cd("C:\\Users\\Lenovo\\Documents\\CoRa")
       using Pkg; 
       using BenchmarkTools;
       Pkg.activate(".");		# Activate local environment (requiere '.toml' files)
	iARG = (mm = "ATFv2",  # Label for motif file
       ex = "Fig2B",      # Label for parameters file
       pp = :mY,         # Label for perturbation type
       ax = :mY,         # Label for condition/environment
       an = "ExSSs");    # Chose analysis type (Options: ExSSs, ExDyn, CoRams, OptDY)
	include("Mass_CoRa_Main.jl");
