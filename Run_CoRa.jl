## Running in julia terminal
	cd("C:\\Users\\Lenovo\\Documents\\CoRa")
       #using Pkg; 
       #Pkg.activate(".");		# Activate local environment (requiere '.toml' files)
	iARG = (mm = "ATFv2",  # Label for motif file
       ex = "FigS2",      # Label for parameters file
       pp = :mY,         # Label for perturbation type
       ax = :mY,         # Label for condition/environment
       an = "CoRams");    # Chose analysis type (Options: ExSSs, ExDyn, CoRams, OptDY)
	include("CoRa_Main_v2.5.jl")