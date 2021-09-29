## Running in julia terminal
	cd("C:\\Users\\Lenovo\\Documents\\CoRa")
	iARG = (mm = "ATFv2",  # Label for motif file
       ex = "Fig2B",      # Label for parameters file
       pp = :mY,         # Label for perturbation type
       ax = :mY,         # Label for condition/environment
       an = "ExSSs");    # Chose analysis type (Options: ExSSs, ExDyn, DYms, OptDY)
	include("CoRa_Main_v2.jl")
