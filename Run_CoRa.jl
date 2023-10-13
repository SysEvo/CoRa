## Running in julia terminal
	#cd("C:\\MyLibrary\\Location\\")
	using Pkg; Pkg.activate(".");
	iARG = (mm = "ATFv2p3",		# Label for motif file
			ex = "Fig2",		# Label for parameters file
			pp = :mY,			# Label for perturbation type
			ax = :mY,			# Label for condition/environment
			an = "ExSSs");		# Chose analysis type (Options: ExSSs, ExDyn, CoRams, OptCoRa)
	include("CoRa_Main.jl")
