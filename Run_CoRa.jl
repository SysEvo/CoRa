## Running in julia terminal
	#cd("C:\\MyLibrary\\Location\\")
	using Pkg; Pkg.activate(".");
	iARG = (mm = "FADv1",		# Label for motif file
			ex = "Fig3",		# Label for parameters file
			pp = :mY,			# Label for perturbation type
			ax = :mY,			# Label for condition/environment
			an = "OptCoRa");		# Chose analysis type (Options: ExSSs, ExDyn, CoRams, OptCoRa)
	include("CoRa_Main.jl")
