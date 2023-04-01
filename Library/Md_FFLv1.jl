# Feedback with feedforward loop (v01)

# Julia v.1.8

module mm
	# Required libraries
	using DifferentialEquations
	using ParameterizedFunctions

	# ODE system
	odeFB = @ode_def begin
		dY  = (mY * (U + W))       - ((g + gY) * Y)
		dU  = (mU * (kD/(Y + kD))) - ((g + gU) * U)
		dW  = (mW * U)             - ((g + gW) * W)
	end g mY gY mU kD gU mW gW mUs;
	# ODE system without feedback
	odeNF = @ode_def begin
		dY  = (mY * (U + W))       - ((g + gY) * Y)
		dU  =    mUs               - ((g + gU) * U)
		dW  = (mW * U)             - ((g + gW) * W)
	end g mY gY mU kD gU mW gW mUs;

	# Define system's output (total Y):
	function outFB(ss)
		return ss[1];
	end;
	function outNF(ss)
		return ss[1];
	end;

	# Define locally analogous system:
	function localNF(p,ss)
		p[:mUs] = p[:mU] * (p[:kD]/(ss[1] + p[:kD]));
	end;
end