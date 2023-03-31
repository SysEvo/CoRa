# Feedback by active degradation (v02)
#   where only W is degraded in the complex
#   and active W in complex form

# Julia v.1.8

module mm
	# Required libraries
	using DifferentialEquations
	using ParameterizedFunctions

	# ODE system
	odeFB = @ode_def begin
		dY = (mY * (W + C)) - ((g + gY) * Y)
		dU = (mU * Y)       - ((g + gU) * U) - (eP * U * W) + ((e0 + gW + eM) * C)
		dW =    mW          - ((g + gW) * W) - (eP * U * W) + ((e0 + gU) * C)
		dC =                - ((g + eM) * C) + (eP * U * W) - ((gU + gW + e0) * C)
	end g mY gY mU gU mW gW e0 eP eM mUs;
	# ODE system without feedback
	odeNF = @ode_def begin
		dY  = (mY * (W + C)) - ((g + gY) * Y)
		dU  =  mUs           - ((g + gU) * U) - (eP * U * W) + ((e0 + gW + eM) * C)
		dW  =    mW          - ((g + gW) * W) - (eP * U * W) + ((e0 + gU) * C)
		dC  =                - ((g + eM) * C) + (eP * U * W) - ((gU + gW + e0) * C)
	end g mY gY mU gU mW gW e0 eP eM mUs;

	# Define system's output (total Y):
	function outFB(ss)
		return ss[1];
	end;
	function outNF(ss)
		return ss[1];
	end;

	# Define locally analogous system:
	function localNF(p,ss)
		p[:mUs] = p[:mU] * ss[1];
	end;
end
