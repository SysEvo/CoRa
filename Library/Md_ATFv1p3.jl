# Antithetic feedback (v01 p03)
#   with inactive W in complex form
#   considering a more complex process

# Julia v.1.8

module mm
	# Required libraries
	using DifferentialEquations
	using ParameterizedFunctions

	# ODE system
	odeFB = @ode_def begin
		dY  = (mY * Y1)                         - ((g + gY) * Y)
		dU  = (mU * Y)                          - ((g + gU) * U) - (eP * U * W) + ((e0 + gW) * C)
		dW  =    mW                             - ((g + gW) * W) - (eP * U * W) + ((e0 + gU) * C)
		dC  =                                   - ((g + eM) * C) + (eP * U * W) - ((gU + gW + e0) * C)
		dY0 = (m0 * W)                          - ((g + gY) * Y0)
		dY1 = (m1 * Y0) + (mP * (Y1/(Y1 + kP))) - ((g + gY) * Y1)
	end g mY gY mU gU mW gW e0 eP eM m0 m1 mP kP mUs;
	# ODE system without feedback
	odeNF = @ode_def begin
		dY  = (mY * Y1)                         - ((g + gY) * Y)
		dU  = (mUs)                             - ((g + gU) * U) - (eP * U * W) + ((e0 + gW) * C)
		dW  =    mW                             - ((g + gW) * W) - (eP * U * W) + ((e0 + gU) * C)
		dC  =                                   - ((g + eM) * C) + (eP * U * W) - ((gU + gW + e0) * C)
		dY0 = (m0 * W)                          - ((g + gY) * Y0)
		dY1 = (m1 * Y0) + (mP * (Y1/(Y1 + kP))) - ((g + gY) * Y1)
	end g mY gY mU gU mW gW e0 eP eM m0 m1 mP kP mUs;

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