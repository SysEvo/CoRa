# Brink motif feedback (v01)
#   with Y inducing B synthesis,

# Julia v.1.1.1

module mm
	# Required libraries
	using DifferentialEquations
	using ParameterizedFunctions

	# ODE system
	odeFB = @ode_def begin
		dY  = (mY * U)  - ((g + gY) * Y)
		dA  =  mA       - (g * A) - (eP * A * B) + (e0 * C) - (bA * A * Us)
		dB  = (mB * Y)  - (g * B) - (eP * A * B) + (e0 * C) - (bI * B * U)
		dC  =           - (g * C) + (eP * A * B) - (e0 * C) + (bA * A * Us)
		dU  =  mU       - (g * U) + (bA * A * Us) - (bI * B * U)
		dUs =          - (g * Us) - (bA * A * Us) + (bI * B * U)
	end g mY gY mA mB mU eP e0 bA bI mBs;
	# ODE system without feedback
	odeNF = @ode_def begin
		dY  = (mY * U)  - ((g + gY) * Y)
		dA  =  mA       - (g * A) - (eP * A * B) + (e0 * C) - (bA * A * Us)
		dB  = (mBs) - (g * B) - (eP * A * B) + (e0 * C) - (bI * B * U)
		dC  =           - (g * C) + (eP * A * B) - (e0 * C) + (bA * A * Us)
		dU  =  mU       - (g * U) + (bA * A * Us) - (bI * B * U)
		dUs =          - (g * Us) - (bA * A * Us) + (bI * B * U)
	end g mY gY mA mB mU eP e0 bA bI mBs;

	# Define system's output (total Y):
	function outFB(ss)
		return ss[1];
	end;
	function outNF(ss)
		return ss[1];
	end;

	# Define locally analogous system:
	function localNF(p,ss)
		p[:mBs] = p[:mB] * ss[1];
	end;
end