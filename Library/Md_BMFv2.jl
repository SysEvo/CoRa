# Brink motif feedback (v02)
#   with Y represses A synthesis,

# Julia v.1.8

module mm
	# Required libraries
	using DifferentialEquations
	using ParameterizedFunctions

	# ODE system
	odeFB = @ode_def begin
		dY  = (mY * U)              - ((g + gY) * Y)
		dA  = (mA * (kD/(kD + Y)))  - (g * A) - (eP * A * B) + (e0 * C) - (bA * A * Us)
		dB  =  mB                   - (g * B) - (eP * A * B) + (e0 * C) - (bI * B * U)
		dC  =                       - (g * C) + (eP * A * B) - (e0 * C) + (bA * A * Us)
		dU  =  mU                   - (g * U) + (bA * A * Us) - (bI * B * U)
		dUs =                       - (g * Us) - (bA * A * Us) + (bI * B * U)
	end g mY gY mA kD mB mU eP e0 bA bI mAs;
	# ODE system without feedback
	odeNF = @ode_def begin
		dY  = (mY * U)              - ((g + gY) * Y)
		dA  =  mAs                  - (g * A) - (eP * A * B) + (e0 * C) - (bA * A * Us)
		dB  =  mB                   - (g * B) - (eP * A * B) + (e0 * C) - (bI * B * U)
		dC  =                       - (g * C) + (eP * A * B) - (e0 * C) + (bA * A * Us)
		dU  =  mU                   - (g * U) + (bA * A * Us) - (bI * B * U)
		dUs =                       - (g * Us) - (bA * A * Us) + (bI * B * U)
	end g mY gY mA kD mB mU eP e0 bA bI mAs;

	# Define system's output (total Y):
	function outFB(ss)
		return ss[1];
	end;
	function outNF(ss)
		return ss[1];
	end;

	# Define locally analogous system:
	function localNF(p,ss)
		p[:mAs] = p[:mA] * (p[:kD]/(p[:kD] + ss[1]));
	end;
end
