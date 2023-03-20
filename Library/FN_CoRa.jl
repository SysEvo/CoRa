# CoRa related functions
#	SS - Steady state function for a given system
#	?
#
# Julia v.1.8

module fn
	# Required libraries
	using DelimitedFiles
	using Distributions
	using DifferentialEquations

	# AllNaN - Assign NaN values to output vectors
	# INPUT: a,b - Vectors to turn into all NaN
	# OUPUT: a,b - All NaN vectors
	function AllNaN(a, b)
		return a.+NaN, b.+NaN
	end

	# NaNCheck - Check for NaN values in the reference vector, and apply AllNan() 
	#            in both vectors if any NaN is found.
	# INPUT: a - Reference vector to be checked for NaN values
	#		 b - Additional output vector to be turned in all NaN values
	# OUPUT: a,b
	function NaNCheck(a, b)
		if(any(isnan.(a)))
			a, b = fn.AllNaN(a, b);
			return a, b
		else
			return "Valid"
		end
	end

	# SS - Steady state function for a given system
	# INPUT: syst - Handle for the ODE system (@ode_def)
	#        p    - Dictionary function with the ODE parameters & values
	#        x0   - Vector of initial state of the ODE system
	#        rtol - Tolerance value for ODE solver
	# OUPUT: ss   - Vector of steady state of the ODE system
	function SS(syst, p, x0, rtol)
		x = copy(x0);
		pV = [p[eval(Meta.parse(string(":",i)))] for i in syst.sys.ps];	# Parameter values as a vector.
		iS = 0;			# Accumulated simulation time.
		dX = rtol + 1;	# Maximum relative change in the ODE simulation. (Note: arbitrary starting value).
		# Simulate ODE system until the error tolerance is reached:
		while(dX > rtol)
			ss = try
				# Try using standard ODE solver:
				fn.solve(fn.ODEProblem(syst,x,1e6,pV); reltol=rtol);
			catch
			end
			if ss.retcode != :Success
				try
					# Try using stiff ODE solver:
					ss = fn.solve(fn.ODEProblem(syst,x,1e6,pV), alg_hints=[:stiff]; reltol=rtol);
				catch err
					# If both failed, return an error message and NaN vextor:
					println("WARNING: Error in ODE simulation: <<",err,">>. ss --> NaN")
					x = zeros(length(syst.syms)).+NaN;
					break
				end
			end;
			# Calculate maximum relative change in the ODE simulation:
			dX = maximum(abs.(big.(ss(1e6))-big.(ss(1e6-0.01)))./big.(ss(1e6)));
			# Update vector of current state of the ODE system:
			x = ss(1e6);
			# Count iterations & stop if larger than:
			iS += 1;
			if(iS>10)
				println("WARNING: Maximum iteration reached (10 times). Max relative Delta: ",dX)
				break
			end
		end
		return x
	end;

	# Check - Check that the NF system is truly locally analogous to FB
	# INPUT: ssFB - Steady state of the feedback (original) system
	#        ssNF - Steady state of the no-feedback (analogous) system
	#        rtol - Tolerance value for ODE solver
	#        syst - Handle for the ODE system (@ode_def)
	# OUPUT: 
	function Check(ssFB, ssNF, rtol, syst)
		if(any.(isnan.(syst.outFB(ssFB))) || any.(isnan.(syst.outFB(ssNF))) || (abs(syst.outFB(ssFB) - syst.outNF(ssNF)) > 1e-4))
			rtol *= 1e-3
			if(rtol < 1e-24)
				println("ERROR: Check NF system (reltol=",rtol,").")
				println(vcat(pert.p,i,[p[eval(Meta.parse(string(":",i)))] for i in syst.sys.ps],syst.outFB(ssFB),syst.outNF(ssNF)))
				if(abs(syst.outFB(ssFB) - syst.outNF(ssNF))/syst.outFB(ssFB) > 0.01)
					ssFB, ssNF = AllNaN(ssFB, ssNF);
					println("Error too large. SS results excluded!")
				end
			end
			return ssFB, ssNF, rtol, "Insufficient"
		else
			return ssFB, ssNF, rtol, "Sufficient"
		end
	end;

	# SSandCheck - Call SS() and Check() accordingly.
	# INPUT: p    - Dictionary function with the ODE parameters & values
	#        x0   - Vector of initial state of the ODE system
	#        rtol - Initial tolerance value for ODE solver
	#        syst - Handle for the ODE system (@ode_def)
	# OUPUT: 
	function SSandCheck(p, x0, rtol, syst)
		flag = "Insufficient"
		ssFB, ssNF = fn.AllNaN(x0, x0);
		while(rtol >= 1e-24 && flag == "Insufficient")
			# Reference steady state:
			ssFB = fn.SS(syst.odeFB, p, x0, rtol);
			if(fn.NaNCheck(ssFB, ssNF) != "Valid")
				println("Condition excluded! ssFB --> NaN")
				break
			end
			# Locally analogous system reference steady state:
			syst.localNF(p,ssFB);
			###Of note here, instead of using the initial condition x0, we use ssFB.
			ssNF = fn.SS(syst.odeNF, p, ssFB, rtol);
			if(fn.NaNCheck(ssNF, ssFB) != "Valid")
				println("Condition excluded! ssNF --> NaN")
				break
			end
			ssFB, ssNF, rtol, flag = Check(ssFB, ssNF, rtol, syst)
		end
		return ssFB, ssNF, rtol
	end

	# Perturbation - Apply perturbation and recalculate steady states
	# INPUT: ssR  - Steady state of the feedback (original) system before perturbation (reference)
	#        soR  - Steady state of the no-feedback (analogous) system before perturbation (reference)
	#        p    - Dictionary function with the ODE parameters & values
	#        rtol - Tolerance value for ODE solver
	#        syst - Handle for the ODE system (@ode_def)
	#        pert  - Handle for the perturbation details
	# OUPUT: ssD  - Steady state of the feedback (original) system after perturbation (disturbed)
	#        soD  - Steady state of the no-feedback (analogous) system after perturbation (disturbed)
	function Perturbation(ssR, soR, p, rtol, mm, pert)
		p[pert.p] *= pert.d;
		ssD = fn.SS(mm.odeFB, p, ssR, rtol);
		soD = fn.SS(mm.odeNF, p, soR, rtol);
		p[pert.p] /= pert.d;
		if(fn.NaNCheck(ssD, soD) != "Valid")
			println("Condition excluded! ssD --> NaN")
		elseif(fn.NaNCheck(soD, ssD) != "Valid")
			println("Condition excluded! SoD --> NaN")
		end
		return ssD, soD
	end

	# CoRa function
	# INPUT: ssR  - Output of full system before perturbation
	#        ssD  - Output of full system after perturbation
	#        soR  - Output of non-feedback system before perturbation
	#        soD  - Output of non-feedback system after perturbation
	# OUPUT:      - CoRa value
	function CoRa(ssR, ssD, soR, soD)
		if abs(log10(soD/soR)) < 1e-4
			return NaN
		end
		return log10(ssD/ssR)/log10(soD/soR)
	end


	# ODE dynamics for a given system
	# INPUT: syst  - Handle for the ODE system (@ode_def)
	#        p     - Dictionary function with the ODE parameters & values
	#        x0    - Vector of initial state of the ODE system
	#        tspan - Time to simulate
	# OUPUT: xD    - Vector of steady state of the ODE system
	function Dyn(syst, p, x0, tspan)
		pV = [p[eval(Meta.parse(string(":",i)))] for i in syst.sys.ps];	# Parameter values as a vector.
		xD = try
			# Try using standard ODE solver:
			fn.solve(fn.ODEProblem(syst,x0,tspan,pV); reltol=1e-6);
		catch
		end
		if xD.retcode != :Success
			try
				# Try using stiff ODE solver:
				ss = fn.solve(fn.ODEProblem(syst,x,1e6,pV), alg_hints=[:stiff]; reltol=rtol);
			catch err
				# If both failed, return an error message and NaN vextor:
				println("WARNING: Error in ODE simulation: <<",err,">>. xD --> NaN")
				xD = zeros(length(syst.syms)).+NaN;
			end
		end;
		return xD
	end;

	# CoRa curve
	# INPUT: p     - Dictionary function with the ODE parameters & values
	#        pert  - Handle for the perturbation details
	#        motif - Handle for the considered motif
	#        x0    - Vector of initial state of the ODE system
	# OUPUT: CoRas - Vector of CoRa values for the range of parameters
	function CoRac(p, pert, mm, x0)
		r = 10 .^ collect(pert.r[1]:pert.s:pert.r[2]);
		CoRas = ones(length(r)) .+ Inf;
		for i in 1:lastindex(r)
			p[pert.p] *= r[i];
			ssR, soR, rtol = SSandCheck(p, x0, 1e-12, mm)
			ssD, soD = fn.Perturbation(ssR, soR, p, rtol, mm, pert)
			CoRas[i] = fn.CoRa(mm.outFB(ssR), mm.outFB(ssD), mm.outNF(soR), mm.outNF(soD))
			p[pert.p] /= r[i];
		end
		return CoRas
	end;

	# CoRa "metrics"
	# INPUT: CoRas - Vector of CoRa values for the range of parameters
	#        pert  - Handle for the perturbation details
	# OUPUT:       - [|CoRa<=eps|, min(CoRa)]
	function CoRam(CoRas, pert)
		return [sum(CoRas .<= pert.eps)/length(CoRas), minimum(filter(!isnan,CoRas))]
	end;
end