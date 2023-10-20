# CoRa related functions
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
			if(any(x.<0))
				println("WARNING: Error in ODE simulation (negative values obtained). ss --> NaN")
				x = zeros(length(syst.syms)).+NaN;
				break
			end
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
				#println(vcat(i,[p[eval(Meta.parse(string(":",i)))] for i in syst.sys.ps],syst.outFB(ssFB),syst.outNF(ssNF)))
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
	#        mm    - Handle for the considered motif
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
	# INPUT: p     - Dictionary function with the ODE parameters & values
	#        pert  - Handle for the perturbation details
	#        mm    - Handle for the considered motif
	#        x0    - Vector of initial state of the ODE system
	# OUPUT: CoRas 		 - Vector of CoRa values for the range of parameters
	#   	 |CoRa<=eps| - Proportion of the CoRas vector less or equal to eps
	#		 min(CoRa)	 - Minimum CoRa value in CoRas
	function CoRam(p,pert,mm,x0)
		CoRas = fn.CoRac(p,pert,mm,x0);  # Calculate CoRa curve
		return [CoRas,sum(CoRas .<= pert.eps)/length(CoRas), minimum(filter(!isnan,CoRas))]
	end;

	# MRW random initial parameters
	# INPUT: mrw - Handle for the optimization (MRW) details
	#        p   - Dictionary function with the ODE parameters & values
	# OUPUT:     - Updated p
	function mrwR(mrw,p)
		for i in 1:length(mrw.pOp)
			p[mrw.pOp[i]] = 10 .^ (rand(Uniform(mrw.pMin[i], mrw.pMax[i])));
		end
	end

	# MRW iteration
	# INPUT: mrw   - Handle for the optimization (MRW) details
	#        p     - Dictionary function with the ODE parameters & values
	#        pert  - Handle for the perturbation details
	#        mm    - Handle for the considered motif
	#        x0    - Vector of initial state of the ODE system
	#		 CoRa0 - Vector of CoRa values for the range of parameters in reference system
	#   	 op0   - Proportion of the CoRas vector less or equal to eps in reference system
	#		 mi0   - Minimum CoRa value in CoRas in reference system
	# OUPUT: CoRa0 - Vector of CoRa values for the range of parameters in updated reference system
	#   	 op0   - Proportion of the CoRas vector less or equal to eps in updated reference system
	#		 mi0   - Minimum CoRa value in CoRas in updated reference system
	function mrwI(mrw,p,pert,mm,x0,CoRa0,op0,mi0)
		# Choose new parameter values:
		r0 = zeros(length(mrw.pOp));			# Vector of parameters to optimize in reference system
		rI = rand(MvNormal(zeros(length(mrw.pOp)), zeros(length(mrw.pOp)) .+ mrw.cov)); # Random values to update parameters
		for pI in 1:length(mrw.pOp)
			r0[pI] = p[mrw.pOp[pI]];			 # Save previous value
			p[mrw.pOp[pI]] *= (mrw.M .^ rI[pI]); # Update value
			# Exclude values outside regime of exploration:
			if p[mrw.pOp[pI]] < (10.0 ^ mrw.pMin[pI])
				p[mrw.pOp[pI]] = (10.0 ^ mrw.pMin[pI])
			elseif p[mrw.pOp[pI]] > (10.0 ^ mrw.pMax[pI])
				p[mrw.pOp[pI]] = (10.0 ^ mrw.pMax[pI])
			end
		end
		# Calculate new CoRa and metrics:
		CoRa1,op1,mi1 = fn.CoRam(p,pert,mm,x0);  # New values of properties to optimize, proportion of CoRas<=eps and min(CoRas) value)
		# Evaluate if accept new parameter values or not:
		## Only accept in the regime of interest, i.e. CoRa>=0:
		c1 = (mi1>=0);
		## If CoRa>eps for all conditions, evaluate the min(CoRa) for both sets:
		### NOTE: As mi0,mi1=[0,1], correct exponential with the expected variance of ~U(0,1)
			xiC = (mi0 ^ 2) / (2 * 0.083);
			xiP = (mi1 ^ 2) / (2 * 0.083);
		c2 = (minimum([op0,op1])==1.0) && (rand() < exp(xiC - xiP));
		## If CoRa>=eps for some conditions, evaluate the |CoRa<=eps| for both sets:
		### NOTE: As op0,op1=[0,1], correct exponential with the expected variance of ~U(0,1)
			xiC = ((1 - op0)^2) / (2 * 0.083);
			xiP = ((1 - op1)^2) / (2 * 0.083);
		c3 = rand() < exp(xiC - xiP);
		if(c1 && (c2 || c3))
			# If yes, update "reference" system
			CoRa0 = CoRa1;
			op0 = op1;
			mi0 = mi1;
		else
			# If not, revert to previous parameter values
			for pI in 1:length(mrw.pOp)
				p[mrw.pOp[pI]] = r0[pI];
			end
		end
		return CoRa0,op0,mi0
	end;
end