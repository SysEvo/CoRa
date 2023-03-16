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

	function CoRa(ssR, ssD, soR, soD)
		if abs(log10(soD/soR)) < 1e-4
			return NaN
		end
		return log10(ssD/ssR)/log10(soD/soR)
	end


	# ODE dynamics for a given system
	# INPUT: syst - Handle for the ODE system (@ode_def)
	#        p    - Dictionary function with the ODE parameters & values
	#        x0   - Vector of initial state of the ODE system
	#        tspan- Time to simulate
	# OUPUT: xD   - Vector of steady state of the ODE system
	function Dyn(syst, p, x0, tspan)
		pV = [p[i] for i in syst.params];
		xD = solve(ODEProblem(syst,x0,tspan,pV),AutoTsit5(Rosenbrock23()),reltol=1e-6);
		return xD
	end;

	# DY function
	# INPUT: Y    - Output of full system before perturbation
	#        YD   - Output of full system after perturbation
	#        Ynf  - Output of non-feedback system before perturbation (Ynf:=Y)
	#        YnfD - Output of non-feedback system after perturbation
	# OUPUT:      - DY value
	function DY(Y,YD,Ynf,YnfD)
		if abs(log10(YnfD/Ynf)) < 1e-4
			return NaN
		end
		return log10(YD/Y)/log10(YnfD/Ynf);
	end;

	# DY curve
	# INPUT: p     - Dictionary function with the ODE parameters & values
	#        pert  - Handle for the perturbation details
	#        motif - Handle for the considered motif
	#        uns  - 1 to use a slower, more stable ODE solver
	# OUPUT: DYs   - Vector of DY values for the range of parameters
	function DYc(p, pert, motif, uns)
		r = 10 .^ collect(pert.r[1]:pert.s:pert.r[2]);
		DYs = Array{Float64}(undef,length(r));
		DYs /= 0;
		p[pert.p] = pert.c;
		for i in 1:length(r)
			p[pert.p] *= r[i];
			rtol = 1e-6;
			flg1 = 1;
			ssR = ones(length(motif.odeFB.syms));
			soR = ones(length(motif.odeNF.syms));
			while(rtol >= 1e-24)
				# Reference steady state:
				ssR = SS(motif.odeFB, p, ssR, rtol, uns);
				# Locally analogous system reference steady state:
				motif.localNF(p,ssR);
				soR = SS(motif.odeNF, p, soR, rtol, uns);
				if(abs(motif.outFB(ssR) - motif.outNF(soR)) > 1e-4)
					rtol *= 1e-3;
					if(rtol < 1e-24)
						println("ERROR: Check NF system (reltol=",rtol*1e3,").")
						println(vcat(pert.p,i,[p[i] for i in motif.odeFB.params],motif.outFB(ssR),motif.outNF(soR)))
						#throw(DomainError("x-("))
						if(abs(motif.outFB(ssR) - motif.outNF(soR))/motif.outFB(ssR) > 0.01)
							flg1 = 0;
							println("SS results excluded!")
						end
					end
				else
					break
				end
			end
			# Perturbation:
			p[pert.p] *= pert.d;
			if(flg1==1)
				ssD = SS(motif.odeFB, p, ssR, rtol, uns);
				soD = SS(motif.odeNF, p, soR, rtol, uns);
				DYs[i] = DY(motif.outFB(ssR), motif.outFB(ssD), motif.outNF(soR), motif.outNF(soD));
			end
			p[pert.p] /= pert.d;
			p[pert.p] /= r[i];
		end
		return DYs
	end;

	# DY "metrics"
	# INPUT: DYs   - Vector of DY values for the range of parameters
	#        pert  - Handle for the perturbation details
	# OUPUT:       - [Range of DY<=eps,start,end,min(DY),optimal rho]
	function DYm(DYs, pert)
		r = 10 .^ collect(pert.r[1]:pert.s:pert.r[2]);
		i = DYs .<= pert.eps;
		j = findall(i);
		x = copy(DYs);
		x[x .=== NaN] .= Inf;
		if isempty(j)
			return [sum(DYs[i])./length(DYs[i]), NaN, NaN, minimum(x), r[argmin(x)]]
		end
		return [sum(DYs[i])./length(DYs[i]), pert.c * r[j[1]], pert.c * r[j[end]], minimum(x), r[argmin(x)]]
	end;

	# SSs for "optimal" control
	# INPUT: p     - Dictionary function with the ODE parameters & values
	#        pert  - Handle for the perturbation details
	#        motif - Handle for the considered motif
	#        DYs   - Vector of DY values for the range of parameters
	#        uns  - 1 to use a slower, more stable ODE solver
	# OUPUT:       - [Optimal rho, steady state for the full system before & after perturbation, and for the non-feedback system before & after perturbation]
	function SSopt(p, pert, motif, DYs, uns)
		r = 10 .^ collect(pert.r[1]:pert.s:pert.r[2]);
		x = copy(DYs);
		x[x .=== NaN] .= Inf;
		p[pert.p] = pert.c;
		p[pert.p] *= r[argmin(x)];
			rtol = 1e-6;
			ssR = zeros(length(motif.odeFB.syms));
			soR = zeros(length(motif.odeNF.syms));
			while(rtol >= 1e-24)
				# Reference steady state:
				ssR = SS(motif.odeFB, p, ssR, rtol, uns);
				# Locally analogous system reference steady state:
				motif.localNF(p,ssR);
				soR = SS(motif.odeNF, p, soR, rtol, uns);
				if(abs(motif.outFB(ssR) - motif.outNF(soR)) > 1e-4)
					rtol *= 1e-3;
					if(rtol < 1e-24)
						println("ERROR: Check NF system (reltol=",rtol*1e3,").")
						println(vcat(pert.p,i,[p[i] for i in motif.odeFB.params],motif.outFB(ssR),motif.outNF(soR)))
						#throw(DomainError("x-("))
					end
				else
					break
				end
			end
			# Perturbation:
			p[pert.p] *= pert.d;
			ssD = SS(motif.odeFB, p, ssR, rtol, uns);
			soD = SS(motif.odeNF, p, soR, rtol, uns);
		p[pert.p] /= pert.d;
		p[pert.p] /= r[argmin(x)];
		return [vcat(p[pert.p] * r[argmin(x)],ssR,ssD,soR,soD)]
	end;
end