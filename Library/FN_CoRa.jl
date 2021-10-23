# Steady state & DY calculation functions
#		Julia v.1.6

module fn
	# Required libraries
	using DelimitedFiles
	using Distributions
	using DifferentialEquations
	using DelimitedFiles
	using Distributions

	###This function is rather simple, but it's for readability in the main code purposes, and prevention of typos
	function Restart(a, b)
		a = zeros(length(a)).+NaN;
		b = zeros(length(b)).+NaN;
		return a, b
	end;

	###Can't have a steady state & check function without the steady state and the checking, now can we?
	function SS(syst, p, x0, rtol)
		pV = [p[eval(Meta.parse(string(":",i)))] for i in syst.sys.ps];
		tS = 0;
		dXrm = 1;
		while(dXrm > rtol)
			ss = try
				solve(ODEProblem(syst,x0,1e6,pV); reltol=rtol,save_everystep = false);
			catch
				try
					solve(ODEProblem(syst,x0,1e6,pV),alg_hints=[:stiff]; reltol=rtol,save_everystep = false);
				catch err
					println("WARNING: Error in ODE simulation: <<",err,">>. ss --> NaN")
					x0 = zeros(length(syst.syms)).+NaN;
					break
				end
			end;
			dXrm = maximum(abs.(big.(ss(1e6))-big.(ss(1e6-0.01)))./big.(ss(1e6)));
			x0 = ss(1e6);
			tS += 1e6;
			if(tS>=1e12)
				println("WARNING: Maximum iteration reached (simulated time 1e18). Max relative Delta: ",dXrm)
				break
			end
		end
		return x0
	end;

	###A small function to check for NaNs, again made for the purposes of readability within the code
	function NaNCheck(a, b)
		if(any(isnan.(a)))
			a, b = fn.Restart(a, b)
			return a, b
		else
			return "Valid"
		end
	end;

	###Now we create the Check function!
	function Check(ssR, soR, rtol, mm)
		if(abs(mm.outFB(ssR) - mm.outNF(soR)) > 1e-4)
			rtol *= 1e-3;
			if(rtol < 1e-24)
				println("ERROR: Check NF system (reltol=",rtol,").")
				println(vcat(pert.p,i,[p[eval(Meta.parse(string(":",i)))] for i in syst.sys.ps],mm.outFB(ssR),mm.outNF(soR)))
				if(abs(mm.outFB(ssR) - mm.outNF(soR))/mm.outFB(ssR) > 0.01)
					ssR = zeros(length(mm.odeFB.syms)).+NaN;
					soR = zeros(length(mm.odeNF.syms)).+NaN;
					println("Error too large. SS results excluded!")
				end
			end
			return ssR, soR, rtol, "Insufficient"
		else
			return ssR, soR, rtol, "Sufficient"
		end
	end;

	###Here it is, the SS&Check function. It will, of course, be built upon an SS() and a Check() funciton previously defined within this document
	function SSandCheck(ssR, soR, p, x0, rtol, mm)
		###The rtol value basically states that this will be attempted 5 times: This is because of the comparison made in "if(abs(mm.outFB(ssR) - mm.outNF(soR)) > 1e-4)". If that is successful
		###Then the value of rtol is multiplied by 1e-3, until it is no longer greater or equal than 1e-24, as stated in our while() condition
		flag = "Insufficient"
		while(rtol >= 1e-24 && flag == "Insufficient")
			# Reference steady state:
			ssR = fn.SS(mm.odeFB, p, x0, rtol);
			if(NaNCheck(ssR, soR) != "Valid")
				println("Condition excluded! ssR --> NaN");
				break;
			end
			# Locally analogous system reference steady state:
			mm.localNF(p,ssR);
			###Of note here, instead of using the initial condition x0, we use ssR.
			soR = fn.SS(mm.odeNF, p, ssR, rtol);
			if(NaNCheck(soR, ssR) != "Valid")
				println("Condition excluded! soR --> NaN");
				break;
			end
			ssR, soR, rtol, flag = Check(ssR, soR, rtol, mm)
		end
		return ssR, soR, rtol
	end;

	function Perturbation(ssR, soR, p, rtol, mm)
		ssD = fn.SS(mm.odeFB, p, ssR, rtol);
		soD = fn.SS(mm.odeNF, p, soR, rtol);
		return ssD, soD
	end;

	function CoRa(ssR, ssD, soR, soD)
		if abs(log10(soD/soR)) < 1e-4
			return NaN
		end
		return log10(ssD/ssR)/log10(soD/soR);
	end;

end