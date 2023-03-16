## CoRa ANALYSIS
#	Mariana Gómez-Schiavon
#	March, 2023
#		Julia v.1.8
#		Required libraries:
#			DifferentialEquations
#			ParameterizedFunctions
#			Statistics
#			Distributions
#			DelimitedFiles

# Load functions & parameters:
using DelimitedFiles
using Distributions
## INPUTS:
mm = include(string("Library/Md_",iARG.mm,".jl"));
fn = include(string("Library/FN_CoRa.jl"));
pO = copy(p);
# iARG = (mm : Label for motif file, ex : Label for parameters file, 
#   pp : Label for perturbation type, an : Chose analysis type);
include(string("InputFiles/ARGS_",iARG.mm,"_Pert_",iARG.ex,".jl"))	# Perturbation details
include(string("InputFiles/ARGS_",iARG.mm,"_Par_",iARG.ex,".jl"))	# Core parameters

## Run analysis
# Calculate CoRa curve for a range of parameters:
if(iARG.an=="ExSSs")
	p = copy(pO);
	open(string("OUT_ExSSs_",iARG.mm,"_",iARG.ex,"_",iARG.pp,"_",iARG.ax,".txt"), "w") do io
		# Generate headers for the data output:
        writedlm(io, [vcat(iARG.ax,[string("FbR_",i) for i in mm.odeFB.syms],[string("FbD_",i) for i in mm.odeFB.syms],[string("NfR_",i) for i in mm.odeNF.syms],[string("NfD_",i) for i in mm.odeNF.syms],string("CoRa(",iARG.pp,")"))],'\t');
		# Define the range of conditions to evaluate, as stated in 
		#   the .*_Pert_.*.jl file:
        r = 10 .^ collect(pert.r[1]:pert.s:pert.r[2]);
        # For each condition:
		for i in 1:length(r)
			# Update condition by multiplying the choosen parameter 
			#   by the corresponding value in our range of conditions:
			p[pert.c] *= r[i];
            # Find the steady state of the reference system (ssR), 
			#   update the analogous system accordingly and find 
			#   its steady state (soR), and finally confirm that 
			#   both systems are locally analogous (ssR=soR):
            ssR, soR, rtol = fn.SSandCheck(p, x0, 1e-12, mm)
            # Perform the perturbation and find thesteady states 
			#   for both systems:
            ssD, soD = fn.Perturbation(ssR, soR, p, rtol, mm, pert)
            # Calculate the resulting CoRa metric:
            CoRa = fn.CoRa(mm.outFB(ssR), mm.outFB(ssD), mm.outNF(soR), mm.outNF(soD))
            # Print output:
            writedlm(io, [vcat(p[pert.c],ssR,ssD,soR,soD,CoRa)],'\t');
			# Reset condition to pre-perturbation state:
            p[pert.c] /= r[i];
        end
	end
# Calculate dynamic response after a perturbation:
elseif(iARG.an=="ExDyn")
	p = copy(pO);
	open(string("OUT_ExDyn_",iARG.mm,"_",iARG.ex,"_",iARG.pp,"_",iARG.ax,".txt"), "w") do io
		writedlm(io, [vcat("FB","rho","time",[string(i) for i in mm.odeNF.syms])],'\t');
		# Reference steady state:
		ssR, soR, rtol = fn.SSandCheck(p, x0, 1e-12, mm)
		# Feedback system:
		x = fn.Dyn(mm.odeFB, p, ssR, 500.0);
		for i in 1:length(x.t)
			writedlm(io, [vcat(1,p[iARG.pp],x.t[i],x.u[i],"NaN")],'\t');
		end
		p[pert.p] *= pert.d;
		x = fn.Dyn(mm.odeFB, p, last(x.u), 9500.0);
		for i in 1:length(x.t)
			writedlm(io, [vcat(1,p[iARG.pp],x.t[i]+500.0,x.u[i],"NaN")],'\t');
		end
		ssD = fn.SS(mm.odeFB, p, ssR, rtol);
		writedlm(io, [vcat(1,p[iARG.pp],"Inf",ssD,"NaN")],'\t');
		p[pert.p] /= pert.d;
		# No-Feedback system:
		x = fn.Dyn(mm.odeNF, p, soR, 500.0);
		for i in 1:length(x.t)
			writedlm(io, [vcat(0,p[iARG.pp],x.t[i],x.u[i])],'\t');
		end
		p[pert.p] *= pert.d;
		x = fn.Dyn(mm.odeNF, p, last(x.u), 9500.0);
		for i in 1:length(x.t)
			writedlm(io, [vcat(0,p[iARG.pp],x.t[i]+500.0,x.u[i])],'\t');
		end
		soD = fn.SS(mm.odeNF, p, soR, rtol);
		writedlm(io, [vcat(0,p[iARG.pp],"Inf",soD)],'\t');
		p[pert.p] /= pert.d;
	end
# Calculate CoRa curve for a range of parameters as another parameter varies:
elseif(iARG.an=="CoRams")
	include(string("InputFiles/ARGS_",iARG.mm,"_CoRams_",iARG.ex,".jl"))	# Parameters to vary
	open(string("OUT_CoRams_",iARG.mm,"_",iARG.ex,"_",iARG.pp,"_",iARG.ax,".txt"), "w") do io
		writedlm(io, [vcat([string(i) for i in keys(pN)],10 .^ collect(pert.r[1]:pert.s:pert.r[2]))],'\t')
		for pI = pN
			for i = pI[2]
				p = copy(pO);
				p[pI[1]] *= (10. ^i);
				CoRas = fn.CoRac(p, pert, mm, x0);
				writedlm(io, [vcat([p[j[1]] for j in pN], CoRas)],'\t')
				p[pI[1]] /= (10. ^i);
			end
		end
	end
# Optimize CoRa curve for a range of parameters:
elseif(iARG.an=="OptDY")
	include(string("InputFiles/ARGS_",iARG.mm,"_OptDY_",iARG.ex,".jl"))
	open(string("OUT_OptDY_",iARG.mm,"_",iARG.ex,"_",iARG.pp,"_",iARG.ax,".txt"), "w") do io
		if(mrw.prtD==1)
			writedlm(io, [vcat("Run","Iteration",[string(param) for param in mrw.pOp],string("|DY<=",pert.eps,"|"),"min(DY)",10 .^ collect(pert.r[1]:pert.s:pert.r[2]))], '\t')
		else
			writedlm(io, [vcat("Run","Iteration",[string(param) for param in mrw.pOp],string("|DY<=",pert.eps,"|"),"min(DY)")], '\t')
		end
		for ruN in 1:mrw.runs
			println("RUN #",ruN)
			uns = 0;
			p = copy(pO);
			####### Uncomment the next three lines for random initial conditions: #######
			#for i in 1:length(mrw.pOp)
			#	p[mrw.pOp[i]] = 10 .^ (rand(Uniform(mrw.pMin[i], mrw.pMax[i])));
			#end
			#############################################################################
			## Temperature function for simulated annealing:
			if(mrw.temp==1)
				mrwT = collect(mrw.iter:-1:1) ./ mrw.iter;
			else
				mrwT = ones(mrw.iter); # NOTE: For MRW, make T=1.
			end
			## Initialize system
			DY0 = fn.DYc(p,pert,mm,uns);    # Calculate DY curve
			DYm = fn.DYm(DY0, pert);    # Calculate metrics of DY curve
			op0 = log10(DYm[3]/DYm[2]);   # Property to optimize (e.g. DY<=eps range length)
			mi0 = DYm[4];                 # Secondary property to optimize (e.g. min(DY) value)
			r0 = zeros(length(mrw.pOp));
			if(mrw.prtD==1)
				writedlm(io, [vcat(ruN,0,[p[i] for i in mrw.pOp],op0,mi0,DY0)],'\t')
			else
				writedlm(io, [vcat(ruN,0,[p[i] for i in mrw.pOp],op0,mi0)],'\t')
			end
			# Optimization iterations
			println("I: minDY = ",mi0,"\t |DY| = ",op0)
			for i in 1:mrw.iter
				rI = rand(MvNormal(zeros(length(mrw.pOp)), zeros(length(mrw.pOp)) .+ mrw.cov)); # Random values to update parameters
				for pI in 1:length(mrw.pOp)  # Update parameter values
					r0[pI] = p[mrw.pOp[pI]]; # Save previous value
					p[mrw.pOp[pI]] *= (mrw.M .^ rI[pI]); # Update value
					# Exclude values outside regime of exploration:
					if p[mrw.pOp[pI]] < (10.0 ^ mrw.pMin[pI])
						p[mrw.pOp[pI]] = (10.0 ^ mrw.pMin[pI])
					elseif p[mrw.pOp[pI]] > (10.0 ^ mrw.pMax[pI])
						p[mrw.pOp[pI]] = (10.0 ^ mrw.pMax[pI])
					end
				end
				DYs = fn.DYc(p,pert,mm,uns);  # Calculate new DY curve
				DYm = fn.DYm(DYs, pert);  # Calculate new metrics of DY curve
				op1 = log10(DYm[3]/DYm[2]); # New value of property to optimize (e.g. DY<=eps range length)
				mi1 = DYm[4];               # New value of secondary property to optimize (e.g. min(DY) value)
				# Evaluate if accept new parameter values or not:
				## Only accept in the regime of interest, i.e. DY>=0:
				c1 = (mi1>=0);
				## If DY>eps for all rho, evaluate the min(DY) for both sets:
				### NOTE: As mi0,mi1=[0,1], correct exponential with the expected variance of ~U(0,1)
					xiC = (mi0 ^ 2) / (2 * 0.083);
					xiP = (mi1 ^ 2) / (2 * 0.083);
				c2 = isnan(op0+op1) && (rand() < exp((xiC - xiP) / mrwT[i]));
				## If DY>=eps for some rho, evaluate the |DY<=eps| for both sets:
				### NOTE: As op0,op1=[0,rrO], but still correct exponential with the expected variance of ~U(0,1)
				###       !! ~U(0,1)*(rrO^2) variance resulted in very noisy runs...
					rrO = pert.r[2] - pert.r[1];
					xiC = (rrO - op0) / (2 * 0.083);
					xiP = (rrO - op1) / (2 * 0.083);
				c3 = rand() < exp((xiC - xiP) / mrwT[i]);
				if(c1 && (c2 || c3))
					# If yes, update "reference" system
					op0 = op1;
					mi0 = mi1;
					DY0 = DYs;
				else
					# If not, revert to previous parameter values
					for pI in 1:length(mrw.pOp)
						p[mrw.pOp[pI]] = r0[pI];
					end
				end
				if(mrw.prtW==1 || i==mrw.iter)
					if(mrw.prtD==1)
						writedlm(io, [vcat(ruN,i,[p[i] for i in mrw.pOp],op0,mi0,DY0)],'\t')
					else
						writedlm(io, [vcat(ruN,i,[p[i] for i in mrw.pOp],op0,mi0)],'\t')
					end
				end
				if(op1==rrO)
					println("Optimal value (|DY<=eps|=",op1,") reached at iteration ",i)
					if(mrw.prtD==1)
						writedlm(io, [vcat(ruN,i,[p[i] for i in mrw.pOp],op0,mi0,DY0)],'\t')
					else
						writedlm(io, [vcat(ruN,i,[p[i] for i in mrw.pOp],op0,mi0)],'\t')
					end
					break;
				end
			end
			println("F: minDY = ",mi0,"\t |DY| = ",op0,"\n")
		end
	end
else
	println("ERROR: Undetermined analysis. Options: ExSSs, ExDyn, DYms, OptDY")
end
