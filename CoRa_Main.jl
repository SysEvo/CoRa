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
# iARG = (mm : Label for motif file, ex : Label for parameters file, 
#   pp : Label for perturbation type, an : Chose analysis type);
include(string("InputFiles/ARGS_",iARG.mm,"_Pert_",iARG.ex,".jl"))	# Perturbation details
include(string("InputFiles/ARGS_",iARG.mm,"_Par_",iARG.ex,".jl"))	# Core parameters
pO = copy(p);

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
elseif(iARG.an=="OptCoRa")
	include(string("InputFiles/ARGS_",iARG.mm,"_OptCoRa_",iARG.ex,".jl"))
	open(string("OUT_OptCoRa_",iARG.mm,"_",iARG.ex,"_",iARG.pp,"_",iARG.ax,".txt"), "w") do io
		if(mrw.prtD==1)
			writedlm(io, [vcat("Run","Iteration",[string(param) for param in mrw.pOp],string("|CoRa<=",pert.eps,"|"),"min(CoRa)",10 .^ collect(pert.r[1]:pert.s:pert.r[2]))], '\t')
		else
			writedlm(io, [vcat("Run","Iteration",[string(param) for param in mrw.pOp],string("|CoRa<=",pert.eps,"|"),"min(CoRa)")], '\t')
		end
		for ruN in 1:mrw.runs
			println("RUN #",ruN)
			p = copy(pO);
			# If random initial parameters specified:
			if(mrw.ranP==1)
				fn.mrwR(mrw,p)
			end
			# Initialize system:
			CoRa0,op0,mi0 = fn.CoRam(p,pert,mm,x0); # Properties to optimize, proportion of CoRas<=eps and min(CoRas) value)
			# Print initial state:
			if(mrw.prtD==1)	# If printing full CoRa curve specified:
				writedlm(io, [vcat(ruN,0,[p[i] for i in mrw.pOp],op0,mi0,CoRa0)],'\t')
			else			# Else:
				writedlm(io, [vcat(ruN,0,[p[i] for i in mrw.pOp],op0,mi0)],'\t')
			end
			# Optimization iterations
			println("I: minCoRa = ",mi0,"\t |CoRa<eps| = ",op0)
			for i in 1:mrw.iter
				CoRa0,op0,mi0 = fn.mrwI(mrw,p,pert,mm,x0,CoRa0,op0,mi0);
				# Print if "printing each walk step" or if last iteration:
				if(mrw.prtW==1 || i==mrw.iter)
					if(mrw.prtD==1)
						writedlm(io, [vcat(ruN,i,[p[i] for i in mrw.pOp],op0,mi0,CoRa0)],'\t')
					else
						writedlm(io, [vcat(ruN,i,[p[i] for i in mrw.pOp],op0,mi0)],'\t')
					end
				end
				# End if optimal CoRa in all conditions and print final values:
				if(op0==1.0)
					println("Optimal value (|CoRa<=eps|=",op1,") reached at iteration ",i)
					if(mrw.prtD==1)
						writedlm(io, [vcat(ruN,i,[p[i] for i in mrw.pOp],op0,mi0,CoRa0)],'\t')
					else
						writedlm(io, [vcat(ruN,i,[p[i] for i in mrw.pOp],op0,mi0)],'\t')
					end
					break;
				end
			end
			println("F: minCoRa = ",mi0,"\t |CoRa| = ",op0,"\n")
		end
	end
else
	println("ERROR: Undetermined analysis. Options: ExSSs, ExDyn, CoRams, OptCoRa")
end
