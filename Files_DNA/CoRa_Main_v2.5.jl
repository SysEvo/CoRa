###CoRa_v2.5

###Alright, once and for all
###In the spirit of proper documentation, I'll be as explicit as I can with these comments
###Taken from v2, no need to modify this just yet
## Load functions & parameters:
using DelimitedFiles
using Distributions

#This includes the model to study itself, it must be prepared before running CoRa, with the explicit constraint that the number of d/dt are the same in the FB and NF equations
mm = include(string("Library\\Md_",iARG.mm,".jl"));
###This next line is changed to use the "updated" functions .jl
fn = include(string("Library\\FN_CoRa.jl"));
## INPUTS:
# iARG = (mm : Label for motif file, ex : Label for parameters file, pp : Label for perturbation type, an : Chose analysis type);
include(string("InputFiles\\ARGS_",iARG.mm,"_Pert_",iARG.ex,".jl"))	# Perturbation details
include(string("InputFiles\\ARGS_",iARG.mm,"_Par_",iARG.ex,".jl"))	# Core parameters
p0 = copy(p);

if(iARG.an == "ExSSs")
    ###The tag was replaced from "io" to "outfile1" with the intention of creating a "report card" output file down the line, which would necessitate the existence of multiple output files
    open(string("OUT_ExSSs_",iARG.mm,"_",iARG.ex,"_",iARG.pp,"_",iARG.ax,".txt"), "w") do outfile1
        ###This generates the headers for our data output, we're preparing the file beforehand
        writedlm(outfile1, [vcat(iARG.ax,[string("FbR_",i) for i in mm.odeFB.syms],[string("FbD_",i) for i in mm.odeFB.syms],[string("NfR_",i) for i in mm.odeNF.syms],[string("NfD_",i) for i in mm.odeNF.syms],string("CoRa(",iARG.pp,")"))],'\t');
        ###This creates our "steps" for the data production, as stated in the .*_Pert_.*.jl file
        r = 10 .^ collect(pert.r[1]:pert.s:pert.r[2]);
        ###So, then, this will repeat the for loop in however many steps for the parameters
        for i in 1:length(r)
            ###The parameter to change is multiplied by the corresponding value in our steps collection, and the error tolerance is set to 1e-12 (arbitrarily)
            p[pert.c] *= r[i];
            ###Up next, we must find the steady states of both ssR and soR, while also checking that the process itself didn't fail. Let's make that a single, "fn.SSandCheck()" function
            ssR, soR, rtol = fn.SSandCheck(p, x0, 1e-12, mm)
            ###Now, we have to do the Perturbation itself!
            ssD, soD = fn.Perturbation(ssR, soR, p, rtol, mm, pert)
            ###And we must return CoRa now, too
            CoRa = fn.CoRa(mm.outFB(ssR), mm.outFB(ssD), mm.outNF(soR), mm.outNF(soD))
            ###With everything done, it's time to output them into the file!
            writedlm(outfile1, [vcat(p[pert.c],ssR,ssD,soR,soD,CoRa)],'\t');
            ###And we reset our perturbation for the next run of the loop
            p[pert.c] /= r[i];
        end
    end
elseif(iARG.an == "CoRams")
    include(string("InputFiles\\ARGS_",iARG.mm,"_DYms_",iARG.ex,".jl"))	# Parameters to vary
	open(string("OUT_DYms_",iARG.mm,"_",iARG.ex,"_",iARG.pp,"_",iARG.ax,".txt"), "w") do outfile1
		writedlm(outfile1, [vcat([string(i) for i in keys(pN)],10 .^ collect(pert.r[1]:pert.s:pert.r[2]))],'\t')
		for pI = pN
			for i = pI[2]
				p = copy(p0);
				p[pI[1]] *= (10. ^i);
				writedlm(outfile1, [vcat([p[j[1]] for j in pN], fn.CoRac(p,pert,mm, x0))],'\t')
			end
		end
	end
end