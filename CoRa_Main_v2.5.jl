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
pO = copy(p);

if(iARG.an == "ExSSs")
    ###The tag was replaced from "io" to "outfile1" with the intention of creating a "report card" output file down the line, which would necessitate the existence of multiple output files
    open(string("OUT_ExSSs_",iARG.mm,"_",iARG.ex,"_",iARG.pp,"_",iARG.ax,".txt"), "w") do outfile1
        ###This generates the headers for our data output, we're preparing the file beforehand
        writedlm(outfile1, [vcat(iARG.ax,[string("FbR_",i) for i in mm.odeFB.syms],[string("FbD_",i) for i in mm.odeFB.syms],[string("NfR_",i) for i in mm.odeNF.syms],[string("NfD_",i) for i in mm.odeNF.syms],string("CoRa(",iARG.pp,")"))],'\t');
        ###This creates our "steps" for the data production, as stated in the .*_Pert_.*.jl file
        r = 10 .^ collect(pert.r[1]:pert.s:pert.r[2]);
        ###Before we actually go into the for loop, we need to define ssR and soR so within the for loop, we can iterate a cleaner "restart" function on them instead
        ssR = zeros(length(mm.odeFB.syms)).+NaN;
        soR = zeros(length(mm.odeNF.syms)).+NaN;
        ###So, then, this will repeat the for loop in however many steps for the parameters
        for i in 1:length(r)
            ###The parameter to change is multiplied by the corresponding value in our steps collection, and the error tolerance is set to 1e-12 (arbitrarily)
            p[pert.c] *= r[i];
            rtol = 1e-12;
            ###Now, ssR & soR must be reset into a zeroes vector of the length of their respective model. We'll do this with the fn.Eestart() function
            ssR, soR = fn.Restart(ssR, soR)
            ###Up next, we must find the steady states of both ssR and soR, while also checking that the process itself didn't fail. Let's make that a single, "fn.SSandCheck()" function
            ssR, soR, rtol = fn.SSandCheck(ssR, soR, p, x0, rtol, mm)
            ###Now, we have to do the Perturbation itself!
            p[pert.p] *= pert.d;
            ssD, soD = fn.Perturbation(ssR, soR, p, rtol, mm)
            ###And we must return CoRa now, too
            CoRa = fn.CoRa(mm.outFB(ssR), mm.outFB(ssD), mm.outNF(soR), mm.outNF(soD))
            ###With everything done, it's time to output them into the file!
            writedlm(outfile1, [vcat(p[pert.c],ssR,ssD,soR,soD,CoRa)],'\t');
            ###And we reset our perturbation for the next run of the loop
            p[pert.p] /= pert.d;
            p[pert.c] /= r[i];
        end
    end
end