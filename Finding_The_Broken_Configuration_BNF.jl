cd("C:\\Users\\ese_1\\OneDrive\\Documentos\\CoRa")
using Pkg;
using CSV;
using DelimitedFiles;
using Distributions;

Pkg.activate(".");		# Activate local environment (requiere '.toml' files)
iARG = (mm = "BNFv1",  # Label for motif file
ex = "1250Set1",      # Label for parameters file
pp = :mY,         # Label for perturbation type
ax = :mY);    # Label for condition/environment
a = readdlm("OUT_ExSSs_BNFv1_1250Set1_mY_mY_From_498.txt")
last_index = floor(Int, a[size(a,1), 1]) 

print(string(last_index, "\n"))

pars = CSV.File("InputFiles\\ARGS_BNFv1_Mass_Par_1250Set1.csv") # Core parameters
mm = include(string("Library\\Md_",iARG.mm,".jl"));
fn = include(string("Library\\FN_CoRa.jl"));
include(string("InputFiles\\ARGS_",iARG.mm,"_Pert_",iARG.ex,".jl")) # Perturbation details
key_names = (:g, :mY, :gY, :mU, :kD, :gU, :b, :bs, :mUs);
x0 = zeros(length(mm.odeFB.syms));

open(string("OUT_ExSSs_",iARG.mm,"_",iARG.ex,"_",iARG.pp,"_",iARG.ax,"_LastLine.txt"), "w") do outfile1
    r = 10 .^ collect(pert.r[1]:pert.s:pert.r[2]);
    for i in last_index:pars.rows
        p = Dict();
        for h in 1:pars.cols
            push!(p, key_names[h] => pars[i][h])
        end
        p0 = copy(p)
        ###The tag was replaced from "io" to "outfile1" with the intention of creating a "report card" output file down the line, which would necessitate the existence of multiple output files
        ###This generates the headers for our data output, we're preparing the file beforehand
        ###This creates our "steps" for the data production, as stated in the .*_Pert_.*.jl file
        CoRa = zeros(length(r)).+ NaN;
        ###So, then, this will repeat the for loop in however many steps for the parameters
        k = 1;
        try
            for k in 1:length(r)
                ###The parameter to change is multiplied by the corresponding value in our steps collection, and the error tolerance is set to 1e-12 (arbitrarily)
                p[pert.c] *= r[k];
                ###Up next, we must find the steady states of both ssR and soR, while also checking that the process itself didn't fail. Let's make that a single, "fn.SSandCheck()" function
                ssR, soR, rtol = fn.SSandCheck(p, x0, 1e-12, mm)
                ###Now, we have to do the Perturbation itself!
                ssD, soD = fn.Perturbation(ssR, soR, p, rtol, mm, pert)
                ###And we must return CoRa now, too
                CoRa[k] = fn.CoRa(mm.outFB(ssR), mm.outFB(ssD), mm.outNF(soR), mm.outNF(soD))
                ###With everything done, it's time to output them into the file!
                ###Check this one out
                ###And we reset our perturbation for the next run of the loop
                p[pert.c] = p0[pert.c];
            end
        catch
        end
        writedlm(string("OUT_ExSSs_",iARG.mm,"_",iARG.ex,"_",iARG.pp,"_",iARG.ax,"_LastLine.txt"), [vcat(i, CoRa)],'\t');
        print(string("Line ", i, " done! \n"))
    end
end