cd("C:\\Users\\ese_1\\OneDrive\\Documentos\\CoRa")
using Pkg;
using CSV;
using DelimitedFiles;
using Distributions;

a = readdlm("OUT_ExSSs_ATFv2_1250Set2_mY_mY.txt")
start = floor(Int, a[size(a,1), 1]) + 1

Pkg.activate(".\\CoRa");		# Activate local environment (requiere '.toml' files)
iARG = (mm = "ATFv2",  # Label for motif file
ex = "1250Set2",      # Label for parameters file
pp = :mY,         # Label for perturbation type
ax = :mY);    # Label for condition/environment
pars = CSV.File("InputFiles\\ARGS_ATFv2_Mass_Par_1250Set2.csv") # Core parameters
strict = true
gap_tol = 10 # This is how wide we accept gaps to be in lines
# print = "short"      # Flag for the output file including/excluding the steady states (Options: "short", "all")

#This includes the model to study itself, it must be prepared before running CoRa, with the explicit constraint that the number of d/dt are the same in the FB and NF equations
mm = include(string("Library\\Md_",iARG.mm,".jl"));
###This next line is changed to use the "updated" functions .jl
fn = include(string("Library\\FN_CoRa_v2.jl"));
## INPUTS:
# iARG = (mm : Label for motif file, ex : Label for parameters file, pp : Label for perturbation type, an : Chose analysis type);
include(string("InputFiles\\ARGS_",iARG.mm,"_Pert_",iARG.ex,".jl")) # Perturbation details
key_names = (:g, :mY, :gY, :mU, :gU, :mW, :gW, :e0, :eP, :eM, :mUs);
x0 = zeros(length(mm.odeFB.syms));
open(string("OUT_ExSSs_",iARG.mm,"_",iARG.ex,"_",iARG.pp,"_",iARG.ax, "_From_", start, ".txt"), "w") do outfile1
    r = 10 .^ collect(pert.r[1]:pert.s:pert.r[2]);
    r = round.(r, sigdigits = 10)
    writedlm(outfile1, [vcat(string("Row"), r)],'\t');
    for i in start:pars.rows
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
            for k in 1:lastindex(r)
                if ((k > 2*gap_tol) && ((sum(isnan.(CoRa[(k-2*gap_tol):k])))/(2*gap_tol)) >= 0.5)
                    println("Set", i, " has failed in", gap_tol, " out of ", gap_tol * 2, "continuous points, so it will no longer be computed")
                    CoRa = zeros(length(r)).+ NaN;
                    break
                end
                ###The parameter to change is multiplied by the corresponding value in our steps collection, and the error tolerance is set to 1e-10 (arbitrarily)
                p[pert.c] *= r[k];
                ###Up next, we must find the steady states of both ssR and soR, while also checking that the process itself didn't fail. Let's make that a single, "fn.SSandCheck()" function
                ssR, soR, rtol = fn.SSandCheck(p, x0, 1e-10, mm, strict, k)
                if any(isnan.(ssR)) || any(isnan.(soR))
                    continue
                end
                ###Now, we have to do the Perturbation itself!
                ssD, soD = fn.Perturbation(ssR, soR, p, rtol, mm, pert, strict, k)
                if any(isnan.(ssD)) || any(isnan.(soD))
                    continue
                end
                ###And we must return CoRa now, too
                CoRa[k] = fn.CoRa(mm.outFB(ssR), mm.outFB(ssD), mm.outNF(soR), mm.outNF(soD))
                ###With everything done, it's time to output them into the file!
                ###Check this one out
                ###And we reset our perturbation for the next run of the loop
                p[pert.c] = p0[pert.c];
            end
        catch
        end
        writedlm(outfile1, [vcat(i, CoRa)],'\t');
        print(string("Line ", i, " done! \n"))
    end
end