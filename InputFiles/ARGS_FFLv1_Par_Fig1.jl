# Kinetic parameters
p = Dict([
    :g   => 0.01,      # Dilution rate (e.g. [0.01,0.24] 1/min)
    :mY  => 0.125,     # Y synthesis rate dependent of W (nM/min)
    :gY  => 0.1,       # Y degradation rate (e.g. [0.0301,0.0408] 1/min, with degODC)
    :mU  => 0.0334,    # U maximum synthesis rate repressed by Y (nM/min)
    :kD  => 1,         # Y repression KD (nM)
    :gU  => 0.01,      # U degradation rate (e.g. [0.0301,0.0408] 1/min, with degODC)
    :mW  => 0.125,     # W synthesis rate dependent of U (1/min)
    :gW  => 0.01,      # W degradation rate (e.g. [0.0301,0.0408] 1/min, with degODC)
    :mUs => NaN,       # LOCAL: U constitutive synthesis rate in the locally analogous system (mU*Yss)
]);

#Inital conditions
x0 = zeros(length(mm.odeFB.syms));
