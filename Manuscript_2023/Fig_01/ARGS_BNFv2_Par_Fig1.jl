# Kinetic parameters
p = Dict([
    :g   => 0.01,      # Dilution rate (e.g. [0.01,0.24] 1/min)
    :mY  => 1.0,       # Y synthesis rate dependent of W (nM/min)
    :gY  => 0.1,       # Y degradation rate (e.g. [0.0301,0.0408] 1/min, with degODC)
    :mU  => 1.0,       # U maximum synthesis rate repressed by Y (nM/min)
    :kD  => 1.0,       # Y repression KD (nM)
    :gU  => 0.0001,    # U degradation rate (e.g. [0.0301,0.0408] 1/min, with degODC)
    :b   => 0.0108,    # U to Us transition rate (1/min)
    :bs  => 0.068,    # Us to U transition rate (1/min)
    :mUs => NaN,       # LOCAL: U constitutive synthesis rate in the locally analogous system (mU * (kD/(Yss + kD)))
]);

#Inital conditions
x0 = zeros(length(mm.odeFB.syms));
