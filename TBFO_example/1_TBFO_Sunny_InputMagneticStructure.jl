using Sunny, Optim, StaticArrays, DataFrames, CSV, LinearAlgebra, Plots, GLM, Statistics, StatsBase, Printf, GLMakie

if Sys.iswindows()
    SunnyAssist_path = "C:\\Users\\vdp\\Dropbox (ORNL)\\Sunny\\SunnyAssist\\"
else
    SunnyAssist_path = "/home/vdp/Dropbox (ORNL)/Sunny/SunnyAssist/"
end

include(joinpath(SunnyAssist_path, "SunnyAssist.jl"))

xtal    = Crystal("C:\\Users\\vdp\\Dropbox (ORNL)\\Sunny\\SunnyAssist\\TBFO_example\\TBFO_CollCode246645.cif"; symprec=1e-4)
magxtal = subcrystal(xtal,"Fe1", "Fe2")

view_crystal(magxtal, 5.0)

# Assign local Hilbert spaces
lhs_Fe3 = SpinInfo(1, S=2, g=2)
formfactor_Fe3 = [FormFactor("Fe3")];
lhs_Fe2 = SpinInfo(3, S=2.5, g=2)
formfactor_Fe2 = [FormFactor("Fe2")];


# Create `System`
sunmode = :dipole
latsize = (1,2,1)
sys     = System(magxtal, latsize, [lhs_Fe3, lhs_Fe2], sunmode; seed=1)

original_spins = Array{SVector{3, Float64}, 4}(undef, 1, 2, 1, 4)

original_spins[1, 1, 1, 1] = SVector(0.0, 2.0, 0.0)
original_spins[1, 1, 1, 2] = SVector(0.0, 2.0, 0.0)
original_spins[1, 1, 1, 3] = SVector(0.0, -2.5, 0.0)
original_spins[1, 1, 1, 4] = SVector(0.0, -2.5, 0.0)

original_spins[1, 2, 1, 1] = SVector(0.0, -2.0, 0.0)
original_spins[1, 2, 1, 2] = SVector(0.0, -2.0, 0.0)
original_spins[1, 2, 1, 3] = SVector(0.0, 2.5, 0.0)
original_spins[1, 2, 1, 4] = SVector(0.0, 2.5, 0.0)

set_initial_spin_configuration!(sys, original_spins)
