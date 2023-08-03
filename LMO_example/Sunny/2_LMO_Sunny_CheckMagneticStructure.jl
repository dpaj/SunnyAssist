using Sunny, Optim, StaticArrays, DataFrames, CSV, LinearAlgebra, Plots, GLM, Statistics, StatsBase, Printf, GLMakie

if Sys.iswindows()
    SunnyAssist_path = "C:\\Users\\vdp\\Dropbox (ORNL)\\Sunny\\"
else
    SunnyAssist_path = "/home/vdp/Dropbox (ORNL)/Sunny/"
end

include(joinpath(SunnyAssist_path, "SunnyAssist.jl"))


a = 5.5333
b = 5.7461
c = 7.6637

latvecs = lattice_vectors(a, b, c, 90, 90, 90)   # A 3x3 matrix of lattice vectors that
                                                 ## define the conventional unit cell

positions = [[-0.0095, 0.0513, 1/4], [1/2, 0, 0] , [0.0777, 0.48493, 1/4], [0.7227, 0.3085, 0.0408]]  # Positions of atoms in fractions
                                                         ## of lattice vectors

types = ["La", "Mn", "O", "O"]

spacegroup = "Pbnm"

LMO = Crystal(latvecs, positions, spacegroup; types)

cryst = subcrystal(LMO, "Mn")

sys = System(cryst, (1,1,1), [SpinInfo(1,S=2, g=2)], :dipole)

original_spins = [[0, -1, 0], [0, -1, 0], [0, 1, 0], [0, 1, 0]]
set_initial_spin_configuration!(sys, original_spins)

GLMakie.activate!(inline=false)

plot_spins(sys; arrowlength=2.5, linewidth=0.75, arrowsize=1.5)