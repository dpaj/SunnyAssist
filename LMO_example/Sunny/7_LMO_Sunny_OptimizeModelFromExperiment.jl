using Sunny, Optim, StaticArrays, DataFrames, CSV, LinearAlgebra, Plots, GLM, Statistics, StatsBase, Printf, GLMakie

if Sys.iswindows()
    SunnyAssist_path = "C:\\Users\\vdp\\Dropbox (ORNL)\\Sunny\\"
else
    SunnyAssist_path = "/home/vdp/Dropbox (ORNL)/Sunny/"
end

include(joinpath(SunnyAssist_path, "SunnyAssist.jl"))


# Here is the data path
if Sys.iswindows()
    data_path = "C:\\Users\\vdp\\Dropbox (ORNL)\\Sunny\\"
else
    data_path = "/home/vdp/Dropbox (ORNL)/Sunny/"
end

# Here are your data files
data_files = ["LaMnO3_dispersion_a.csv", "LaMnO3_dispersion_b.csv", "LaMnO3_dispersion_c.csv"]

# Prepend the data path to each file
data_files = [joinpath(data_path, file) for file in data_files]

LMO_data = load_group(data_files)

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

function set_system_params!(sys::System, parameters::Dict{String, Float64})
    S = spin_operators(sys, 1)

    set_exchange!(sys, [parameters["Jab"]*2.0   0.0   0.0;
                        0.0   parameters["Jab"]*2.0   0.0;
                        0.0   0.0   parameters["Jab"]*2.0], Bond(2,1,[0,0,0]))

    set_exchange!(sys, [parameters["Jc"]*2.0   0.0   0.0;
                        0.0   parameters["Jc"]*2.0   0.0;
                        0.0   0.0   parameters["Jc"]*2.0], Bond(1,3,[0,0,0]))

    set_onsite_coupling!(sys, parameters["D"]*S[2]^2, 1)
end

my_parameters = Dict("Jab" => -0.2, "Jc" => 0.14, "D" => -0.04, "S" => 2.0)
my_parameters = Dict("Jab" => -0.6, "Jc" => 0.3, "D" => -0.5, "S" => 2.0)
param_to_optimize = Dict("Jab" => true, "Jc" => true, "D" => true, "S" => false)
param_bounds = Dict("Jab" => (-1.0, 0.0), "Jc" => (0.0, 1.0), "D" => (-1.0, 0.0))

swt = SpinWaveTheory(sys)

sunny_model = SunnyModel("Sunny", my_parameters, param_to_optimize, param_bounds, sunny_dispersion, sys, swt)

set_system_params!(sys, sunny_model.parameters.values)

refine_model!(sunny_model, LMO_data, true)

update_dispersion_and_residuals!(LMO_data, sunny_model, true)

plot_group!(LMO_data, sunny_model)