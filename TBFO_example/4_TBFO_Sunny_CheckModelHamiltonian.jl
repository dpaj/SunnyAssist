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

function set_system_params!(sys::System, parameters::Dict{String, Float64})
    S_Fe1 = spin_operators(sys, 1)
    S_Fe2 = spin_operators(sys, 3)

    set_exchange!(sys, [parameters["JcS"]   0.0   0.0;
                        0.0   parameters["JcS"]   0.0;
                        0.0   0.0   parameters["JcS"]], Bond(1, 4, [0, 0, 0]))

    set_exchange!(sys, [parameters["Jb"]   0.0   0.0;
                        0.0   parameters["Jb"]   0.0;
                        0.0   0.0   parameters["Jb"]], Bond(1, 1, [0, 1, 0]))

    set_exchange!(sys, [parameters["Jb"]   0.0   0.0;
                        0.0   parameters["Jb"]   0.0;
                        0.0   0.0   parameters["Jb"]], Bond(3, 3, [0, 1, 0]))

    set_exchange!(sys, [parameters["JcL"]   0.0   0.0;
                        0.0   parameters["JcL"]   0.0;
                        0.0   0.0   parameters["JcL"]], Bond(2, 3, [0, 0, 1]))

    set_exchange!(sys, [parameters["Ja"]   0.0   0.0;
                        0.0   parameters["Ja"]   0.0;
                        0.0   0.0   parameters["Ja"]], Bond(1, 3, [0, 0, 0]))

    set_onsite_coupling!(sys, parameters["Db"]*S_Fe1[2]^2, 1)
    set_onsite_coupling!(sys, parameters["Db"]*S_Fe2[2]^2, 3)
end

sunny_parameters = Dict("Ja" => 3.17, "Jb" => 7.31, "JcS" => 0.57, "JcL" => 5.11, "Db" => -0.47)

swt = SpinWaveTheory(sys)

sunny_model = SunnyModel("Sunny", sunny_parameters, sunny_dispersion, sys, swt)

set_system_params!(sys, sunny_model.parameters.values)

q_L = Q_maker(100, [0, 0, 0], [0, 0, 1])
q_H = Q_maker(100, [0, 0.5, 1], [1, 0.5, 1])
empty_experiment = empty_experiment_group([q_L, q_H])

update_dispersion_and_residuals!(empty_experiment, sunny_model, true)

plot_model!(empty_experiment, sunny_model)

#print_wrapped_intensities(sys)

# Select a sequence of wavevectors that will define a piecewise linear
# interpolation in reciprocal lattice units (RLU).

q_points = [[0,0,0], [1,0,0], [0,1,0], [1/2,0,0], [0,1,0], [0,0,0]];

# The function [`reciprocal_space_path`](@ref) will linearly sample a `path`
# between the provided $q$-points with a given `density`. The `xticks` return
# value provides labels for use in plotting.

density = 50
path, xticks = reciprocal_space_path(cryst, q_points, density);

# The [`dispersion`](@ref) function defines the quasiparticle excitation
# energies $ω_i(𝐪)$ for each point $𝐪$ along the reciprocal space path.

disp = dispersion(swt, path);


# In addition to the band energies $ω_i(𝐪)$, Sunny can calculate the inelastic
# neutron scattering intensity $I_i(𝐪)$ for each band $i$ according to an
# [`intensity_formula`](@ref). We choose to apply a polarization correction
# $(1 - 𝐪⊗𝐪)$ by setting the mode argument to `:perp`. Selecting
# `delta_function_kernel` specifies that we want the energy and intensity of
# each band individually.

formula = intensity_formula(swt, :perp; kernel=delta_function_kernel)

# The function [`intensities_bands`](@ref) uses linear spin wave theory to
# calculate both the dispersion and intensity data for the provided path.

disp, intensity = intensities_bands(swt, path, formula);

# These can be plotted in GLMakie.

fig = Figure()
ax = Axis(fig[1,1]; xlabel="𝐪", ylabel="Energy (meV)", xticks, xticklabelrotation=π/6)
GLMakie.ylims!(ax, 0.0, 50)
GLMakie.xlims!(ax, 1, size(disp, 1))
colorrange = extrema(intensity)
for i in axes(disp)[2]
    GLMakie.lines!(ax, 1:length(disp[:,i]), disp[:,i]; color=intensity[:,i], colorrange)
end
fig


# To make comparisons with inelastic neutron scattering (INS) data, it is
# helpful to employ an empirical broadening kernel, e.g., a
# [`lorentzian`](@ref).

γ = 0.15 # width in meV
broadened_formula = intensity_formula(swt, :perp; kernel=lorentzian(γ))

# The [`intensities_broadened`](@ref) function requires an energy range in
# addition to the $𝐪$-space path.

energies = collect(0:0.01:10)  # 0 < ω < 10 (meV).
is1 = intensities_broadened(swt, path, energies, broadened_formula);

# A real FeI$_2$ sample will exhibit competing magnetic domains associated with
# spontaneous symmetry breaking of the 6-fold rotational symmetry of the
# triangular lattice. Note that the wavevectors $𝐪$ and $-𝐪$ are equivalent in
# the structure factor, which leaves three distinct domain orientations, which
# are related by 120° rotations about the $ẑ$-axis. Rather than rotating the
# spin configuration directly, on can rotate the $𝐪$-space path. Below, we use
# [`rotation_in_rlu`](@ref) to average the intensities over all three possible
# orientations.

R = rotation_in_rlu(cryst, [0, 0, 1], 2π/3)
is2 = intensities_broadened(swt, [R*q for q in path], energies, broadened_formula)
is3 = intensities_broadened(swt, [R*R*q for q in path], energies, broadened_formula)
is_averaged = (is1 + is2 + is3) / 3

fig = Figure()
ax = Axis(fig[1,1]; xlabel="(H,0,0)", ylabel="Energy (meV)", xticks, xticklabelrotation=π/6)
GLMakie.heatmap!(ax, 1:size(is_averaged, 1), energies, is_averaged)
fig
