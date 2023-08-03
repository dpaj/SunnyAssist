using Sunny, Optim, StaticArrays, DataFrames, CSV, LinearAlgebra, Plots, GLM, Statistics, StatsBase, Printf

if Sys.iswindows()
    SunnyAssist_path = "C:\\Users\\vdp\\Dropbox (ORNL)\\Sunny\\"
else
    SunnyAssist_path = "/home/vdp/Dropbox (ORNL)/Sunny/"
end

include(joinpath(SunnyAssist_path, "SunnyAssist.jl"))


function analytical_dispersion(parameters::Dict{String, Float64}, q::SVector{3, Float64})
    Jab = parameters["Jab"]
    Jc = parameters["Jc"]
    D = parameters["D"]
    S = parameters["S"]

    A_q = 2*Jab*(2-cos(pi*(q[1]+q[2]))-cos(pi*(q[1]-q[2])))-2*Jc+D
    B_q = -2*Jc*cos(pi*q[3])

    if A_q^2 - B_q^2 >= 0
        return SVector{1, Float64}(2*S*sqrt(A_q^2-B_q^2))
    else
        return SVector{1, Float64}(0)
    end
end


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

my_parameters = Dict("Jab" => -0.2, "Jc" => 0.14, "D" => -0.04, "S" => 2.0)
my_parameters = Dict("Jab" => -0.6, "Jc" => 0.3, "D" => -0.5, "S" => 2.0)
param_to_optimize = Dict("Jab" => true, "Jc" => true, "D" => true, "S" => false)
param_bounds = Dict("Jab" => (-1.0, 0.0), "Jc" => (0.0, 1.0), "D" => (-1.0, 0.0))

analytical_model = AnalyticalModel("analytical", my_parameters, param_to_optimize, param_bounds, analytical_dispersion)

refine_model!(analytical_model, LMO_data, true)

update_dispersion_and_residuals!(LMO_data, analytical_model, true)

plot_group!(LMO_data, analytical_model)

@time begin
    objective_function_maps, param_ranges = objective_function_mapper(analytical_model, LMO_data, 20)
    plot_all_heatmaps(objective_function_maps, param_ranges, analytical_model)
end

@time begin
    objective_function_maps, param_ranges = objective_function_mapper(analytical_model, LMO_data, 2.0, 20)
    plot_all_heatmaps(objective_function_maps, param_ranges, analytical_model)
end
