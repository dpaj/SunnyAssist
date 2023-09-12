using Sunny, Optim, StaticArrays, DataFrames, CSV, LinearAlgebra, Plots, GLM, Statistics, StatsBase, Printf

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
data_path = joinpath(SunnyAssist_path, "LMO_example")

# Here are your data files
data_files = ["LaMnO3_dispersion_a.csv", "LaMnO3_dispersion_b.csv", "LaMnO3_dispersion_c.csv"]

# Prepend the data path to each file
data_files = [joinpath(data_path, file) for file in data_files]

LMO_data = load_group(data_files)
