using Sunny, Optim, StaticArrays, DataFrames, CSV, LinearAlgebra, Plots, GLM, Statistics, StatsBase, Printf

include(joinpath(SunnyAssist_path, "SunnyAssist.jl"))


# Here is the data path
data_path = joinpath(SunnyAssist_path, "LMO_example")

# Here are your data files
data_files = ["LaMnO3_dispersion_a.csv", "LaMnO3_dispersion_b.csv", "LaMnO3_dispersion_c.csv"]

# Prepend the data path to each file
data_files = [joinpath(data_path, file) for file in data_files]

LMO_data = load_group(data_files)

plot_experiment!(LMO_data)
