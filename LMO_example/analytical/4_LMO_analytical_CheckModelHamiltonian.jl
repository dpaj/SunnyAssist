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

my_parameters = Dict("Jab" => -0.2, "Jc" => 0.14, "D" => -0.04, "S" => 2.0)

analytical_model = AnalyticalModel("analytical", my_parameters, analytical_dispersion)

q_L = Q_maker(100, [0, 0, 0], [0, 0, 1])
q_H = Q_maker(100, [0, 0, 0], [1, 0, 0])
empty_experiment = empty_experiment_group([q_L, q_H])


update_dispersion_and_residuals!(empty_experiment, analytical_model, true)

plot_model!(empty_experiment, analytical_model)
