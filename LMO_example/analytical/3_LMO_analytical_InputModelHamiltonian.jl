using Sunny, Optim, StaticArrays, DataFrames, CSV, LinearAlgebra, Plots, GLM, Statistics, StatsBase, Printf

if Sys.iswindows()
    SunnyAssist_path = "C:\\Users\\vdp\\Dropbox (ORNL)\\Sunny\\SunnyAssist\\"
else
    SunnyAssist_path = "/home/vdp/Dropbox (ORNL)/Sunny/SunnyAssist/"
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