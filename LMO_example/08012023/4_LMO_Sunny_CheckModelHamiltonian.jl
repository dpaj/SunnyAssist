using Sunny, Optim, StaticArrays, DataFrames, CSV, LinearAlgebra, Plots, GLM, Statistics, StatsBase, Printf, GLMakie


function Q_maker(N::Int64, start::Vector{<:Real}, stop::Vector{<:Real})
    # Make sure that the start and stop vectors have the same length
    if length(start) != length(stop)
        error("The start and stop vectors must have the same length")
    end

    # Make sure that the start and stop vectors are three-dimensional
    if length(start) != 3
        error("The start and stop vectors must be three-dimensional")
    end

    # Convert the vectors to Vector{Float64}
    start = convert(Vector{Float64}, start)
    stop = convert(Vector{Float64}, stop)

    # Create an empty array of SVectors
    points = Vector{SVector{3, Float64}}(undef, N)

    # Calculate the step size
    step = (stop .- start) / (N - 1)

    # Generate the points
    for i in 1:N
        points[i] = SVector{3, Float64}(start .+ step .* (i - 1))
    end

    return points
end

mutable struct ExpData
    name::String
    Q::Vector{SVector{3, Float64}}
    E::Vector{SVector{1, Float64}}
    direction::Union{Vector{Float64}, String}
end

mutable struct ModelData
    name::String
    model_Q::Vector{SVector{3, Float64}}
    model_E::Vector{SVector{1, Float64}}
    residuals::Vector{SVector{1, Float64}}
    direction::Union{Vector{Float64}, String}
    model_Q_dense::Vector{SVector{3, Float64}}
    model_E_dense::Vector{SVector{1, Float64}}
end

mutable struct ExpDataGroup
    files::Vector{ExpData}
    model::Vector{ModelData}
    all_Q::Vector{SVector{3, Float64}}
    all_E::Vector{SVector{1, Float64}}
    model_all_E::Vector{SVector{1, Float64}}
end

abstract type AbstractCalculatorWrapper end

mutable struct SunnyCalculatorWrapper <: AbstractCalculatorWrapper
    dispersion::Function
    sys::System
    swt::SpinWaveTheory
end

mutable struct AnalyticalCalculatorWrapper <: AbstractCalculatorWrapper
    dispersion::Function
end

mutable struct ParameterWrapper
    values::Dict{String, Float64}
    param_to_optimize::Dict{String, Bool}
    param_bounds::Dict{String, Tuple{Float64, Float64}}
end

abstract type AbstractModel end

struct AnalyticalModel <: AbstractModel
    name::String
    parameters::ParameterWrapper
    calculator::AnalyticalCalculatorWrapper
end

struct SunnyModel <: AbstractModel
    name::String
    parameters::ParameterWrapper
    calculator::SunnyCalculatorWrapper
end

function AnalyticalModel(name::String, values::Dict{String, Float64}, param_to_optimize::Dict{String, Bool}, param_bounds::Dict{String, Tuple{Float64, Float64}}, dispersion_function::Function)
    parameters = ParameterWrapper(values, param_to_optimize, param_bounds)
    calculator = AnalyticalCalculatorWrapper(dispersion_function)
    AnalyticalModel(name, parameters, calculator)
end

function AnalyticalModel(name::String, values::Dict{String, Float64}, dispersion_function::Function)
    param_to_optimize = Dict("nothing" => false)
    param_bounds = Dict("nothing" => (0.0, 0.0))
    parameters = ParameterWrapper(values, param_to_optimize, param_bounds)
    calculator = AnalyticalCalculatorWrapper(dispersion_function)
    AnalyticalModel(name, parameters, calculator)
end

function SunnyModel(name::String, values::Dict{String, Float64}, param_to_optimize::Dict{String, Bool}, param_bounds::Dict{String, Tuple{Float64, Float64}}, dispersion_function::Function, sys::System, swt::SpinWaveTheory)
    parameters = ParameterWrapper(values, param_to_optimize, param_bounds)
    calculator = SunnyCalculatorWrapper(dispersion_function, sys, swt)
    SunnyModel(name, parameters, calculator)
end

function SunnyModel(name::String, values::Dict{String, Float64}, dispersion_function::Function, sys::System, swt::SpinWaveTheory)
    param_to_optimize = Dict("nothing" => false)
    param_bounds = Dict("nothing" => (0.0, 0.0))
    parameters = ParameterWrapper(values, param_to_optimize, param_bounds)
    calculator = SunnyCalculatorWrapper(dispersion_function, sys, swt)
    SunnyModel(name, parameters, calculator)
end

function dispersion(model::AbstractModel, q::SVector{3, Float64})
    if model.calculator isa SunnyCalculatorWrapper
        return model.calculator.dispersion(model.calculator.sys, model.calculator.swt, model.parameters.values, q)
    elseif model.calculator isa AnalyticalCalculatorWrapper
        return model.calculator.dispersion(model.parameters.values, q)
    end
end

function load_datafile(path::String, n_dense::Int64=100)
    df = CSV.read(path, DataFrame)
    Q = [SVector{3, Float64}(df[i, :H], df[i, :K], df[i, :L]) for i in 1:nrow(df)]
    E = [SVector{1, Float64}(df[i, :E]) for i in 1:nrow(df)]
    direction = find_direction(df)
    ExpData(path, Q, E, direction)
end

function load_modeldata(path::String, n_dense::Int64=100)
    df = CSV.read(path, DataFrame)
    model_Q = [SVector{3, Float64}(df[i, :H], df[i, :K], df[i, :L]) for i in 1:nrow(df)]
    model_E = [SVector{1, Float64}(0) for _ in 1:nrow(df)]
    residuals = [SVector{1, Float64}(0) for _ in 1:nrow(df)]
    direction = find_direction(df)
    model_Q_dense = [model_Q[1] + t*(model_Q[end] - model_Q[1]) for t in range(0, 1, length=n_dense)]
    model_E_dense = [SVector{1, Float64}(0) for _ in 1:n_dense]
    ModelData(path, model_Q, model_E, residuals, direction, model_Q_dense, model_E_dense)
end

function find_direction(df::DataFrame)
    data = hcat(df.H, df.K, df.L)
    centered_data = data .- mean(data, dims=1)
    covariance_matrix = cov(centered_data)
    eigenvalues, eigenvectors = eigen(covariance_matrix)
    _, max_index = findmax(eigenvalues)
    direction = eigenvectors[:, max_index]
    if any(isnan, direction)
        return "no well-defined direction vector"
    else
        direction = direction / abs(direction[findfirst(x -> x ≠ 0, direction)])
        return round.(direction, digits=2)
    end
end

function find_direction(q_points::Vector{SVector{3, Float64}})
    data = hcat([q[1] for q in q_points], [q[2] for q in q_points], [q[3] for q in q_points])
    centered_data = data .- mean(data, dims=1)
    covariance_matrix = cov(centered_data)
    eigenvalues, eigenvectors = eigen(covariance_matrix)
    _, max_index = findmax(eigenvalues)
    direction = eigenvectors[:, max_index]
    if any(isnan, direction)
        return "no well-defined direction vector"
    else
        direction = direction / abs(direction[findfirst(x -> x ≠ 0, direction)])
        return round.(direction, digits=2)
    end
end


function load_group(paths::Vector{String}, n_dense::Int64=100)
    files = [load_datafile(path, n_dense) for path in paths]
    model = [load_modeldata(path, n_dense) for path in paths]
    all_Q = reduce(vcat, [file.Q for file in files])
    all_E = reduce(vcat, [file.E for file in files])
    model_all_E = [SVector{1, Float64}(0) for _ in 1:length(all_E)]
    ExpDataGroup(files, model, all_Q, all_E, model_all_E)
end

function directional_parameter(file::ExpData, HKL::SVector{3, Float64})
    H, K, L = HKL
    dot(HKL - file.Q[1], file.direction)
end

function format_label(d::Vector{Float64}, hkl::SVector{3, Float64})
    hkl = Vector{Float64}(hkl) # Convert SVector to Vector
    direction = [x == 1.0 ? "ξ" : x == 0.0 ? "" : "$(x)ξ" for x in d]
    shift = [x == 0.0 ? "" : x == 1.0 ? "1" : @sprintf("%.2f", x) for x in hkl]
    labels = ["$(shift[i])" * (isempty(shift[i]) ? "" : " + ") * "$(direction[i])" for i in 1:3]
    labels = [isempty(l) ? "0" : l for l in labels]
    return labels
end

function plot_group!(group::ExpDataGroup, model::AbstractModel)
    plt = Plots.plot(layout = (1, length(group.files)), legend = false)
    for (i, file) in enumerate(group.files)
        x_values = [directional_parameter(file, file.Q[j]) for j in 1:length(file.Q)]
        y_values = [file.E[j][1] for j in 1:length(file.E)]
        
        # Adjusted to refer to group.model[i] instead of file
        model_y_values = [group.model[i].model_E[j][1] for j in 1:length(group.model[i].model_E)]
        residuals_values = [group.model[i].residuals[j][1] for j in 1:length(group.model[i].residuals)] 
        
        Plots.scatter!(plt[i], x_values, y_values, markershape = :circle, label="E")
        Plots.scatter!(plt[i], x_values, model_y_values, markershape = :square, color = :red, label="model_E")
        Plots.plot!(plt[i], x_values, residuals_values, color = :blue, linewidth=2, label="residuals") 
        
        # Overplot model_E_dense as a function of model_Q_dense
        x_values_dense = [directional_parameter(file, group.model[i].model_Q_dense[j]) for j in 1:length(group.model[i].model_Q_dense)]
        model_y_values_dense = [group.model[i].model_E_dense[j][1] for j in 1:length(group.model[i].model_E_dense)]
        Plots.plot!(plt[i], x_values_dense, model_y_values_dense, color = :black, linewidth=2)
    end

    for (i, file) in enumerate(group.files)
        xlabel = format_label(file.direction, file.Q[1])
        Plots.plot!(plt[i], xlabel=join(xlabel, ", "), ylabel="ℏω")
    end
    
    # Create a title with the model's optimized parameters
    title_str = "Optimized Parameters: " * join(["$(key) = $(round(value, digits=3))" for (key, value) in model.parameters.values], ", ")
    Plots.plot!(plt, title=title_str, subplot = 2, titlefont=font(12))
    
    display(plt)
end

function plot_model!(group::ExpDataGroup, model::AbstractModel)
    plt = Plots.plot(layout = (1, length(group.files)), legend = false)
    for (i, file) in enumerate(group.files)
        # X values are the directional parameters of each Q point
        x_values = [directional_parameter(file, file.Q[j]) for j in 1:length(file.Q)]

        # Y values are the model_E for each Q point
        model_y_values = [group.model[i].model_E[j][1] for j in 1:length(group.model[i].model_E)]
        
        # Scatter plot of model's E values
        Plots.scatter!(plt[i], x_values, model_y_values, markershape = :square, color = :red, label="model_E")
        
        # Overplot model_E_dense as a function of model_Q_dense
        x_values_dense = [directional_parameter(file, group.model[i].model_Q_dense[j]) for j in 1:length(group.model[i].model_Q_dense)]
        model_y_values_dense = [group.model[i].model_E_dense[j][1] for j in 1:length(group.model[i].model_E_dense)]
        Plots.plot!(plt[i], x_values_dense, model_y_values_dense, color = :black, linewidth=2)
        
        # Add labels to the plot
        xlabel = format_label(file.direction, file.Q[1])
        Plots.plot!(plt[i], xlabel=join(xlabel, ", "), ylabel="ℏω")
    end
    
    # Create a title with the model's parameters
    title_str = "Parameters: " * join(["$(key) = $(round(value, digits=3))" for (key, value) in model.parameters.values], ", ")
    Plots.plot!(plt, title=title_str, subplot = 2, titlefont=font(12))
    
    display(plt)
end


function plot_experiment!(group::ExpDataGroup)
    plt = Plots.plot(layout = (1, length(group.files)), legend = false)
    
    for (i, file) in enumerate(group.files)
        # X values are the directional parameters of each Q point
        x_values = [directional_parameter(file, file.Q[j]) for j in 1:length(file.Q)]
        
        # Y values are the E values for each Q point
        y_values = [file.E[j][1] for j in 1:length(file.E)]
        
        # Scatter plot of experimental E values
        Plots.scatter!(plt[i], x_values, y_values, markershape = :circle, label="E")
        
        # Add labels to the plot
        xlabel = format_label(file.direction, file.Q[1])
        Plots.plot!(plt[i], xlabel=join(xlabel, ", "), ylabel="ℏω")
    end
    
    display(plt)
end

function update_dispersion_simulation!(q::Vector{SVector{3, Float64}}, model::AnalyticalModel)
    for (i, q_val) in enumerate(q)
        model.model_E[i] = dispersion(model, q_val)
    end

    for (i, q_val) in enumerate(model.model_Q_dense)
        model.model_E_dense[i] = dispersion(model, q_val)
    end

end


function update_dispersion_and_residuals!(group::ExpDataGroup, model::AnalyticalModel, update_dense_model::Bool=false)
    for (i, q) in enumerate(group.all_Q)
        group.model_all_E[i] = dispersion(model, q)
    end

    for (i, m_file) in enumerate(group.model)
        for (j, q) in enumerate(m_file.model_Q)
            m_file.model_E[j] = dispersion(model, q)
            m_file.residuals[j] = group.files[i].E[j] - m_file.model_E[j] # Update residuals
        end
        if update_dense_model
            for (j, q) in enumerate(m_file.model_Q_dense)
                m_file.model_E_dense[j] = dispersion(model, q)
            end
        end
    end
end

function update_dispersion_and_residuals!(group::ExpDataGroup, model::SunnyModel, update_dense_model::Bool=false)
    set_system_params!(model.calculator.sys, model.parameters.values)
    model.calculator.swt = SpinWaveTheory(model.calculator.sys, energy_ϵ = 1e-5, energy_tol = 1e-6)
    for (i, q) in enumerate(group.all_Q)
        group.model_all_E[i] = dispersion(model, q)
    end

    for (i, m_file) in enumerate(group.model)
        for (j, q) in enumerate(m_file.model_Q)
            m_file.model_E[j] = dispersion(model, q)
            m_file.residuals[j] = group.files[i].E[j] - m_file.model_E[j] # Update residuals
        end
        if update_dense_model
            for (j, q) in enumerate(m_file.model_Q_dense)
                m_file.model_E_dense[j] = dispersion(model, q)
            end
        end
    end
end


function objective_function(params::Vector{Float64}, param_names::Vector{String}, model::AbstractModel, group::ExpDataGroup)
    for (i, name) in enumerate(param_names)
        model.parameters.values[name] = params[i]
    end
    
    update_dispersion_and_residuals!(group, model, false)  # This will dispatch to the correct function based on the model type
    
    return sum(sum(m_file.residuals[i][1]^2 for i in 1:length(m_file.residuals)) for m_file in group.model)
end

function refine_model!(model::AbstractModel, group::ExpDataGroup, verbose::Bool = false)
    param_names = String[]
    initial_params = Float64[]
    lower_bounds = Float64[]  # Vector to hold lower bounds
    upper_bounds = Float64[]  # Vector to hold upper bounds

    for (name, should_optimize) in model.parameters.param_to_optimize
        if should_optimize
            push!(param_names, name)
            push!(initial_params, model.parameters.values[name])
            bounds = get(model.parameters.param_bounds, name, (-Inf, Inf))  # Use (-Inf, Inf) if no bounds specified
            push!(lower_bounds, bounds[1])
            push!(upper_bounds, bounds[2])
        end
    end
    
    function_counter = Ref(0)
    obj_fn(params) = (function_counter[] += 1; objective_function(params, param_names, model, group))


    callback_fn(optimizer_state) = (verbose && println(optimizer_state); return false)

    #callback_fn(optimizer_state) = (verbose && (println("Iteration: ", optimizer_state.iteration); println("Current parameter values: ", optimizer_state.metadata["x"]); println("Current objective function value: ", optimizer_state.value)); return false)

    inner_optimizer = LBFGS()
    
    optim_options = Optim.Options(
        g_tol = 1e-12,
        x_tol = 1e-12,
        f_tol = 1e-12,
        show_trace = verbose,
        show_every = verbose ? 1 : 0,  # 0 means never show
        iterations = 1000,
        callback = callback_fn,
    )

    optimizer = Fminbox(inner_optimizer)
    
    start_time = time()
    elapsed_time = @elapsed optimal_solution = optimize(obj_fn, lower_bounds, upper_bounds, initial_params, optimizer, optim_options)
    end_time = time()

    optimized_parameters = Optim.minimizer(optimal_solution)
    for (i, name) in enumerate(param_names)
        model.parameters.values[name] = optimized_parameters[i]
    end

    println("Total number of function calls: ", function_counter[])
    println("Total time taken: ", elapsed_time, " seconds.")
end

function empty_experiment_group(q_points_list::Vector{Vector{SVector{3, Float64}}}, n_dense::Int64=100)
    # Create one ExpData and ModelData object for each vector of q-points
    files = []
    model = []
    for q_points in q_points_list
        direction = find_direction(q_points)
        push!(files, ExpData("Empty Experiment", q_points, [SVector{1, Float64}(0) for _ in q_points], direction))
        push!(model, ModelData("Empty Model", q_points, [SVector{1, Float64}(0) for _ in q_points], [SVector{1, Float64}(0) for _ in q_points], direction, q_points, [SVector{1, Float64}(0) for _ in q_points]))
    end

    # The combined all_Q, all_E, and model_all_E
    all_Q = reduce(vcat, [file.Q for file in files])
    all_E = reduce(vcat, [file.E for file in files])
    model_all_E = [SVector{1, Float64}(0) for _ in 1:length(all_E)]

    # Return the "empty" ExpDataGroup
    return ExpDataGroup(files, model, all_Q, all_E, model_all_E)
end


function objective_function_mapper(model::AbstractModel, group::ExpDataGroup, scalar::Float64, n_points::Int64 = 100)
    # Make a deep copy of the model to prevent modifying the original
    model_copy = deepcopy(model)

    param_names = [name for (name, should_optimize) in model_copy.parameters.param_to_optimize if should_optimize]
    n_params = length(param_names)
    
    # Initial parameters
    initial_params = [model_copy.parameters.values[name] for name in param_names]

    # Grid for each parameter
    param_grids = [range(param - scalar*abs(param), param + scalar*abs(param), length=n_points) for param in initial_params]
    
    # Initialize heatmap data storage
    heatmaps = Dict{Tuple{String, String}, Matrix{Float64}}()
    param_ranges = Dict{String, AbstractRange}()

    # Iterate over each pair of parameters
    for i in 1:n_params
        for j in (i+1):n_params

            model_copy = deepcopy(model)

            name1 = param_names[i]
            name2 = param_names[j]
            heatmap = zeros(n_points, n_points)

            # Evaluate the objective function at each point in the grid
            for (i1, param1) in enumerate(param_grids[i])
                for (i2, param2) in enumerate(param_grids[j])
                    # Construct the parameter vector with the new parameter values
                    param_values = copy(initial_params)
                    param_values[i] = param1
                    param_values[j] = param2

                    # Compute the objective function value
                    heatmap[i1, i2] = objective_function(param_values, param_names, model_copy, group)
                end
            end
            
            # Store the heatmap
            heatmaps[(name1, name2)] = heatmap
        end

        # Store the parameter ranges
        param_ranges[param_names[i]] = param_grids[i]
    end

    return heatmaps, param_ranges
end


function objective_function_mapper(model::AbstractModel, group::ExpDataGroup, n_points::Int64 = 100)
    param_names = [name for (name, should_optimize) in model.parameters.param_to_optimize if should_optimize]
    n_params = length(param_names)

    # Grid for each parameter
    param_grids = [range(model.parameters.param_bounds[name][1], model.parameters.param_bounds[name][2], length=n_points) for name in param_names]

    # Initialize heatmap data storage
    heatmaps = Dict{Tuple{String, String}, Matrix{Float64}}()
    param_ranges = Dict{String, AbstractRange}()

    # Iterate over each pair of parameters
    for i in 1:n_params
        for j in (i+1):n_params
            # Make a deep copy of the model to prevent modifying the original
            model_copy = deepcopy(model)
            
            name1 = param_names[i]
            name2 = param_names[j]
            heatmap = zeros(n_points, n_points)

            # Evaluate the objective function at each point in the grid
            for (i1, param1) in enumerate(param_grids[i])
                for (i2, param2) in enumerate(param_grids[j])
                    # Construct the parameter vector with the new parameter values
                    param_values = [model_copy.parameters.values[name] for name in param_names]
                    param_values[i] = param1
                    param_values[j] = param2

                    # Compute the objective function value
                    heatmap[i1, i2] = objective_function(param_values, param_names, model_copy, group)
                end
            end

            # Store the heatmap
            heatmaps[(name1, name2)] = heatmap
        end

        # Store the parameter ranges
        param_ranges[param_names[i]] = param_grids[i]
    end

    return heatmaps, param_ranges
end



function plot_all_heatmaps(heatmaps::Dict{Tuple{String, String}, Array{Float64, 2}}, param_ranges::Dict{String, AbstractRange}, model::AbstractModel)
    for ((param1, param2), heatmap_data) in heatmaps
        p = heatmap(param_ranges[param1], param_ranges[param2], heatmap_data,
                    title="Objective Function Heatmap for $param1 and $param2",
                    xlabel=param1, ylabel=param2, color=:viridis)
                    
        # Add a point for the optimal solution
        optimal_params = model.parameters.values
        scatter!([optimal_params[param1]], [optimal_params[param2]], color=:red, markersize=10, label="Optimal solution")
                    
        # Save the plot to a file. The filename is generated from the parameter names.
        #savefig(p, "$(param1)_$(param2)_heatmap.png")
        display(p)
    end
end


################################################################

function set_initial_spin_configuration!(sys::System, original_spins::Vector{Vector{Int}})
    # Set the spins
    for i in 1:length(sys.dipoles)
        polarize_spin!(sys, original_spins[i], CartesianIndex(1, 1, 1, i))# Cartesian indices for the spins
    end
end

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

sys = System(cryst, (1,1,1), [SpinInfo(1,S=2)], :dipole)

original_spins = [[0, -1, 0], [0, -1, 0], [0, 1, 0], [0, 1, 0]]
set_initial_spin_configuration!(sys, original_spins)

print_symmetry_table(cryst, 5.0)

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

function sunny_dispersion(sys::System, swt::SpinWaveTheory, parameters::Dict{String, Float64}, q::SVector{3, Float64})

    function get_index_of_max_norm(Sαβ)
        norms = [norm(Sαβ[i]) for i in 1:length(Sαβ)]
        index_of_max_norm = argmax(norms)
        return index_of_max_norm
    end

    disp, Sαβ = dssf(swt, [q])
    index = get_index_of_max_norm(Sαβ)

    return SVector{1, Float64}(disp[index])
end

my_parameters = Dict("Jab" => -0.6, "Jc" => 0.3, "D" => -0.5, "S" => 2.0)

swt = SpinWaveTheory(sys)

sunny_model = SunnyModel("Sunny", my_parameters, sunny_dispersion, sys, swt)

set_system_params!(sys, sunny_model.parameters.values)

q_L = Q_maker(100, [0, 0, 0], [0, 0, 1])
q_H = Q_maker(100, [0, 0, 0], [1, 0, 0])
empty_experiment = empty_experiment_group([q_L, q_H])

update_dispersion_and_residuals!(empty_experiment, sunny_model, true)

plot_model!(empty_experiment, sunny_model)