"""
	steady_state(problem::NonLinearOptics{MeanField}; tmax = 250.0, reltol = 1e-11, abstol = 1e-10, m = 90, ode_solver = false)
"""
function steady_state(problem::NonLinearOptics{T}; tmax = 250.0, reltol = 1e-11, abstol = 1e-11, m = 90, ode_solver = false)  where {T<:NonLinear}

    G = copy(interaction_matrix(problem))
    Ωₙ = laser_field(problem.laser, problem.atoms)
    parameters = get_evolution_params(problem, G, Ωₙ)
    _u₀ = default_initial_condition(problem)

    u₀::Vector{ComplexF64} = time_evolution(problem, _u₀, (0, 25.0); reltol = reltol, abstol = abstol, save_on = false).u[end]
    problemFunction = get_evolution_function(problem)

	if ode_solver        
		return time_evolution(problem, _u₀, (0, tmax); reltol = reltol, abstol = abstol, save_on = false).u[end]
    end

    try
        my_prob = SteadyStateProblem(problemFunction, u₀, parameters)
        solution = solve(my_prob, NewtonRaphson())
        return solution.u
    catch
        @warn "Steady State may not be accurate. Consider setting `ode_solver = true`."
        return u₀
    end
end