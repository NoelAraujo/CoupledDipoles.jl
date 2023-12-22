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
        my_prob = SteadyStateProblem(problemFunction, u₀, parameters)
        return solve(my_prob,  DynamicSS(VCABM(); reltol = reltol, abstol = abstol, tspan=tmax)).u
		# return time_evolution(problem, _u₀, (0, tmax); reltol = reltol, abstol = abstol, save_on = false).u[end]
    end

    try
        my_prob = SteadyStateProblem(problemFunction, u₀, parameters)
        solution = solve(my_prob, NewtonRaphson())
        return solution.u
    catch
        ## For lower N (N < 100 ?), nlsolve does not converge (i don't know why)
        ## Instead of returning an error, I return the result from time evolution.
        @warn "Steady State may not be accurate. Consider increasing number of particles."
        return u₀
    end
end