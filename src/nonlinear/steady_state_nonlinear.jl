"""
	steady_state(problem::NonLinearOptics{MeanField}; tmax = 250.0, reltol = 1e-11, abstol = 1e-10, m = 90, ode_solver = false)
"""
function steady_state(problem::NonLinearOptics{T}; tmax = 250.0, reltol = 1e-11, abstol = 1e-10, m = 90, ode_solver = false)  where {T<:NonLinear}

	if ode_solver
		tspan = (0.0, tmax)
		u₀ = default_initial_condition(problem)
		return time_evolution(problem, u₀, tspan; reltol = reltol, abstol = abstol, save_on = false).u[end]
    end

    G = copy(interaction_matrix(problem))
    Ωₙ = laser_field(problem.laser, problem.atoms)
    parameters = get_evolution_params(problem, G, Ωₙ)

    ## `nlsolve` convergence is senstive to initial conditions
    ## therefore, i decided to make a small time evoltuion, and use the result as initial condition
    ## NOTE: this trick is usefull only for small N, for now, is N < 1200
    ## but this number was NOT obtained by systematic tests, and could be improved
    u₀ = default_initial_condition(problem)
    problemFunction = get_evolution_function(problem)

    if problem.atoms.N < 1200
        tspan = (0.0, tmax)
        u₀_guess::Vector{ComplexF64} = time_evolution(problem, u₀, tspan; reltol = reltol, abstol = abstol, save_on = false).u[end] # evolve a little bit
    end

    try
        solution = nlsolve((du, u) -> problemFunction(du, u, parameters, 0.0), u₀_guess, method = :anderson, m = m, autodiff = :forward)
        return solution.zero
    catch
        ## For lower N (N < 100 ?), nlsolve does not converge (i don't know why)
        ## Instead of returning an error, I return the result from time evolution.
        @warn "Steady State may not be accurate. Consider increasing number of particles."
        return u₀_guess
    end
end