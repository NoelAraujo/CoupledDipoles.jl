"""
    steady_state(problem::LinearOptics{Scalar})

Solve `x=G\\Ω`, with default `interaction_matrix` and `laser_field`.
"""
function steady_state(problem::LinearOptics{Scalar}; tmax=250.0, reltol=1e-9, abstol=1e-9, ode_solver=false)

    if problem.atoms.N > 1
        if  ode_solver
            tspan = (0.0, tmax)
            u₀ = default_initial_condition(problem)
            βₛ = time_evolution(problem, u₀, tspan; reltol=reltol, abstol=abstol, save_on=false).u[end]
            return βₛ
        else
            G = interaction_matrix(problem)
            Ωₙ = vec(laser_field(problem, problem.atoms.r))
            βₛ = -(G \ Ωₙ)
        end
    else ## single atom
        G = interaction_matrix(problem)

        # the conversion from matrix to vector, avoids an extra function to deal single atom
        Ωₙ = laser_field(problem, problem.atoms.r)
        βₛ =   [-(Ωₙ / G[1])]   # I need a vector for single element
    end
    return βₛ
end


"""
    steady_state(problem::LinearOptics{Vectorial})

Solve `x=-G\\Ω`, with default `interaction_matrix` and `laser_field`. The solution x is reshaped as a 3xN matrix.
"""
function steady_state(problem::LinearOptics{Vectorial}; tmax=250.0, reltol=1e-9, abstol=1e-9, ode_solver=false)

    if problem.atoms.N > 1
        if  ode_solver
            tspan = (0.0, tmax)
            u₀ = default_initial_condition(problem)
            βₛ = time_evolution(problem, u₀, tspan; reltol=reltol, abstol=abstol, save_on=false).u[end] # evolve a little bit
            return _vecAux_longArray_into_Matrix(problem.atoms.N, βₛ)
        else
            G = interaction_matrix(problem)
            Ωₙ = LASER_FACTOR.*laser_field(problem, problem.atoms)
            Ωₙ_eff = _vecAux_Matrix_into_longArray(Ωₙ)

            βₛ = G \ Ωₙ_eff
        end
    else ## single atom
        G = interaction_matrix(problem)
        Ωₙ = laser_field(problem, problem.atoms)
        βₛ = -(Ωₙ / G[1])
    end

    βₛ_eff = _vecAux_longArray_into_Matrix(problem.atoms.N, βₛ)
    return βₛ_eff
end