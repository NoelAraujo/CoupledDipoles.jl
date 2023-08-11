"""
    default_initial_condition(::Scalar) = zeros(ComplexF64, N)
"""
function default_initial_condition(problem::LinearOptics{Scalar})
    return zeros(ComplexF64, problem.atoms.N) # I must use "zeros" and NOT an undef array - with trash data inside
end


"""
    default_initial_condition(::Vectorial) = zeros(ComplexF64, 3N)
"""
function default_initial_condition(problem::LinearOptics{Vectorial})
    return zeros(ComplexF64, 3problem.atoms.N) # I must use "zeros" and NOT an undef array - with trash data inside
end


"""
    steady_state(problem::LinearOptics{Scalar})

Solve `x=G\\Ω`, with default `interaction_matrix` and `laser_field`.
"""
function steady_state(problem::LinearOptics{Scalar}; tmax=250.0, reltol=1e-7, abstol=1e-7, bruteforce=false)
    G = interaction_matrix(problem)
    if problem.atoms.N > 1
        if  bruteforce
            tspan = (0.0, tmax)
            u₀ = default_initial_condition(problem)
            βₛ = time_evolution(problem, u₀, tspan; reltol=reltol, abstol=abstol, save_on=false).u[end] # evolve a little bit
            return βₛ
        else
            Ωₙ = vec(laser_field(problem, problem.atoms.r))
            βₛ = -(G \ Ωₙ)
        end
    else
        # the conversion from matrix to vector, avoids an extra function to deal single atom
        Ωₙ = laser_field(problem, vec(problem.atoms.r))
        βₛ =   [-(Ωₙ / G[1])]   # I need a vector for single element
    end
    return βₛ
end


"""
    steady_state(problem::LinearOptics{Vectorial})

Solve `x=-G\\Ω`, with default `interaction_matrix` and `laser_field`. The solution x is reshaped as a 3xN matrix.
"""
function steady_state(problem::LinearOptics{Vectorial}; tmax=250.0, reltol=1e-7, abstol=1e-7, bruteforce=false)
    G = interaction_matrix(problem)
    Ωₙ = laser_field(problem, problem.atoms.r)
    Ωₙ_eff = _vecAux_Matrix_into_longArray(Ωₙ)

    if problem.atoms.N > 1
        if  bruteforce
            tspan = (0.0, tmax)
            u₀ = default_initial_condition(problem)
            βₛ = time_evolution(problem, u₀, tspan; reltol=reltol, abstol=abstol, save_on=false).u[end] # evolve a little bit
            return βₛ
        else
            βₛ = -(G \ Ωₙ_eff)
        end
    else
        βₛ = -(Ωₙ / G[1])
    end

    βₛ_eff = _vecAux_longArray_into_Matrix(problem.atoms.N, βₛ)
    return βₛ_eff
end

function time_evolution(
    problem::LinearOptics{T},
    u₀,
    tspan::Tuple;
    bruteForce = false,
    kargs...,
) where {T<:Linear}
    ### if time is big, and laser is swithc-off, the 'FORMAL SOLUTION' of the ODE problem is much faster
    if problem.laser.s ≈ 0 && bruteForce==true
        time_interval = get(kargs, :saveat, range(tspan[1], tspan[2], length = 20))
        solution = time_evolution_laser_off(problem, u₀, time_interval)
        return solution
    elseif bruteForce==true
        time_interval = get(kargs, :saveat, range(tspan[1], tspan[2], length = 20))
        solution = time_evolution_laser_on(problem, u₀, time_interval)
        return solution
    end

    ### use default G and Ωₙ
    G = copy(interaction_matrix(problem))
    Ωₙ = laser_field(problem, problem.atoms.r)
    solution = time_evolution(problem, u₀, tspan, Ωₙ, G; kargs...)

    return solution
end
function time_evolution(
    problem::LinearOptics{T},
    u₀,
    tspan::Tuple,
    Ωₙ::VecOrMat;
    kargs...,
) where {T<:Linear}
    ### use default G
    G = copy(interaction_matrix(problem))
    solution = time_evolution(problem, u₀, tspan, Ωₙ, G; kargs...)

    return solution
end
function time_evolution(
    problem::LinearOptics{T},
    u₀,
    tspan::Tuple,
    Ωₙ::VecOrMat,
    G::Matrix;
    kargs...,
) where {T<:Linear}
    @debug "start: time evolution - LinearOptics"

    ### parameters == constant vector and matrices
    parameters = get_evolution_params(problem, G, Ωₙ)

    ### calls for solver
    problemFunction = get_evolution_function(problem)
    prob = ODEProblem(problemFunction, u₀, tspan, parameters)
    solution = OrdinaryDiffEq.solve(
        prob,
        VCABM3();
        kargs...,
    )

    @debug "end  : time evolution - LinearOptics"
    return solution
end

get_evolution_function(problem::LinearOptics{Scalar}) = Scalar!
get_evolution_params(problem::LinearOptics{Scalar}, G, Ωₙ) = view(G, :, :), view(Ωₙ, :)

# MUDAR AQUI
get_evolution_function(problem::LinearOptics{Vectorial}) = Scalar!
function get_evolution_params(problem::LinearOptics{Vectorial}, G, Ωₙ)
    ## Ωₙ_eff  = [all X - all Y - all Z]
    Ωₙ_eff = vcat(view(Ωₙ, 1, :), view(Ωₙ, 2, :), view(Ωₙ, 3, :))
    return view(G, :, :), view(Ωₙ_eff, :)
end
function Scalar!(du, u, p, t)
    G, Ωₙ = p

    #=
    Equivalent to:
        "du[:] = G*u + Ωₙ"
    But using inplace operation
    =#
    du .= muladd(G, u, Ωₙ)
    return nothing
end


## right now, DifferentialEquaions is faster
### TO DO:
## needs optimization with 'laser_contribution' term
## needs to discover the regime where is good to use (best N? best time interval? best number time steps?)
function time_evolution_laser_on(
    problem,
    initial_state::AbstractVector,
    time_interval::AbstractVector;
    useBigFloat=false
)
    N = problem.atoms.N
    G = interaction_matrix(problem)
    Ωₙ = laser_field(problem, problem.atoms.r)

    Λ, P = eigen(G)
    Pinv = inv(P)
    Λinv = inv(Diagonal(Λ))
    PinvΩₙ = Pinv*vec(Ωₙ)

    if useBigFloat
        Λ = BigFloat.(real.(Λ)) + BigFloat.(imag.(Λ))*im
        Pinv = BigFloat.(real.(Pinv)) + BigFloat.(imag.(Pinv))*im
        Λinv = BigFloat.(real.(Λinv)) + BigFloat.(imag.(Λinv))*im
        PinvΩₙ = BigFloat.(real.(PinvΩₙ)) + BigFloat.(imag.(PinvΩₙ))*im
    end
    u = ThreadsX.map(time_interval) do t
        term1 = (P * Diagonal(exp.(+t * Λ)) * Pinv)*initial_state
        term2 = P * Diagonal(exp.(+t * Λ)) * (  Λinv*(Diagonal(exp.(-t * Λ)) - I )*PinvΩₙ    )
        term1 + term2
    end
    return (t=time_interval, u=u)
end

### NOT INTENDED TO BE EXPOSED
function time_evolution_laser_off(
    problem,
    initial_state::AbstractVector,
    time_interval::AbstractVector,
)
    G = interaction_matrix(problem)
    Λ, P = eigen(G)
    Pinv = inv(P)

    u = ThreadsX.map(time_interval) do t
        P * Diagonal(exp.(t * Λ)) * Pinv * initial_state
    end
    return (t=time_interval, u=u)
end

# this function is executed when i compute the vectorial model
function time_evolution_laser_off(
    problem,
    initial_state::AbstractMatrix,
    time_interval::AbstractVector,
)
    initial_state_array =_vecAux_Matrix_into_longArray(initial_state)
    time_evolution_laser_off(problem, initial_state_array, time_interval)
end
