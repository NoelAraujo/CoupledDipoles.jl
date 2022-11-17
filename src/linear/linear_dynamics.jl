"""
    default_initial_condition(::Scalar) = zeros(ComplexF64, N)
"""
function default_initial_condition(problem::LinearOptics{Scalar})
    return zeros(ComplexF64, problem.atoms.N) # I must use "zeros" and NOT an undef array - with trash data inside
end

"""
    steady_state(problem::LinearOptics{Scalar})

Solve `x=G\\Ω`, with default `interaction_matrix` and `laser_field`.
"""
function steady_state(problem::LinearOptics{Scalar})
    G = interaction_matrix(problem)
    Ωₙ = laser_field(problem, problem.atoms.r)
    if problem.atoms.N > 1
        βₛ = -(G \ Ωₙ)
    else
        # The negative sign was NOT forgotten
        # after some math, I verified that it does not exist
        βₛ = Ωₙ / G[1]
    end
end


"""
    steady_state(problem::LinearOptics{Vectorial})

Solve `x=G\\Ω`, with default `interaction_matrix` and `laser_field`. The solution x is reshaped as a 3xN matrix.
"""
function steady_state(problem::LinearOptics{Vectorial})
    G = interaction_matrix(problem)
    Ωₙ = laser_field(problem, problem.atoms.r)
    ## Ωₙ_eff  = [all X - all Y - all Z]
    Ωₙ_eff = vcat(view(Ωₙ, 1, :), view(Ωₙ, 3, :), view(Ωₙ, 3, :))
    if problem.atoms.N > 1
        βₛ = -(G \ Ωₙ_eff)
    else
        # The negative sign was NOT forgotten
        # after some math, I verified that it does not exist
        βₛ = Ωₙ / G[1]
    end
    return reshape(βₛ, 3, problem.atoms.N)
end

function time_evolution(
    problem::LinearOptics{T},
    u₀,
    tspan::Tuple;
    kargs...,
) where {T<:Linear}
    ### if time is big, and laser is swithc-off, the 'FORMAL SOLUTION' of the ODE problem is much faster
    if problem.laser.s ≈ 0
        time_interval = get(kargs, :saveat, range(tspan[1], tspan[2], length = 20))
        states = time_evolution_laser_off(problem, u₀, time_interval)
        solution = (t = time_interval, u = states)
        return solution
    end

    ### use default G and Ωₙ
    G = copy(interaction_matrix(problem))
    Ωₙ = laser_field(problem, problem.atoms)
    solution = time_evolution(problem, u₀, tspan, Ωₙ, G; kargs...)

    return solution
end
function time_evolution(
    problem::LinearOptics{T},
    u₀,
    tspan::Tuple,
    Ωₙ::Vector;
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
    Ωₙ::Vector,
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
        VCABM();
        dt = 1e-10,
        abstol = 1e-10,
        reltol = 1e-10,
        kargs...,
    )

    @debug "end  : time evolution - LinearOptics"
    return solution
end

get_evolution_function(problem::LinearOptics{Scalar}) = Scalar!
get_evolution_params(problem::LinearOptics{Scalar}, G, Ωₙ) = view(G, :, :), view(Ωₙ, :)

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
    Ωₙ = laser_field(problem, problem.atoms)

    Λ, P = eigen(G)
    Pinv = inv(P)
    Λinv = inv(Diagonal(Λ))
    PinvΩₙ = Pinv*Ωₙ

    if useBigFloat
        Λ = BigFloat.(real.(Λ)) + BigFloat.(imag.(Λ))*im
        Pinv = BigFloat.(real.(Pinv)) + BigFloat.(imag.(Pinv))*im
        Λinv = BigFloat.(real.(Λinv)) + BigFloat.(imag.(Λinv))*im
        PinvΩₙ = BigFloat.(real.(PinvΩₙ)) + BigFloat.(imag.(PinvΩₙ))*im
    end
    z = ThreadsX.map(time_interval) do t
        term1 = (P * Diagonal(exp.(+t * Λ)) * Pinv)*initial_state
        term2 = P * Diagonal(exp.(+t * Λ)) * (  Λinv*(Diagonal(exp.(-t * Λ)) - I )*PinvΩₙ    )
        term1 + term2
    end
    return z
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

    z = ThreadsX.map(time_interval) do t
        P * Diagonal(exp.(t * Λ)) * Pinv * initial_state
    end
    return z
end
