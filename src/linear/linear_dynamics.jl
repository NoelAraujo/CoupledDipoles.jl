"""
    default_initial_condition(::Scalar) = zeros(ComplexF64, N)
"""
function default_initial_condition(problem::LinearOptics{Scalar})
    return zeros(ComplexF64, problem.atoms.N) # I must use "zeros" and NOT an undef array - with trash data inside
end

"""
    steady_state(problem)
"""
function steady_state(problem::LinearOptics{Scalar})
    @debug "start: get steady state"

    G = interaction_matrix(problem)
    Ωₙ = laser_field(problem.laser, problem.atoms)
    if problem.atoms.N > 1
        βₛ = -(G \ Ωₙ)
    else
        # The negative sign was NOT forgotten
        # after some math, I verified that it does not exist
        βₛ = Ωₙ / G[1]
    end

    @debug "end  : get steady state"
    return βₛ
end

function time_evolution(problem::LinearOptics{T}, u₀, tspan::Tuple; kargs...) where {T<:Linear}
    ### use default G and Ωₙ
    G = copy(interaction_matrix(problem))
    Ωₙ = laser_field(problem.laser, problem.atoms)
    solution = time_evolution(problem, u₀, tspan, Ωₙ, G; kargs...)

    return solution
end
function time_evolution(problem::LinearOptics{T}, u₀, tspan::Tuple, Ωₙ::Vector; kargs...) where {T<:Linear}
    ### use default G
    G = copy(interaction_matrix(problem))
    solution = time_evolution(problem, u₀, tspan, Ωₙ, G; kargs...)

    return solution
end
function time_evolution(problem::LinearOptics{T}, u₀, tspan::Tuple, Ωₙ::Vector, G::Matrix; kargs...) where {T<:Linear}
    @debug "start: time evolution - LinearOptics"

    ### parameters == constant vector and matrices
    parameters = view(G, :, :), view(Ωₙ, :)

    ### calls for solver
    problemFunction = get_evolution_function(problem)
    prob = ODEProblem(problemFunction, u₀, tspan, parameters)
    solution = DifferentialEquations.solve(prob, VCABM(); dt=1e-10, abstol=1e-10, reltol=1e-10, kargs...)

    @debug "end  : time evolution - LinearOptics"
    return solution
end

get_evolution_function(problem::LinearOptics{Scalar}) = Scalar!

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
