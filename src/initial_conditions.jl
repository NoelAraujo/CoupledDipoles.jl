"""
    default_initial_condition(::Scalar)

    u₀ = [ <σ⁻> ] = [ zeros ]
"""
function default_initial_condition(problem::LinearOptics{Scalar})
    return zeros(ComplexF64, problem.atoms.N) # I must use "zeros" and NOT an undef array - with trash data inside
end


"""
    default_initial_condition(::Vectorial)

    u₀ = [ <σ⁻>_x | <σ⁻>_y | <σ⁻>_z ] = [ zeros | zeros | zeros ]

"""
function default_initial_condition(problem::LinearOptics{Vectorial})
    return zeros(ComplexF64, 3problem.atoms.N) # I must use "zeros" and NOT an undef array - with trash data inside
end


"""
    default_evolution_initial_condition(NonLinearOptics{MeanField})
β₀ = zeros(ComplexF64, atoms.N)
z₀ = 2β₀.*conj.(β₀) .- 1

u₀ = [ <σ⁻> | <σᶻ> ] = [ zeros | -ones ]
"""
function default_initial_condition(problem::NonLinearOptics{MeanField})
    β₀ = zeros(ComplexF64, problem.atoms.N)
    z₀ = -ones(ComplexF64, problem.atoms.N)
    u₀ = vcat(β₀, z₀)
    return u₀
end

"""
    default_evolution_initial_condition(NonLinearOptics{PairCorrelation})

    u₀ = [ <σ⁻> | <σᶻ> | <σᶻσ⁻> | <σ⁺σ⁻> | <σ⁻σ⁻> | <σᶻσᶻ> ] = [ zeros | -ones | zeros | zeros | zeros | ones (zero on diagonal) ]
"""
function default_initial_condition(problem::NonLinearOptics{PairCorrelation})
    N = problem.atoms.N
    u₀ = zeros(ComplexF64, 2*N + 4*N^2)

    # i can't argument/explain/justify the conditions
    u₀[N+1:2*N] .= -1 # σᶻ
    # only 'σᶻσᶻ' has different values than zero
    # it is a matrix with ones, and zeros at diagonal
    u₀[2*N+3*N^2+1:2*N+4*N^2] .= +1
    u₀[2*N+3*N^2+1:N+1:2*N+4*N^2] .= 0.0

    return u₀
end