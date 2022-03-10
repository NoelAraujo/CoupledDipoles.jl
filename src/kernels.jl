"""
    interaction_matrix(::LinearOptics)
"""
function interaction_matrix(problem::LinearOptics)
    @debug "start: interaction_matrix"

    G = get_empty_matrix(problem.physic, problem.atoms)
    problem.kernelFunction(problem.atoms, problem.laser, G)

    @debug "end  : interaction_matrix"
    return G
end

function get_empty_matrix(physic::Scalar, atoms::Atom{<:ThreeD})
    Array{ComplexF64}(undef, atoms.N, atoms.N)
end
function get_empty_matrix(physic::Vectorial, atoms::Atom{<:ThreeD})
    Array{ComplexF64}(undef, 3atoms.N, 3atoms.N)
end

function interaction_matrix(problem::NonLinearOptics{MeanField})
    @debug "start: interaction_matrix"

    if haskey(problem.data, :G)
        return problem.data[:G]
    end

    temp_scalar_problem = LinearOptics(Scalar(), problem.atoms, problem.laser)
    G = get_empty_matrix(temp_scalar_problem.physic, temp_scalar_problem.atoms)
    temp_scalar_problem.kernelFunction(temp_scalar_problem.atoms, problem.laser, G)

    temp_scalar_problem = 1
    GC.gc()
    problem.data[:G] = G
    @debug "end  : interaction_matrix"
    return G
end



"""
    green_scalar!(atoms, laser, G)

Computes:  
    @. G = -(Γ/2)*exp(1im*k₀ * R_jk) / (1im * k₀ * R_jk)  
    G[diagind(G)] .= 1im * laser.Δ - Γ/2
"""
function green_scalar!(atoms, laser, G)
    @debug "start: green_scalar!"

    G[:] = get_pairwise_matrix(atoms.r) # R_jk

    Threads.@threads for j in eachindex(G) # 
        @inbounds G[j] = -(Γ / 2) * cis(k₀ * G[j]) / (1im * k₀ * G[j])
    end
    G[diagind(G)] .= 1im * laser.Δ - Γ / 2

    @debug "end  : green_scalar!"
    return nothing
end


"""
    green_vectorial!(atoms, laser, G)

check bencmark file for details
"""
function green_vectorial!(atoms, laser, G)
    @debug "start: green_vectorial!"

    δ(x, y) = Int(==(x, y))
    N = atoms.N
    Δ = laser.Δ

    Xt, Yt, Zt = atoms.r[1, :], atoms.r[2, :], atoms.r[3, :]

    Xjn = Xt * ones(1, N) - ones(N, 1) * Xt'
    Yjn = Yt * ones(1, N) - ones(N, 1) * Yt'
    Zjn = Zt * ones(1, N) - ones(N, 1) * Zt'

    array_XYZ_jn = [Xjn, Yjn, Zjn]
    Rjn = sqrt.(Xjn .^ 2 + Zjn .^ 2 + Yjn .^ 2)

    P(x) = 1 - 1 / x + 1 / x^2
    Q(x) = -1 + 3 / x - 3 / x^2

    α_range = β_range = [-1, +1, 0]

    term2 = (3 / 2) * exp.(im * k₀ * Rjn) ./ (k₀ * Rjn)
    term2[findall(isnan.(term2))] .= 0
    term2[findall(isinf.(term2))] .= 0

    P_Rjn = P.(im * k₀ * Rjn)
    P_Rjn[findall(isnan.(P_Rjn))] .= 0

    Q_Rjn_over_Rjn_squared = Q.(im * k₀ * Rjn) ./ (k₀ * Rjn) .^ 2
    Q_Rjn_over_Rjn_squared[findall(isnan.(Q_Rjn_over_Rjn_squared))] .= 0

    A = []
    for (α_idx, α) ∈ enumerate(α_range)
        B = []
        for (β_idx, β) ∈ enumerate(β_range)
            term1 = im * I(N) .* δ(α, β)
            # term2 = (3/2)*exp.(im*Rjn)./Rjn  ## Defined outside      
            term3 = P_Rjn .* δ(α, β)
            term4 = Q_Rjn_over_Rjn_squared .* (array_XYZ_jn[α_idx] .* array_XYZ_jn[β_idx])

            K = term1 + term2 .* (term3 .+ term4)
            push!(B, K)
        end
        push!(A, vcat(B[1:length(β_range)]...))
    end

    G[:] = (-Γ / 2) * hcat(A[1:length(α_range)]...) #join all matrices
    G .= -im * G # For my code to work, I need this imaginary number. Credits to Sheila
    G[diagind(G)] .= (-Γ / 2 + 1im * Δ)

    @debug "end  : green_vectorial!"
    return nothing
end
