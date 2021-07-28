"""
    get_interaction_matrix(problem)
"""
function get_interaction_matrix(problem)
    @debug "start: get_interaction_matrix"
    
    G = get_empty_matrix(problem.physic, problem.atoms)
    problem.kernelFunction(problem.atoms, problem.laser, G)

    @debug "end  : get_interaction_matrix"
    return G
end

function get_empty_matrix(physic::Scalar, atoms::Shape{<:ThreeD})
    Array{ComplexF64}(undef, atoms.N, atoms.N)
end
function get_empty_matrix(physic::Vectorial, atoms::Shape{<:ThreeD})
    Array{ComplexF64}(undef, 3atoms.N, 3atoms.N)
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
    experimental: green_vectorial!(atoms, laser, G)    
"""
function green_vectorial!(atoms, laser, G)
    @debug "start: green_vectorial!"

    δ(x,y) = Int(==(x,y))

    Xt, Yt, Zt = r[:,1], r[:,2] ,r[:,3]

    Xjn = Xt*ones(1,N) - ones(N,1)*Xt'
    Yjn = Yt*ones(1,N) - ones(N,1)*Yt'
    Zjn = Zt*ones(1,N) - ones(N,1)*Zt'
    
    array_XYZ_jn = [Xjn, Yjn, Zjn]
    Rjn = sqrt.(Xjn.^2 + Zjn.^2 + Yjn.^2)
    
    P(x) =  1 - 1/x + 1/x^2 
    Q(x) = -1 + 3/x - 3/x^2 

    α_range = β_range = [-1, +1, 0]

    term2 = (3/2)*exp.(im*Rjn)./Rjn
    term2[findall(isnan.(term2))] .= 0
    term2[findall(isinf.(term2))] .= 0

    P_Rjn = P.(im*Rjn)
    P_Rjn[findall(isnan.(P_Rjn))] .= 0

    Q_Rjn_over_Rjn_squared = Q.(im*Rjn)./Rjn.^2
    Q_Rjn_over_Rjn_squared[findall(isnan.(Q_Rjn_over_Rjn_squared))] .= 0

    A = []
    for (α_idx, α) ∈ enumerate(α_range)
        B = []
        for (β_idx, β) ∈ enumerate(β_range)
            term1 = im*I(N).*δ(α,β)
            # term2 = (3/2)*exp.(im*Rjn)./Rjn  ## Defined outside      
            term3 = P_Rjn.*δ(α,β)
            term4 = Q_Rjn_over_Rjn_squared.*(array_XYZ_jn[α_idx].*array_XYZ_jn[β_idx])

            K = term1 + term2.*(term3 .+ term4)
            push!(B, K)
        end
        push!(A, vcat(B[1:length(β_range)]...)  )
    end

    G[:] = hcat(A[1:length(α_range)]...)
    
    @debug "end  : green_vectorial!"
    return nothing
end

# """
#     get_interaction_matrix(problem)
# """
# function get_interaction_matrix(problem::SimulationScalar)
#     H = zeros(ComplexF64, problem.atoms.N, problem.atoms.N)
#     problem.KernelFunction!(problem.atoms, problem.laser, H)
#     return H
# end
# """
#     get_interaction_matrix(problem, H)
# """
# function get_interaction_matrix(problem::SimulationScalar, H)    
#     problem.KernelFunction!(problem.atoms, problem.laser, H)
#     return H
# end
# """
#     get_energy_shift_and_linewith(problem::SimulationScalar)
# """
# function get_energy_shift_and_linewith(problem::SimulationScalar)
#     spectrum = get_spectrum(problem)
#     ωₙ, Γₙ = imag.(spectrum.λ), -real.(spectrum.λ)
#     if any(Γₙ .< 0 )
#         @warn "some Γₙ were negatives and were ignored"
#         Γₙ = abs.(Γₙ)
#     end
#     return ωₙ, Γₙ
# end
# """
#     get_ψ²(problem::SimulationScalar, n::Integer)
# """
# function get_ψ²(problem::SimulationScalar, n::Integer)
#     return abs2.(problem.ψ[:, n])
# end


### --------------- Mean Field---------------
# """
#     get_interaction_matrix(problem::SimulationMeanField)
# returns the Scalar Problem matrix
# """
# function get_interaction_matrix(problem::SimulationMeanField)
#     H = zeros(ComplexF64, problem.atoms.N, problem.atoms.N)
#     green_scalar!(problem.atoms, problem.laser, H)
#     return H
# end