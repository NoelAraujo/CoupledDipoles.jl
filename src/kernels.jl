### --------------- Scalar---------------
"""
    green_scalar!(atoms, laser, G)

Computes:  
    @. G = -(Γ/2)*exp(1im*k₀ * R_jk) / (1im * k₀ * R_jk)  
    G[diagind(G)] .= 1im * laser.Δ - Γ/2
"""
function green_scalar!(atoms, laser, G)
    G[:] = get_pairwise_matrix(atoms.r)

    # Threads.@threads for j in eachindex(G)
    #     @inbounds G[j] = -(Γ / 2) * cis(k₀ * G[j]) / (1im * k₀ * G[j])
    # end
    LazyArrays.@~ G .= -(Γ / 2) * cis.(k₀ * G) ./ (1im * k₀ * G);
    G[diagind(G)] .= 1im * laser.Δ - Γ / 2

    return nothing
end
"""
    get_interaction_matrix(problem)
"""
function get_interaction_matrix(problem::SimulationScalar)
    H = zeros(ComplexF64, problem.atoms.N, problem.atoms.N)
    problem.KernelFunction!(problem.atoms, problem.laser, H)
    return H
end
"""
    get_interaction_matrix(problem, H)
"""
function get_interaction_matrix(problem::SimulationScalar, H)    
    problem.KernelFunction!(problem.atoms, problem.laser, H)
    return H
end
"""
    get_energy_shift_and_linewith(problem::SimulationScalar)
"""
function get_energy_shift_and_linewith(problem::SimulationScalar)
    spectrum = get_spectrum(problem)
    ωₙ, Γₙ = imag.(spectrum.λ), -real.(spectrum.λ)
    if any(Γₙ .< 0 )
        @warn "some Γₙ were negatives and were ignored"
        Γₙ = abs.(Γₙ)
    end
    return ωₙ, Γₙ
end
"""
    get_ψ²(problem::SimulationScalar, n::Integer)
"""
function get_ψ²(problem::SimulationScalar, n::Integer)
    return abs2.(problem.ψ[:, n])
end


### --------------- Mean Field---------------
"""
    get_interaction_matrix(problem::SimulationMeanField)
returns the Scalar Problem matrix
"""
function get_interaction_matrix(problem::SimulationMeanField)
    H = zeros(ComplexF64, problem.atoms.N, problem.atoms.N)
    green_scalar!(problem.atoms, problem.laser, H)
    return H
end