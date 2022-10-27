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
    return Array{ComplexF64}(undef, atoms.N, atoms.N)
end
function get_empty_matrix(physic::Vectorial, atoms::Atom{<:ThreeD})
    return Array{ComplexF64}(undef, 3atoms.N, 3atoms.N)
end

function interaction_matrix(problem::NonLinearOptics{MeanField})
    @debug "start: interaction_matrix"

    # if haskey(problem.data, :G)
    #     return problem.data[:G]
    # end

    temp_scalar_problem = LinearOptics(Scalar(), problem.atoms, problem.laser)
    G = get_empty_matrix(temp_scalar_problem.physic, temp_scalar_problem.atoms)
    temp_scalar_problem.kernelFunction(temp_scalar_problem.atoms, problem.laser, G)

    temp_scalar_problem = 1
    # problem.data[:G] = G
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

    get_pairwise_matrix!(atoms.r, G) # R_jk

    Threads.@threads for j in eachindex(G)
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

    N = atoms.N
    Δ = laser.Δ

    Xt, Yt, Zt = atoms.r[1, :], atoms.r[2, :], atoms.r[3, :]

    Xjm = Xt * ones(1, N) - ones(N, 1) * Xt'
    Yjm = Yt * ones(1, N) - ones(N, 1) * Yt'
    Zjm = Zt * ones(1, N) - ones(N, 1) * Zt'
    Rjm = sqrt.(Xjm .^ 2 + Zjm .^ 2 + Yjm .^ 2)

    Xjm = view(Xjm./Rjm, :, :)
    Yjm = view(Yjm./Rjm, :, :)
    Zjm = view(Zjm./Rjm, :, :)

    temp1 = (3cis.(k₀*Rjm))./(2im*k₀.*Rjm)
    temp2 = ( im./(k₀.*Rjm) - 1.0./(k₀.*Rjm).^2)
    onesTemp = fill(one(eltype(G)) ,N,N)

    ## fill matriz by collumns, because Julia matrices are column-major
    ## G[:, 1:N] = [Gxx; Gyx; Gzx]
    G[1:N,         1:N] .= temp1.*( (onesTemp - Xjm.*Xjm) .+ (onesTemp - 3.0.*Xjm.*Xjm).*temp2)
    G[(N+1):(2N),  1:N] .= temp1.*( ( - Yjm.*Xjm) .+ ( - 3.0.*Yjm.*Xjm).*temp2)
    G[(2N+1):(3N), 1:N] .= temp1.*( ( - Zjm.*Xjm) .+ ( - 3.0.*Zjm.*Xjm).*temp2)

    ## G[:, (N+1):(2N)] = [Gxy; Gyy; Gzy]
    G[1:N,         (N+1):(2N)] .= temp1.*( ( - Xjm.*Yjm) .+ ( - 3.0.*Xjm.*Yjm).*temp2)
    G[(N+1):(2N),  (N+1):(2N)] .= temp1.*( (onesTemp - Yjm.*Yjm) .+ (onesTemp - 3.0.*Yjm.*Yjm).*temp2)
    G[(2N+1):(3N), (N+1):(2N)] .= temp1.*( ( - Zjm.*Yjm) .+ ( - 3.0.*Zjm.*Yjm).*temp2)

    ## G[:, (2N+1):(3N)] = [Gxz; Gyz; Gzz]
    G[1:N,        (2N+1):(3N)]  .= temp1.*( ( - Xjm.*Zjm) .+ ( - 3.0.*Xjm.*Zjm).*temp2)
    G[(N+1):(2N),  (2N+1):(3N)] .= temp1.*( ( - Yjm.*Zjm) .+ ( - 3.0.*Yjm.*Zjm).*temp2)
    G[(2N+1):(3N), (2N+1):(3N)] .= temp1.*( (onesTemp - Zjm.*Zjm) .+ (onesTemp - 3.0.*Zjm.*Zjm).*temp2)

    # DO NOT CHANGE THE ORDER OF THE NEXT TWO LINES
    G .= -(Γ/2).*G
    G[findall( isnan.(G) )] .= im*Δ - Γ/3

    # force clean variables
    Xjm = Yjm = Zjm = temp1 = temp2 = onesTemp = 1

    @debug "end  : green_vectorial!"
    return nothing
end