"""
    interaction_matrix(::LinearOptics)

The interaction matrix dependends on the dimension and model.

For Scalar/Vectorial in 3D check Eq 2.23 and Eq 2.26 of the Thesis: `https://doi.org/10.11606/T.76.2024.tde-26012024-114225`
"""
function interaction_matrix(problem::LinearOptics)
    @debug "start: interaction_matrix"

    G = get_empty_matrix(problem.physic, problem.atoms)
    problem.kernelFunction!(problem.atoms, problem.laser, G)

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

    temp_scalar_problem = LinearOptics(Scalar(), problem.atoms, problem.laser)
    G = get_empty_matrix(temp_scalar_problem.physic, temp_scalar_problem.atoms)
    temp_scalar_problem.kernelFunction!(temp_scalar_problem.atoms, problem.laser, G)

    temp_scalar_problem = 1
    @debug "end  : interaction_matrix"
    return G
end
function interaction_matrix(problem::NonLinearOptics{PairCorrelation})
	@debug "start: interaction_matrix"

	temp_scalar_problem = LinearOptics(Scalar(), problem.atoms, problem.laser)
	G = get_empty_matrix(temp_scalar_problem.physic, temp_scalar_problem.atoms)
	temp_scalar_problem.kernelFunction!(temp_scalar_problem.atoms, problem.laser, G)


	temp_scalar_problem = 1

	@debug "end  : interaction_matrix"
	return G
end

"""
    green_scalar!(atoms, laser, G)

Computes:
    @. G = +1im*(Γ/2)*exp(1im*k₀ * R_jk) / (k₀ * R_jk)
    G[diagind(G)] .= 1im * laser.Δ - Γ/2
"""
function green_scalar!(atoms, laser, G)
    @debug "start: green_scalar!"

    N = atoms.N
    if N == 1
        fill!(G, im*laser.Δ - Γ/2)
        return nothing
    end

    get_pairwise_matrix!(atoms.r, G) # R_jm
    G .*= k₀

    Threads.@threads for j in eachindex(G)
        @inbounds G[j] = +im*(Γ / 2) * cis(G[j]) / (G[j])
    end
    G[diagind(G)] .= 1im * laser.Δ - Γ / 2

    @debug "end  : green_scalar!"
    return nothing
end

function Rjm_parallel!(Rjm, Xjm, Yjm, Zjm)
    Threads.@threads for j in eachindex(Rjm)
        Rjm[j] = sqrt(Xjm[j]^ 2 + Zjm[j]^ 2 + Yjm[j]^ 2)
    end
    return nothing
end
function temp1_parallel!(temp1, Rjm)
    Threads.@threads for j in eachindex(temp1)
        temp1[j] = im*(Γ/2)*(3cis(k₀*Rjm[j]))/(2*k₀*Rjm[j])
    end
    return nothing
end
function temp2_parallel!(temp2, Rjm)
    Threads.@threads for j in eachindex(temp2)
        temp2[j] = ( im/(k₀*Rjm[j]) - 1.0/(k₀*Rjm[j])^2)
    end
    return nothing
end

function Gx_parallel!(Gxx, Gyx, Gzx, temp1, Xjm, Yjm, Zjm, temp2)
    Threads.@threads for j in eachindex(Gxx)
        _temp1 = temp1[j]
        _temp2 = temp2[j]
        _Xjm = Xjm[j]
        _Yjm = Yjm[j]
        _Zjm = Zjm[j]
        Gxx[j] = _temp1*( (1.0 - _Xjm*_Xjm) + (1.0 - 3.0*_Xjm*_Xjm)*_temp2)
        Gyx[j] = _temp1*( ( - _Yjm*_Xjm) + ( - 3.0*_Yjm*_Xjm)*_temp2)
        Gzx[j] = _temp1*( ( - _Zjm*_Xjm) + ( - 3.0*_Zjm*_Xjm)*_temp2)
    end
    return nothing
end

function Gy_parallel!(Gxy, Gyy, Gzy, temp1, Xjm, Yjm, Zjm, temp2)
    Threads.@threads for j in eachindex(Gxy)
        _temp1 = temp1[j]
        _temp2 = temp2[j]
        _Xjm = Xjm[j]
        _Yjm = Yjm[j]
        _Zjm = Zjm[j]
        Gxy[j] = _temp1*( ( - _Xjm*_Yjm) + ( - 3.0*_Xjm*_Yjm)*_temp2)
        Gyy[j] = _temp1*( (1.0 - _Yjm*_Yjm) + (1.0 - 3.0*_Yjm*_Yjm)*_temp2)
        Gzy[j] = _temp1*( ( - _Zjm*_Yjm) + ( - 3.0*_Zjm*_Yjm)*_temp2)
    end
    return nothing
end

function Gz_parallel!(Gxz, Gyz, Gzz, temp1, Xjm, Yjm, Zjm, temp2)
    Threads.@threads for j in eachindex(Gxz)
        _temp1 = temp1[j]
        _temp2 = temp2[j]
        _Xjm = Xjm[j]
        _Yjm = Yjm[j]
        _Zjm = Zjm[j]
        Gxz[j] = _temp1*( ( - _Xjm*_Zjm) + ( - 3.0*_Xjm*_Zjm)*_temp2)
        Gyz[j] = _temp1*( ( - _Yjm*_Zjm) + ( - 3.0*_Yjm*_Zjm)*_temp2)
        Gzz[j] = _temp1*( (1.0 - _Zjm*_Zjm) + (1.0 - 3.0*_Zjm*_Zjm)*_temp2)
    end
     return nothing
end


function G_removeNaN!(G, Δ)
    Threads.@threads for j in eachindex(G)
        if isnan(G[j])
            G[j] = 0.0im
        end
    end
    return nothing
end

function my_pairwise(a, b=a)
    na = length(a)
    r = Array{eltype(a)}(undef, na, na)

    @inbounds for (j, bj) in enumerate(b), (i, ai) in enumerate(a)
        r[i, j] = ai - bj
    end
    return r
end
function green_vectorial!(atoms, laser, G)
    @debug "start: green_vectorial!"

    N = atoms.N
    Δ = laser.Δ
    if N==1
        fill!(G, im*Δ - Γ/2)
        return nothing
    end

    Xt, Yt, Zt = atoms.r[1, :], atoms.r[2, :], atoms.r[3, :]
    Xjm = my_pairwise(Xt)
    Yjm = my_pairwise(Yt)
    Zjm = my_pairwise(Zt)

    Rjm = Array{Float64,2}(undef, N, N)
    Rjm_parallel!(Rjm, Xjm, Yjm, Zjm)

    Xjm = view(Xjm./Rjm, :, :)
    Yjm = view(Yjm./Rjm, :, :)
    Zjm = view(Zjm./Rjm, :, :)


    temp1 = Array{ComplexF64,2}(undef, N, N)
    temp2 = Array{ComplexF64,2}(undef, N, N)
    temp1_parallel!(temp1, Rjm)
    temp2_parallel!(temp2, Rjm)


    ## fill matriz by collumns, because Julia matrices are column-major
    ## G[:, 1:N] = [Gxx; Gyx; Gzx]
    Gxx, Gyx, Gzx = Array{ComplexF64,2}(undef, N, N), Array{ComplexF64,2}(undef, N, N), Array{ComplexF64,2}(undef, N, N)
    Gx_parallel!(Gxx, Gyx, Gzx, temp1, Xjm, Yjm, Zjm, temp2)

    @inbounds @. G[1:N,         1:N] = copy(Gxx)
    @inbounds @. G[(N+1):(2N),  1:N] = copy(Gyx)
    @inbounds @. G[(2N+1):(3N), 1:N] = copy(Gzx)



    ## G[:, (N+1):(2N)] = [Gxy; Gyy; Gzy]
    Gxy, Gyy, Gzy = Gxx, Gyx, Gzx
    Gy_parallel!(Gxx, Gyx, Gzy, temp1, Xjm, Yjm, Zjm, temp2)

    @inbounds @. G[1:N,         (N+1):(2N)] = copy(Gxy)
    @inbounds @. G[(N+1):(2N),  (N+1):(2N)] = copy(Gyy)
    @inbounds @. G[(2N+1):(3N), (N+1):(2N)] = copy(Gzy)


    ## G[:, (2N+1):(3N)] = [Gxz; Gyz; Gzz]
    Gxz, Gyz, Gzz = Gxx, Gyx, Gzx
    Gz_parallel!(Gxz, Gyz, Gzz, temp1, Xjm, Yjm, Zjm, temp2)

    @inbounds @. G[1:N,        (2N+1):(3N)]  = copy(Gxz)
    @inbounds @. G[(N+1):(2N),  (2N+1):(3N)] = copy(Gyz)
    @inbounds @. G[(2N+1):(3N), (2N+1):(3N)] = copy(Gzz)

    # I don't know how to justify that this works
    # 'Γ' vectorial is different than 'Γ' scalar (Γ_vec = 2/3 Γ_sca)
    G[diagind(G)] .= im*Δ - Γ/2 # diagonals have the single atom solution
    G_removeNaN!(G, Δ)

    # force clean variables
    Xjm = Yjm = Zjm = temp1 = temp2 = 1

    @debug "end  : green_vectorial!"
    return nothing
end