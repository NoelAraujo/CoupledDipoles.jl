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

@parallel_indices (x,y) function Rjm_parallel!(Rjm, Xjm, Yjm, Zjm)
        Rjm[x,y] = sqrt(Xjm[x,y]^ 2 + Zjm[x,y]^ 2 + Yjm[x,y]^ 2)
    return nothing
end
@parallel_indices (x,y) function temp1_parallel!(temp1, Rjm)
        temp1[x,y] = (3cis(k₀*Rjm[x,y]))/(2im*k₀*Rjm[x,y])
    return nothing
end
@parallel_indices (x,y) function temp2_parallel!(temp2, Rjm)
        temp2[x,y] = ( im/(k₀*Rjm[x,y]) - 1.0/(k₀*Rjm[x,y])^2)
    return nothing
end

@parallel_indices (x,y) function Gx_parallel!(Gxx, Gyx, Gzx, temp1, Xjm, Yjm, Zjm, temp2)
       Gxx[x,y] = temp1[x,y]*( (1 - Xjm[x,y]*Xjm[x,y]) + (1 - 3.0*Xjm[x,y]*Xjm[x,y])*temp2[x,y])
       Gyx[x,y] = temp1[x,y]*( ( - Yjm[x,y]*Xjm[x,y]) + ( - 3.0*Yjm[x,y]*Xjm[x,y])*temp2[x,y])
       Gzx[x,y] = temp1[x,y]*( ( - Zjm[x,y]*Xjm[x,y]) + ( - 3.0*Zjm[x,y]*Xjm[x,y])*temp2[x,y])
    return nothing
end

@parallel_indices (x,y) function Gy_parallel!(Gxy, Gyy, Gzy, temp1, Xjm, Yjm, Zjm, temp2)
    Gxy[x,y] = temp1[x,y]*( ( - Xjm[x,y]*Yjm[x,y]) + ( - 3.0*Xjm[x,y]*Yjm[x,y])*temp2[x,y])
    Gyy[x,y] = temp1[x,y]*( (1 - Yjm[x,y]*Yjm[x,y]) + (1 - 3.0.*Yjm[x,y]*Yjm[x,y])*temp2[x,y])
    Gzy[x,y] = temp1[x,y]*( ( - Zjm[x,y]*Yjm[x,y]) + ( - 3.0.*Zjm[x,y]*Yjm[x,y])*temp2[x,y])
    return nothing
end

@parallel_indices (x,y) function Gz_parallel!(Gxz, Gyz, Gzz, temp1, Xjm, Yjm, Zjm, temp2)
     Gxz[x,y] = temp1[x,y]*( ( - Xjm[x,y]*Zjm[x,y]) + ( - 3.0*Xjm[x,y]*Zjm[x,y])*temp2[x,y])
     Gyz[x,y] = temp1[x,y]*( ( - Yjm[x,y]*Zjm[x,y]) + ( - 3.0*Yjm[x,y]*Zjm[x,y])*temp2[x,y])
     Gzz[x,y] = temp1[x,y]*( (1 - Zjm[x,y]*Zjm[x,y]) + (1 - 3.0*Zjm[x,y]*Zjm[x,y])*temp2[x,y])
     return nothing
end


@parallel_indices (x,y) function G_diagonals!(G, Δ)
    G[x,y] = -(Γ/2).*G[x,y]
    if isnan(G[x,y])
        G[x,y] = im*Δ - Γ/3
    end
     return nothing
end


function green_vectorial!(atoms, laser, G)
    @debug "start: green_vectorial!"

    N = atoms.N
    Δ = laser.Δ

    Xt, Yt, Zt = atoms.r[1, :], atoms.r[2, :], atoms.r[3, :]
    Xjm = Array{Float64,2}(undef, N, N)

    Xjm = Distances.pairwise(Euclidean(), Xt', Xt'; dims=2)
    Yjm = Distances.pairwise(Euclidean(), Yt', Yt'; dims=2)
    Zjm = Distances.pairwise(Euclidean(), Zt', Zt'; dims=2)

    Rjm = Array{Float64,2}(undef, N, N)
    @parallel Rjm_parallel!(Rjm, Xjm, Yjm, Zjm)

    Xjm = view(Xjm./Rjm, :, :)
    Yjm = view(Yjm./Rjm, :, :)
    Zjm = view(Zjm./Rjm, :, :)


    temp1 = Array{ComplexF64,2}(undef, N, N)
    temp2 = Array{ComplexF64,2}(undef, N, N)
    @parallel temp1_parallel!(temp1, Rjm)
    @parallel temp2_parallel!(temp2, Rjm)


    ## fill matriz by collumns, because Julia matrices are column-major
    ## G[:, 1:N] = [Gxx; Gyx; Gzx]
    Gxx, Gyx, Gzx = Array{ComplexF64,2}(undef, N, N), Array{ComplexF64,2}(undef, N, N), Array{ComplexF64,2}(undef, N, N)
    @parallel Gx_parallel!(Gxx, Gyx, Gzx, temp1, Xjm, Yjm, Zjm, temp2)

    @inbounds @. G[1:N,         1:N] = copy(Gxx)
    @inbounds @. G[(N+1):(2N),  1:N] = copy(Gyx)
    @inbounds @. G[(2N+1):(3N), 1:N] = copy(Gzx)



    ## G[:, (N+1):(2N)] = [Gxy; Gyy; Gzy]
    Gxy, Gyy, Gzy = Gxx, Gyx, Gzx
    @parallel Gy_parallel!(Gxx, Gyx, Gzy, temp1, Xjm, Yjm, Zjm, temp2)

    @inbounds @. G[1:N,         (N+1):(2N)] = copy(Gxy)
    @inbounds @. G[(N+1):(2N),  (N+1):(2N)] = copy(Gyy)
    @inbounds @. G[(2N+1):(3N), (N+1):(2N)] = copy(Gzy)


    ## G[:, (2N+1):(3N)] = [Gxz; Gyz; Gzz]
    Gxz, Gyz, Gzz = Gxx, Gyx, Gzx
    @parallel Gz_parallel!(Gxz, Gyz, Gzz, temp1, Xjm, Yjm, Zjm, temp2)

    @inbounds @. G[1:N,        (2N+1):(3N)]  = copy(Gxz)
    @inbounds @. G[(N+1):(2N),  (2N+1):(3N)] = copy(Gyz)
    @inbounds @. G[(2N+1):(3N), (2N+1):(3N)] = copy(Gzz)


    # DO NOT CHANGE THE ORDER OF THE NEXT TWO LINES
    @parallel G_diagonals!(G, Δ)

    # force clean variables
    Xjm = Yjm = Zjm = temp1 = temp2 = 1

    @debug "end  : green_vectorial!"
    return nothing
end