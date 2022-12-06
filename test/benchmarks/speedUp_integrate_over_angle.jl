@views function v0(N, r, atoms_states, θ)
    β = view(atoms_states, 1:N)
    βp = conj.(β)

    x, y, z = r[1,:], r[2,:], r[3,:]
    intensity = zero(ComplexF64)
    for j=1:N
        for jp=1:N
            xjj = x[j] - x[jp]
            yjj = y[j] - y[jp]
            zjj = z[j] - z[jp]

            intensity += myKernel2(β[j], βp[jp], xjj, yjj, zjj, θ)
        end
    end
    return real(intensity)
end
function myKernel2(beta_j, beta_jp, xjj, yjj, zjj, θ)
    k₀sinθ = abs(k₀*sin(θ))
    k₀cosθ = k₀*cos(θ)

    v = beta_j*beta_jp*exp(-im*zjj*k₀cosθ)*Bessels.besselj0(k₀sinθ*sqrt(xjj^2+yjj^2))
    return v
end

@views function v0_2(N, r, atoms_states, θ)

    βₙ = view(atoms_states, 1:N)
    βₘ = conj.(βₙ)

    k₀sinθ = abs(k₀*sin(θ))
    k₀cosθ = k₀*cos(θ)
    number_configurations = ((N^2) ÷ 2 - N ÷ 2)

    xₙₘ = Array{Float64}(undef, number_configurations)
    yₙₘ, zₙₘ = similar(xₙₘ), similar(xₙₘ)
    count = 1
    for m in 1:N
        rm = r[:,m]
        for n in (m+1):N
            xₙₘ[count] = r[1,n] - rm[1]
            yₙₘ[count] = r[2,n] - rm[2]
            zₙₘ[count] = r[3,n] - rm[3]
            count += 1
        end
    end

    βₙₘ = Array{ComplexF64}(undef, number_configurations)
    count = 1
    for m in 1:N
        b_m = βₘ[m]
        for n in (m+1):N
            b_n = βₙ[n]
            βₙₘ[count] = b_n*b_m
            count += 1
        end
    end

    total_intensity = ThreadsX.mapreduce(+, 1:number_configurations) do ii
        myKernel3(βₙₘ[ii], xₙₘ[ii], yₙₘ[ii], zₙₘ[ii], k₀sinθ, k₀cosθ)
    end
    total_intensity += sum(abs2, βₙ)/2
    return 2real(total_intensity)
end
function myKernel3(beta_jj, xjj, yjj, zjj, k₀sinθ, k₀cosθ)
    v = beta_jj*exp(-im*zjj*k₀cosθ)*Bessels.besselj0(k₀sinθ*sqrt(xjj^2+yjj^2))
    return v
end

using CoupledDipoles
using LinearAlgebra, ThreadsX, SpecialFunctions
using Random; Random.seed!(33)
const k₀ = 1

N = 200; kR = 3
atoms = Atom(Sphere(gaussian=true), N, kR)
r = atoms.r
s, Δ = 1e-6, 0.5
laser = Laser(PlaneWave3D(), s, Δ)

problem = LinearOptics(Scalar(), atoms, laser)
atoms_states = steady_state(problem)

θ = deg2rad(1)

sol1 = v0(N, r, atoms_states, θ);
sol2 = v0_2(N, r, atoms_states, θ);

@show sol1, sol2


sol1, sol2 = [], []
θ_range = range(4π/10, 16π/10, length=100)
for θ in θ_range
    push!(sol1, v0(N, r, atoms_states, θ))
    push!(sol2, v0_2(N, r, atoms_states, θ))
end

using Plots
plot(θ_range, sol1, lw=3)
plot!(θ_range, sol2, lw=3, linestyle=:dash)


t12 = []
for i=1:1000
    t1 =@elapsed v0(N, r, atoms_states, θ);
    rand(1000);
    t2 =@elapsed v0_2(N, r, atoms_states, θ);
    rand(1000);
    push!(t12, t1/t2)
end

histogram(t12)