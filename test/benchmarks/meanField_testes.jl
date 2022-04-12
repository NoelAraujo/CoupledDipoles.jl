"""
    G DOES has "-Γ/2" inside
"""
function dev_MeanField!(du, u, p, t)
    G, diagG, Ωₙ, Wₙ, Δ, N, G_βₙ, temp1, temp2 = p

    βₙ = @view u[1:N]
    zₙ = @view u[N+1:end]

    mul!(G_βₙ, G, βₙ) # == G_βₙ = G*βₙ
    Wₙ .= Ωₙ / 2 .+ im*G_βₙ
    temp1 .= (im * Δ - Γ / 2) * βₙ .+ im * Wₙ .* zₙ
    temp2 .= -Γ * (1 .+ zₙ) - 4 * imag.(βₙ .* conj.(Wₙ))
    du[:] .= vcat(temp1, temp2)


    return nothing
end
"""
    G does NOT has "-Γ/2" inside
"""
function dev_MeanField!_v2(du, u, p, t)
    G, diagG, Ωₙ, Wₙ, Δ, N, G_βₙ, temp1, temp2 = p

    βₙ = @view u[1:N]
    zₙ = @view u[N+1:end]

    for j=1:N
        temp1[j] = (im*Δ - Γ/2)*βₙ[j] + 0.5*im*Ωₙ[j]*zₙ[j]  + (Γ/2)*sum( G[j,m]*βₙ[m] for m =1:N if j ≠ m  )*zₙ[j]
    end
    for j=1:N
        temp2[j] = (-im*Ωₙ[j]*conj(βₙ[j]) + im*conj(Ωₙ[j])*βₙ[j])  - Γ*(1 + zₙ[j])  - (sum(G[j,m]*βₙ[m]*conj(βₙ[j]) for m =1:N if m≠j) +  conj.( sum(G[j,m]*βₙ[m]*conj(βₙ[j]) for m =1:N if m≠j)  )  )
    end 
    
    du[:] .= vcat(temp1, temp2)

    return nothing
end
"""
    G DOES has "-Γ/2" inside
"""
function dev_MeanField!_v3(du, u, p, t)
    G, diagG, Ωₙ, Wₙ, Δ, N, G_βₙ, temp1, temp2 = p

    βₙ = @view u[1:N]
    zₙ = @view u[N+1:end]

    for j=1:N
        temp1[j] = (im*Δ - Γ/2)*βₙ[j] + 0.5*im*Ωₙ[j]*zₙ[j]  + (-2/Γ)*(Γ/2)*sum( G[j,m]*βₙ[m] for m =1:N if j ≠ m  )*zₙ[j]
    end
    for j=1:N
        temp2[j] = (-im*Ωₙ[j]*conj(βₙ[j]) + im*conj(Ωₙ[j])*βₙ[j])  - Γ*(1 + zₙ[j])  - (-2/Γ)*(sum(G[j,m]*βₙ[m]*conj(βₙ[j]) for m =1:N if m≠j) +  conj.( sum(G[j,m]*βₙ[m]*conj(βₙ[j]) for m =1:N if m≠j)  )  )
    end 
    
    du[:] .= vcat(temp1, temp2)
    return nothing
end


using LinearAlgebra
const Γ = 1

N = 5
rjk = rand(N,N)
G1 = -(Γ/2)*(cos.(rjk)./rjk + im*(sin.(rjk)./rjk))
G2 = (cos.(rjk)./rjk + im*(sin.(rjk)./rjk))
G3 = -(Γ/2)*(cos.(rjk)./rjk + im*(sin.(rjk)./rjk))
G1[diagind(G1)] .= 0


Ωₙ = rand(ComplexF64, N)
Wₙ = zeros(ComplexF64, N)
Δ = rand()
G_βₙ = similar(Ωₙ)
temp1 = similar(Ωₙ)
temp2 = similar(Ωₙ)

p1 = G1, diag(G1), Ωₙ, Wₙ, Δ, N, G_βₙ, temp1, temp2
p2 = G2, diag(G2), Ωₙ, Wₙ, Δ, N, G_βₙ, temp1, temp2
p3 = G3, diag(G3), Ωₙ, Wₙ, Δ, N, G_βₙ, temp1, temp2

du1 = zeros(ComplexF64, 2N)
du2 = zeros(ComplexF64, 2N)
du3 = zeros(ComplexF64, 2N)

u = rand(ComplexF64, 2N)
t = 0
dev_MeanField!(du1, u, p1, t)
dev_MeanField!_v2(du2, u, p2, t)
dev_MeanField!_v3(du3, u, p3, t)

[du1 du2 du3 du1.≈du2 du1.≈du3] # |> vscodedisplay