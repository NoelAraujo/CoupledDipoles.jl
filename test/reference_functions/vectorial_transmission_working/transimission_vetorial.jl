using Random
using LinearAlgebra
using CairoMakie
using ProgressBars

include("kerne3c.jl")
include("Scattering.jl")
include("beam_profile.jl")
mean(x) = sum(x)/length(x)

Base.@kwdef mutable struct Constantes
    lambda =  2*pi
    k = 2*pi/lambda
    Gamma = 1.0
    hbar = 1.0

    epsilon0 = hbar*k^3/(12*pi*Gamma) # consequencia de supor 'ddip/hbar = 1/2'
    ddip = sqrt( (hbar*3*pi)*epsilon0*Gamma/k^3 )
end
c = Constantes()

@assert c.ddip/c.hbar == 0.5 # essa condição é imposta
@assert c.ddip * c.k^3 / (4 * pi * c.epsilon0) == 3/2 # essa condição é consequencia

# ---------------- Simulation parameters ----------------

nRepetitions = 2  # Number of realizations
delta = range(-10, stop=10, length=25)  # Detuning
Delta = delta .* c.Gamma

kR = 17.71830997686217301634315 # 1.3e-6*c.k
R = kR/c.k

b0 = 6
N = round(Int, (kR)^2*b0/3)  # Number of particles

kL = 8.177681527782540982229875
L = kL/c.k #
rholambda3 = N/(kR)^2/kL*(2*pi)^(3/2)

# Observation grid = single point on the (0,0,Ro) direction
Phi = 0
ThObs = 0*pi
Ntheta = length(ThObs)

Ro = 100*kR  # Distance of observation
Xobs = Ro .* sin.(ThObs)' .* cos.(Phi)
Yobs = Ro .* sin.(ThObs)' .* sin.(Phi)
Zobs = Ro .* cos.(ThObs)' .* ones(1, length(Phi))

# Incident beam
E0 = 2000
w0 = R/3  # Beam waist
z0 = pi*w0^2/c.lambda  # Rayleigh length
Elaser(x, y, z) = laser_beam_profile(x, y, z, E0, w0, z0, c.k)

# Initialize arrays to store results
Esx = zeros(ComplexF64, nRepetitions, length(delta))
Esy = zeros(ComplexF64, nRepetitions, length(delta))
Esz = zeros(ComplexF64, nRepetitions, length(delta))
Imean = zeros(length(delta))

function compute_transmission!(delta, nRepetitions, Delta, Esx, Esy, Esz, Imean, c, R, L)
    Gamma, k, hbar, epsilon0, ddip = c.Gamma,  c.k, c.hbar, c.epsilon0, c.ddip

    identity_3N = Diagonal(ones(3*N))
    M3 = zeros(ComplexF64, 3N, 3N)
    M31 = similar(M3)
    c1 = zeros(ComplexF64, N, N)
    c2 = similar(c1)

    ## conver this into a normal loop if needed
    for ib in ProgressBars.ProgressBar(1:length(delta))
        for j in 1:nRepetitions
            Random.seed!(j)
            X = randn(N) .* R
            Y = randn(N) .* R
            Z = randn(N) .* L

            # Calculate scattered light intensity and transmission coefficient
            kernel3c!(X, Y, Z, M31, c1, c2, k);
            @. M3 = (-Gamma/2 + im*Delta[ib])*identity_3N - Gamma*3/4*M31

            Om = (im*ddip*[Elaser.(X, Y, Z); zeros(N); zeros(N)]/hbar)
            # Om = (-im*[Elaser.(X, Y, Z); zeros(N); zeros(N)]/2)
            Beta30 = M3\Om

            _Esx, _Esy, _Esz = Escattered3(X, Y, Z, Beta30, Xobs, Yobs, Zobs, k)
            Esx[j, ib], Esy[j, ib], Esz[j, ib] = _Esx[1], _Esy[1], _Esz[1]
        end

        # Calculate mean intensity
        factor = im * ddip * k^3 / (4 * pi * epsilon0)
        # factor = -im*3Gamma/2
        I = abs.(factor .* Esx[:, ib] .+ Elaser.(Xobs, Yobs, Zobs)).^2 .+
            abs.(factor .* Esy[:, ib]).^2 .+
            abs.(factor .* Esz[:, ib]).^2

        Imean[ib] = mean(I) / abs(Elaser.(Xobs, Yobs, Zobs)[1]).^2
    end
end

let
    @time compute_transmission!(delta, nRepetitions, Delta, Esx, Esy, Esz, Imean, c, R, L);
    # @profview compute_transmission!(delta, nRepetitions, Delta, Esx, Esy, Esz, Imean, c, R, L)

    fig = Figure()
    # ax = Axis(fig[1, 1], xlabel="Detuning", ylabel="Intensity", yscale=log10)
    ax = Axis(fig[1, 1], xlabel="Detuning", ylabel="Intensity")
    lines!(ax, delta, exp.(-b0 ./ (1 .+ 4 .* delta.^2)))
    lines!(ax, collect(delta), Imean)
    # xlims!(ax, -5, 5)
    fig
end