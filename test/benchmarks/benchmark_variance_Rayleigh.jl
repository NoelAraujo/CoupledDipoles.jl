using CoupledDipole, Revise
using Random
Random.seed!(1134)

using Statistics: mean, var
using StatsBase: fit, normalize, Histogram
using ProgressMeter, Plots, LinearAlgebra

### ------------ ATOMS SPECS ---------------------
ρλ³ = 45
kL = 32.4;
ρ_over_k₀³ = ρλ³ / (2π)^3 
N = floor(Int, ρ_over_k₀³ * kL^3)

### ------------ LASER SPECS ---------------------
Δ = 1
s = 1e-6

### -------- PRODUCE INTENSITIES -----------------
sensors = ring_on_space(; num_pts=360, kR=1.5kL, θ=5π / 12)
many_intensities = Float64[]
maxRep = 10

@showprogress for rep in 1:maxRep
    atoms = Cube(N, kL)
    laser = Gaussian_3D(estimate_waist(atoms), s, Δ)
    simulation = ScalarProblem(atoms, laser)

    βₙ = get_steady_state(simulation)
    intensities = get_scattered_intensity(simulation, βₙ, sensors)

    for i in intensities
        push!(many_intensities, copy(i))
    end
end
all_intensities_over_mean = many_intensities ./ mean(many_intensities)

### ------------ CREATE HISTOGRAM ---------------------
bins = 10.0 .^ range(log10(1e-10), log10(100); length=100)
h = fit(Histogram, all_intensities_over_mean, bins)

h_norm = normalize(h; mode=:pdf)
bins_edges = collect(h_norm.edges[1])
bins_centers = [
    sqrt(bins_edges[i] * bins_edges[i + 1]) for i in 1:(length(bins_edges) - 1)
]


### ------------ FINAL PLOT ---------------------
## theoretical curve
x_ray = range(0.01, 100; step=0.15)
y_ray = exp.(-x_ray)
plot(x_ray, y_ray; linestyle=:dash, label="Rayleigh", lw=4)

scatter!(
    bins_centers, h_norm.weights; label="N=$(N), kL=$(kL)", markershape=:circle, markersize=5
)
plot!(; guidefont=15, tickfont=15, legendfont=10, size=(1000, 500))
plot!(;
    left_margin=5Plots.mm,
    right_margin=5Plots.mm,
    top_margin=5Plots.mm,
    bottom_margin=5Plots.mm,
    c=:blue,
    gridalpha=0.8,
    minorgrid=true,
    minorgridalpha=0.5,
    minorgridstyle=:dashdot,
    yscale=:log10,
    xscale=:log10,
    xlims=(10^-1, 10^2),
    xticks=[10^-1, 10^0, 10^1, 10^2],
    ylims=(10^-6, 10^1),
    yticks=[10^1, 10^0, 10^-2, 10^-4, 10^-6],
)

plot!(; ylims=(1e-6, 10), xlims=(1e-1, 100), scale=:log10)
xlabel!("I")
ylabel!("P(I)")
title!("variance = $(  round(var(all_intensities_over_mean),digits=3 ))")
