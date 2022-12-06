using CoupledDipoles, Revise
using Random


using Statistics: mean, var
using StatsBase: fit, normalize, Histogram
using ProgressMeter, Plots, LinearAlgebra; plotly()

const λ = 2π
Δ_range=range(-1.5, 1.5; length=45)
ρ_range=range(5, 40; length=40)
kR=3λ
kh=4λ

### ------------ ATOMS SPECS ---------------------
ρλ³ = ρ_range[4] # ρ_range[24]
kL = 32.4;
ρ_over_k₀³ = ρλ³ / (2π)^3
N = floor(Int, ρ_over_k₀³ * kh * (π * kR^2))

### ------------ LASER SPECS ---------------------
Δ = Δ_range[26]
s = 1e-5
w₀ = 2λ

### -------- PRODUCE INTENSITIES -----------------
sensors = get_sensors_ring(; num_pts=360, kR=300λ, θ=deg2rad(76.1111111))
# scatt_func = scattering_fuction(:farField, :ThreeD)


maxRep = 64
many_intensities = Float64[]

p = Progress(maxRep; showspeed = true)
for rep = 1:maxRep
    Random.seed!(1134 + rep)
    atoms = Atom(Cylinder(), N, kR, kh)
    laser = Laser(Gaussian3D(w₀), s, Δ)
    # simulation = LinearOptics(Scalar(), atoms, laser)
    simulation = LinearOptics(Vectorial(), atoms, laser)

    βₙ = steady_state(simulation)
    intensities = scattered_intensity(simulation, βₙ, sensors; regime=:near_field)
    append!(many_intensities, intensities)

    ProgressMeter.next!(p)
end
all_intensities_over_mean = many_intensities ./ mean(many_intensities);

### ------------ CREATE HISTOGRAM ---------------------
bins = 10.0 .^ range(log10(1e-6), log10(100); length = 100)
h = fit(Histogram, all_intensities_over_mean, bins)

h_norm = normalize(h; mode = :pdf)
bins_edges = collect(h_norm.edges[1])
bins_centers = [sqrt(bins_edges[i] * bins_edges[i+1]) for i = 1:(length(bins_edges)-1)]


### ------------ FINAL PLOT ---------------------
## theoretical curve
x_ray = range(0.01, 100; step = 0.15)
y_ray = exp.(-x_ray)
plot(x_ray, y_ray; linestyle = :dash, label = "Rayleigh", lw = 4, size=(600,400))

notNull = findall( h_norm.weights .> 0)
scatter!(
    bins_centers[notNull],
    h_norm.weights[notNull];
    label = "N=$(N), kL=$(kL)",
    markershape = :circle,
    markersize = 5,
    scale=:log10
)
plot!(; guidefont = 15, tickfont = 15, legendfontsize = 10, size = (1000, 500))

plot!(; ylims = (1e-6, 10), xlims = (1e-1, 100))
xlabel!("I")
ylabel!("P(I)")
title!(
    "w₀=$(round(w₀,digits=2)),variance = $(  round(var(all_intensities_over_mean),digits=3 ))",
)
