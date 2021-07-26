function get_laser_over_atoms(laser::Laser{Gaussian3D}, atoms)
    @debug "applying Gaussian 3D over atoms"
    w₀, s, Δ = laser.pump.w₀, laser.s, laser.Δ

    E = zeros(Complex{eltype(atoms.r)}, atoms.N)
    E₀ = estimate_E₀(laser)
    v_atoms = view(Float64.(atoms.r), :, :)

    for n in 1:(atoms.N)  # Threads.@threads
        oneAtom = v_atoms[:,n]
        @inbounds E[n] = _core_Gaussian3D(oneAtom, E₀, w₀)
    end
    return E
end

function get_laser_over_oneAtom(laser::Laser{Gaussian3D}, position::AbstractArray)
    w₀, s, Δ = laser.pump.w₀, laser.s, laser.Δ

    E₀ = estimate_E₀(laser)
    E = _core_Gaussian3D(position, E₀, w₀)
    
    return E
end


estimate_E₀(laser) = √(laser.s * (1 + 4(laser.Δ / Γ)^2) / 2)

function _core_Gaussian3D(oneAtom, E₀, w₀)
    x, y, z = oneAtom[1], oneAtom[2], oneAtom[3]

    ## This formula is stable for z==0
    # Ref: Eq 3.11 from "CHAPTER 3. PROPAGATION AND FOCUSING OF OPTICAL FIELDS"
    denominator_factor = 1 .+ 2im .* z / (k₀ * w₀^2)
    Eᵢ = E₀ * exp(+im*k₀ * z)
    Eᵢ = Eᵢ .* exp.(-(x .^ 2 + y .^ 2) ./ (denominator_factor .* w₀^2))
    Eᵢ = Eᵢ ./ denominator_factor
    
    return Eᵢ
end

# @code_warntype get_sensor_intensities(laser, cloud, E, detectors_positions);
# @code_warntype laser_over_atoms(laser, cloud)


kL = 32.4;
s = 1e-6;
Γ = 1;
k₀ = 1;

# ρ = 5
ρ = 44;
Δ = 1.0;

ρ_over_λ³ = ρ / (2π)^3 # my "k₀=1", and not "λ=1"
N = floor(Int, ρ_over_λ³ * kL^3)

laser = Laser(Gaussian3D(kL/8), s, Δ)
sensors = get_sensors_ring(360; R=1.5kL, angle=5π/12)


using ProgressMeter
final_result = []
@showprogress for i=1:5
    new_cloud = Shape(Cube(), N, kL)
    
    G = Complex.(get_pairwise_matrix(new_cloud.r))
    LazyArrays.@~ G .= -(Γ / 2) * cis.(k₀ * G) ./ (1im * k₀ * G);
    G[diagind(G)] .= 1im * laser.Δ - Γ / 2
    
    Ω = (im/2)*get_laser_over_atoms(laser, new_cloud)
    
    βₛ = G\Ω
    
    one_data = get_intensities_over_sensors(laser, new_cloud, βₛ, sensors)
    for i ∈ one_data
        push!(final_result, copy(i))
    end
end
using Statistics: mean, var
using StatsBase: fit, normalize, Histogram
all_intensities_over_mean = final_result./mean(final_result)

bins = 10.0.^range(log10(1e-10),log10(100), length=100)
h = fit(Histogram, all_intensities_over_mean, bins)
h_norm = normalize(h, mode=:pdf)
bins_edges = collect(h_norm.edges[1])
bins_centers=[sqrt(bins_edges[i]*bins_edges[i+1]) for i =1:length(bins_edges)-1]

x_ray = range(0.01,10, step=0.15)
y_ray = exp.(-x_ray)

using Plots
plot(x_ray,y_ray, linestyle=:dash, label="Rayleigh", lw=4)
scatter!(bins_centers, h_norm.weights, label="Low Density", markershape=:circle,
markersize=10)
plot!(guidefont=17, tickfont=17, legendfont=15, size=(800, 600))
plot!(left_margin = 5Plots.mm,
    right_margin = 5Plots.mm,
    top_margin = 5Plots.mm,
    bottom_margin = 5Plots.mm,
    c = :blue,
    gridalpha = 0.8,
    minorgrid = true,
    minorgridalpha = 0.5,
    minorgridstyle = :dashdot,
    yscale = :log10,
    xscale = :log10,
    xlims = (10^-1,10^2),
    xticks = [10^-1,10^0,10^1,10^2],
    ylims = (10^-6,10^1),
    yticks = [10^1,10^0,10^-2,10^-4,10^-6]
)
plot!(ylims=(1e-6, 10), xlims=(1e-1, 100), scale=:log10)
xlabel!("I")
ylabel!("P(I)")
