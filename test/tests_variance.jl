@time begin
    using CairoMakie, LinearAlgebra
    using Statistics: mean, var
    using StatsBase: fit, normalize, Histogram
end
@time begin
    using CoupledDipoles
    using ProgressMeter, Random
end

### ------------ ATOMS SPECS ---------------------
L = 32.4
N = [284, 700]


### ------------ LASER SPECS ---------------------
Δ = 1.0
s = 1e-5
w₀ = L / 4

### ------------ SIMULATION SPECS ---------------------
sensors = get_sensors_ring(; num_pts = 72, kR = 300, θ = 5π / 12)
maxRep = 25



# let 
#     Random.seed!(1134)
#     N = 10
#     atoms = Atom(Cube(), N, L)
#     laser = Laser(Gaussian3D(w₀), s, Δ)
#     simulation = LinearOptics(Vectorial(), atoms, laser)

#     βₙ = steady_state(simulation)
# end


### -------- PRODUCE INTENSITIES -----------------
all_intensities = map(N) do N
    many_intensities = @showprogress map(1:maxRep) do rep
        Random.seed!(1134 + rep)

        atoms = Atom(Cube(), N, L)
        laser = Laser(Gaussian3D(w₀), s, Δ)
        simulation = LinearOptics(Vectorial(), atoms, laser)

        βₙ = steady_state(simulation)
        intensities = scattered_intensity(simulation, βₙ, sensors; regime = :near_field, use_sequencial=false)

        intensities
    end

    many_intensities = reduce(vcat, many_intensities)
    all_intensities_over_mean = many_intensities ./ mean(many_intensities)

    all_intensities_over_mean
end;


### ------------ CREATE HISTOGRAM ---------------------
bins = 10.0 .^ range(log10(1e-6), log10(75); length = 30)

xy_data = map(eachindex(N)) do n
    h = fit(Histogram, all_intensities[n], bins)

    h_norm = normalize(h; mode = :pdf)
    bins_edges = collect(h_norm.edges[1])
    bins_centers = [sqrt(bins_edges[i] * bins_edges[i+1]) for i = 1:(length(bins_edges)-1)]
    variance = var(all_intensities[n])

    # x_data_histogram, y_data_histogram, variance
    (bins_centers, h_norm.weights, variance)
end


### ------------ FINAL PLOT ---------------------
begin
    fig = Figure(size = (800, 450))
    ax = Axis(
        fig[1, 1],
        xlabel = "Intensity",
        ylabel = "Probability Distribution",
        title = "",
        xlabelsize = 25,
        ylabelsize = 25,
        xticklabelsize = 20,
        yticklabelsize = 20,
        xscale = log10,
        yscale = log10,
    )

    ## theoretical curve
    x_ray = range(0.01, 50; step = 0.15)
    y_ray = exp.(-x_ray)
    lines!(ax, x_ray, y_ray, linestyle = :dash, label = "Rayleigh", color = :black, lw = 4)

    for n = 1:2
        x = xy_data[n][1]
        y = xy_data[n][2]
        v = xy_data[n][3] # variance
        notNull = findall(y .> 0)
        scatter!(
            ax,
            x[notNull],
            y[notNull];
            label = "N=$(N[n]), Variance = $( round(v,digits=3 ))",
            markershape = :circle,
            markersize = 20,
        )
    end
    ylims!(1e-6, 10)
    xlims!(1e-2, 100)
    axislegend(position = :rt, labelsize = 20)
    fig
end