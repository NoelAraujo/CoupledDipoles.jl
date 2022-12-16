using CairoMakie
using Distributed

addprocs(
    2;
    exeflags = `--project=$(Base.active_project()) --threads 2`,
    topology = :master_worker,
    enable_threaded_blas = true,
) # for GPU
addprocs(
    4;
    exeflags = `--project=$(Base.active_project()) --threads 2`,
    topology = :master_worker,
    enable_threaded_blas = true,
) # for CPU

@time @everywhere begin
    using CoupledDipoles, Random
    using ProgressMeter
    using Statistics: mean, var

    const λ = 2π
    function make_atoms(cloud_shape, N, ρ; w₀ = 2λ)
        if cloud_shape == :cube
            return Atom(CoupledDipoles.Cube(), cube_inputs(N, ρ)...)
        elseif cloud_shape == :sphere
            return Atom(CoupledDipoles.Sphere(), sphere_inputs(N, ρ)...)
        else
            cloud_shape == :cylinder
            return Atom(CoupledDipoles.Cylinder(), cylinder_inputs(N, ρ; R = w₀)...)
        end
    end
end


### ------------ LASER SPECS ---------------------
Δ = 1.0
s = 1e-5
w₀ = 8.0975

laser = Laser(Gaussian3D(w₀), s, Δ)


### ------------ ATOMS SPECS ---------------------
ρ = 44.3 / (2π)^3
N = 1500



### ------------ SIMPLE PARALLEL  ---------------------
function compute_variances_parallel_v1(cloud_shape, N, ρ, laser, θ_range, maxReps; w₀ = 2λ)

    # to get a progress message, we will monitor each angle
    variances = @showprogress map(θ_range) do θ
        sensors = get_sensors_ring(; num_pts = 180, kR = 1000, θ = θ)

        # each configuration is computed in parallel at once
        intensities = @distributed (vcat) for rep = 1:maxReps
            Random.seed!(1134 + rep)
            atoms = make_atoms(cloud_shape, N, ρ; w₀ = w₀)

            simulation = LinearOptics(Scalar(), atoms, laser)
            βₙ = steady_state(simulation)
            _ii = scattered_intensity(simulation, βₙ, sensors; regime = :far_field)
            _ii
        end

        var(intensities ./ mean(intensities))
    end
    return variances
end


### ------------ TWEAKED PARALLEL  ---------------------
@everywhere function get_state_state(input_channel, output_channel)
    while true
        cloud_shape, N, ρ, laser, rep, w₀, sensors = take!(input_channel)

        Random.seed!(1134 + rep)
        atoms = make_atoms(cloud_shape, N, ρ; w₀ = w₀)
        simulation = LinearOptics(Scalar(), atoms, laser)
        βₙ = steady_state(simulation)

        put!(output_channel, (simulation, βₙ, sensors))
    end
end
@everywhere function get_intensities(input_channel, output_channel)
    while true
        simulation, βₙ, sensors = take!(input_channel)

        _ii = scattered_intensity(simulation, βₙ, sensors; regime = :far_field)

        simulation = βₙ = 1 # clear memory
        put!(output_channel, _ii)
    end
end

## channels are the communication link from where one transfer data between process
channelSize = 10nworkers() # 10 times number of workers to avoid botlenecks (is an arbitrary number)
ss_input = RemoteChannel(() -> Channel(channelSize))
ii_input = RemoteChannel(() -> Channel(channelSize))
ii_output = RemoteChannel(() -> Channel(channelSize))

## process 2 and 3 will use GPU, therefore, they will compute steady state
# [@spawnat w ENV["COUPLED_DIPOLES_USE_GPU"] = false for w in workers()]
[@spawnat w use_gpu(true) for w in [2, 3]]

for i in [2, 3]
    @spawnat i get_state_state(ss_input, ii_input)
end

## in principles, all workers can compute the intensity
## but for this example, i will limit only process 4 and 5 to do so
for i in [4, 5, 6, 7] #workers()
    @spawnat i get_intensities(ii_input, ii_output)
end

function compute_variances_parallel_v2(cloud_shape, N, ρ, laser, θ_range, maxReps, ss_input, ii_output; w₀ = 2λ)
    p = Progress(maxReps * length(θ_range))

    ## put all data into the channels
    variances = map(θ_range) do θ
        sensors = get_sensors_ring(; num_pts = 180, kR = 1000, θ = θ)

        ## out configurations in the queue
        for rep = 1:maxReps
            simulations_inputs = cloud_shape, N, ρ, laser, rep, w₀, sensors
            @async put!(ss_input, simulations_inputs)
        end

        ## take the results
        intensities = mapreduce(vcat, 1:maxReps) do rep
            singleResult = take!(ii_output)
            ProgressMeter.next!(p)
            singleResult
        end
        var(intensities ./ mean(intensities))
    end
    return variances
end




### ------------ MAKE SIMULATION  ---------------------
θ_range = range(0, π, length = 15)
maxReps = 20


variances_cylinder_v1 = compute_variances_parallel_v1(:cylinder, N, ρ, laser, θ_range, maxReps; w₀ = w₀);
variances_cylinder_v2 =
    compute_variances_parallel_v2(:cylinder, N, ρ, laser, θ_range, maxReps, ss_input, ii_output; w₀ = w₀);

all(variances_cylinder_v1 .≈ variances_cylinder_v2)



### ------------ PLOT  ---------------------
fig = Figure()
ax = Axis(fig[1, 1], ylabel = "variance", xlabel = "θ",
    yscale=log10, yminorticksvisible = false,
    xticks = MultiplesTicks(5, pi, "π"),
)

lines!(θ_range, variances_cylinder_v1, label="v1", linewidth=4)
lines!(θ_range, variances_cylinder_v2, label="v2", linewidth=4, linestyle=:dash)
hlines!([1], color=:blue, linestyle=:dot, linewidth=4)

ylims!(1e-5, 10^3)
xlims!(0, 1.01π)
axislegend(position = :lt, labelsize=20)
current_figure()