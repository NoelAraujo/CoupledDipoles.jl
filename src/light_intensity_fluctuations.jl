using TestingNewCodingStyle
using Statistics: mean, var
using StatsBase: fit, normalize, Histogram
using ProgressMeter
using Plots, FileIO
# using Revise

using Distributed
addprocs(3)
# addprocs(4, exeflags=`--threads 3`, enable_threaded_blas=true)

# machine3 = [("pc3@192.168.15.7",5)]
# addprocs(machine3;
#     multiplex=true,
#     exeflags=`--threads 2`,
#     enable_threaded_blas=true, dir="/home/pc3")

# @everywhere push!(LOAD_PATH, pwd())
# @everywhere using TestingNewCodes

fullPATH = homedir() * "/Documentos/TestingNewCodingStyle"
for p in workers()
    fetch(@spawnat p cd(fullPATH))
    fetch(@spawnat p push!(LOAD_PATH, pwd()))
end
@everywhere using TestingNewCodingStyle

@everywhere function get_intensity_sample(atoms, s, Δ, sensors)
    # atoms = Cube(N, kL)
    laser = Gaussian_3D(estimate_waist(atoms), s, Δ)
    simulation = ScalarProblem(atoms, laser)

    βₙ = steady_state(simulation)
    intensities = get_scattered_intensity(simulation, βₙ, sensors)
    GC.gc()
    return intensities
end

## data that I will use
ρλ³_range = range(5, 45; length=5)
maxRepAvailable = 10
kL = 32.4;

# each density needs different repetitions to convergence
low_density_rep = maxRepAvailable
high_density_rep = maxRepAvailable / 2
rep_range = round.(Int, range(low_density_rep, high_density_rep; length=length(ρλ³_range)))

s = 1e-6
nSensors = 360
sensors = ring_on_space(; num_pts=nSensors, kR=1.5kL, θ=5π / 12)
Δ_range = range(0, 2; length=3)

variances = zeros(length(ρλ³_range), length(Δ_range))
for (iρ, ρλ³) in enumerate(ρλ³_range)
    ρ_over_k₀³ = ρλ³ / (2π)^3
    N = floor(Int, ρ_over_k₀³ * kL^3)

    atoms_saved = load(fullPATH * "/saved/Cubes_kL=$(kL), rho=$(ρλ³), maxRepAvailable=$(maxRepAvailable).jld2", "atoms_saved")
    maxRep = rep_range[iρ]

    for (iΔ, Δ) in enumerate(Δ_range)
        # println(iρ, " -- ", iΔ)
        println("density: $(iρ)/$(length(ρλ³_range)) -- detuning: $(iΔ)/$(length(Δ_range))")
        # @showprogress
        many_intensities = pmap(1:maxRep) do rep
            intensities = get_intensity_sample(atoms_saved[rep], s, Δ, sensors)
        end

        all_intensities = Float64[]
        for i in many_intensities
            for j in 1:length(i)
                push!(all_intensities, i[j])
            end
        end
        all_intensities_over_mean = all_intensities ./ mean(all_intensities)

        variances[iρ, iΔ] = var(all_intensities_over_mean)
    end
    # heatmap(Δ_range, ρλ³_range, variances,
    #     c=:jet, xlabel="Δ/Γ", ylabel="ρλ³",
    #     title="s=$(s), kL=$(kL)"
    # ) |> display
end

# save("heatmap_variances.jld2", Dict("variances" => variances))

using FileIO
variances = load("src/heatmap_variances.jld2", "variances")

ρλ³_range = range(5, 45; length=5)
Δ_range = range(-1, 2; length=5)

using Plots

#Δ_range, ρλ³_range,
heatmap(variances; c=:jet, xlabel="Δ/Γ", ylabel="ρλ³", title="")
# savefig("florent_middle_resolution.png")
