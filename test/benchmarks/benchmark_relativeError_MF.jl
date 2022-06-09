using Distributed
ip_m2 = "pc2@192.168.15.8"
nprocess_m2 = 4
machine2 = [(ip_m2, nprocess_m2)]
addprocs(machine2; topology=:master_worker, exeflags=`--threads 3`, enable_threaded_blas=true, exename="/home/pc2/julia/bin/julia", dir="/home/pc2")

ip_m3 = "pc3@192.168.15.7"
nprocess_m3 = 4
machine3 = [(ip_m3, nprocess_m3)]
addprocs(machine3; topology=:master_worker, exeflags=`--threads 2`, enable_threaded_blas=true, exename="/home/pc3/julia_program/julia-1.7.2/bin/julia", dir="/home/pc3")

@everywhere using CoupledDipoles
[@spawnat w include(homedir() * "/usar_GPU.jl") for w in workers()]

using Random
using ProgressMeter
const k₀ = 1;
const λ = 2π / k₀;

function makeClouds(atomType, atomInfo, maxRep)
    clouds = Atom[]
    for rep in 1:maxRep
        push!(clouds, Atom(atomType, atomInfo...))
    end
    return clouds
end
function getMap(clouds, laserInfo, Δ_range, time_range)
    maxRep = length(clouds)

    p = Progress(length(Δ_range) * maxRep; showspeed=true)

    finalMap = zeros(length(Δ_range), length(time_range))
    for rep in 1:maxRep
        println("rep: $(rep)/$(maxRep)")
        finalMap += errorsOverTime(clouds[rep], laserInfo, Δ_range, time_range)
    end
    return finalMap / maxRep
end

function errorsOverTime(atoms, laserInfo, Δ_range, time_range)
    w₀, s = laserInfo

    oneCloudManyDetunings = zeros(length(Δ_range), length(time_range))
    xx = @showprogress pmap(eachindex(Δ_range)) do idx
        Δ = Δ_range[idx]
        laser = Laser(Gaussian3D(w₀), s, Δ)
        oneCloudManyDetunings[idx, :] .= relatriveErrors_singleΔ(atoms, laser, time_range)
    end

    [oneCloudManyDetunings[idx, :] .= x for (idx, x) in enumerate(xx)]
    return oneCloudManyDetunings
end
@everywhere function relatriveErrors_singleΔ(atoms, laser, time_range)
    pmf = NonLinearOptics(MeanField(), atoms, laser)
    u₀ = default_initial_condition(pmf)
    tspan = (time_range[1], time_range[end])

    sol = time_evolution(pmf, u₀, tspan; saveat=time_range, progress=false)

    errors_abs2 = [sum(abs2, sol.u[t] - sol.u[end]) / sum(abs2, sol.u[end]) for t in 1:length(sol.u)]

    return errors_abs2
end

Random.seed!(1134)
# atom setup
N = 1000
ρ = 0.181414
factorR = 5
kR = factorR * λ

# laser setup
factorW₀ = 2
w₀ = factorW₀ * λ
s = 1e-4

# simulation setup
laserInfo = w₀, s

atomType = CoupledDipoles.Cylinder()
atomInfo = cylinder_inputs(N, ρ; R=kR)

maxRep = 20
maxΔs = 40
nTimeSteps = 100

time_range = range(0.0, 10_000; length=nTimeSteps)
Δ_range = range(-10, 10; length=maxΔs)

# make simulation
clouds = makeClouds(atomType, atomInfo, maxRep);
mapErrors = getMap(clouds, laserInfo, Δ_range, time_range);

using Plots;
plotly();
plot(; size=(800, 600), guidefont=17, tickfont=15)
heatmap!(time_range, Δ_range, log10.(mapErrors); cbtitle="log10( relative error )", clims=(-12, -7)) # 
xlabel!("t/Γ")
ylabel!("Δ/Γ")
display(title!("$(typeof(atomType)): N=$(N), kR=$(factorR)λ,  w₀=$(factorW₀)λ, s=$(s), $(maxRep) reps"))
