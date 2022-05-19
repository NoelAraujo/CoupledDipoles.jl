using CoupledDipoles
using Random
using Statistics: mean, var
using ProgressMeter, Plots
Random.seed!(1134)

function createAtoms(atomType, atomInfo, laser; maxRep=40)
    simulations, atomicStates = [], []
    p = Progress(maxRep; showspeed=true)
    for rep in 1:maxRep
        simulation, βₛ = getState(atomType, atomInfo, laser)

        push!(simulations, simulation)
        push!(atomicStates, βₛ)

        ProgressMeter.next!(p)
    end
    return simulations, atomicStates
end
function getState(atomType, atomInfo, laser)
    atoms = Atom(atomType, atomInfo...)
    simulation = LinearOptics(Scalar(), atoms, laser)

    βₛ = steady_state(simulation)
    return simulation, βₛ
end

function computeVariances(simulations, atomicStates, sensorRadius, θ_range; maxRep=40)
    variances = zeros(length(θ_range))

    p = Progress(length(θ_range) * maxRep; showspeed=true)
    for (θ_idx, θ) in enumerate(θ_range)
        sensors = get_sensors_ring(; num_pts=3 * 360, kR=sensorRadius, θ=θ)
        variances[θ_idx] = getOneVariance(simulations, atomicStates, sensors, p)
    end
    return variances
end
function getOneVariance(simulations, atomicStates, sensors, p)
    intensities = Float64[]
    for (simulation, βₛ) in zip(simulations, atomicStates)
        oneRepetition = getIntensities(simulation, βₛ, sensors)
        append!(intensities, oneRepetition)

        ProgressMeter.next!(p)
    end
    variance = computeVariance(intensities)
    return variance
end
@inline function getIntensities(simulation, βₙ, sensors)
    scatt_func = scattering_fuction(:farField, :ThreeD)
    intensities = scattering_intensity(simulation, βₙ, sensors, scatt_func)
    return intensities
end
@inline function computeVariance(intensities)
    intensities_over_mean = intensities ./ mean(intensities)
    return var(intensities_over_mean)
end

const k₀ = 1;
const λ = 2π / k₀;

# atoms setup
ρλ³ = 45
ρ_over_k₀³ = ρλ³ / (2π)^3
N = 7500
factorR = 5
kR = factorR * λ

# laser setup
factorλ = 2
w₀, s, Δ = factorλ * λ, 1e-6, 0.77
laserinfo = w₀, s, Δ
laser = Laser(Gaussian3D(w₀), s, Δ)

# simulation setup
maxRep = 50
sensorRadius = 250 * kR
θ_range = range(0, π; length=100)

cloudGeometries = [CoupledDipoles.Cube(), CoupledDipoles.Sphere(), CoupledDipoles.Cylinder()]
cloudParams = [cube_inputs(N, ρ_over_k₀³), sphere_inputs(N, ρ_over_k₀³), cylinder_inputs(N, ρ_over_k₀³; R=kR)]
lstyles = [:dot, :dash, :solid]

# making simulation
toSave = zeros(length(cloudGeometries), length(θ_range))
for (atomType, atomInfo, idx) in zip(cloudGeometries, cloudParams, 1:length(cloudGeometries))
    println("** Cloud: $(typeof(atomType))")

    simulations, atomicStates = createAtoms(atomType, atomInfo, laser; maxRep=maxRep)
    singleVarianceCurve = computeVariances(simulations, atomicStates, sensorRadius, θ_range; maxRep=maxRep)

    toSave[idx, :] .= singleVarianceCurve
end

# making figure
plot(; xlabel="θ", ylabel="Variance", ylims=(1e-4, 50))
for (atomType,  idx) in zip(cloudGeometries, 1:length(cloudGeometries))
    singleVarianceCurve = toSave[idx, :]
    plot!(θ_range, singleVarianceCurve; label="$(typeof(atomType))", lw=4, linestyle=lstyles[idx], yscale=:log10)
end

vline!([π / 4]; linestyle=:dash, c=:black, label="")
vline!([π / 2]; linestyle=:dash, c=:red, label="")
hline!([1]; c=:black, label="")
title!("N=$(N), kR=$(factorR), w₀=$(factorλ)λ₀, Δ=$( round(Δ, digits=2) ), $(maxRep) reps")
plot!(; xticks=([0, π / 4, π / 2, π], ["0", "π/4", "π/2", "π"]))
display(plot!(; yticks=([1e-3, 1e-1, 1, 50], ["1e-3", "1e-1", "1", "50"])))
