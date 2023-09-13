module CoupledDipoles

using LinearAlgebra
using Distances

using LinRegOutliers
using DataFrames
using OrdinaryDiffEq
using Logging: global_logger
using TerminalLoggers: TerminalLogger
global_logger(TerminalLogger())

using ThreadsX
using Printf
using LsqFit: curve_fit, coef
using Optim: minimizer, optimize, Options
using HCubature
## MKL improves performance for Intel, but degrades performance for AMD
## --> only for reference: Intel vendor is called 'GenuineIntel', and AMD is 'AuthenticAMD' on Linux
@static if Sys.iswindows()
    vendor_id = readchomp(`wmic cpu get name`)
    if occursin("Intel", vendor_id)
        using MKL
    end
elseif Sys.islinux()
    vendor_id = readchomp(pipeline(`lscpu`))
    if occursin("GenuineIntel", vendor_id)
        using MKL
    end
end
using NLsolve
using Distributions: MvNormal
using Bessels
using Random

using ProgressMeter
using Tullio

const k₀ = 1
const Γ = 1
const FARFIELD_FACTOR = 100
const LASER_FACTOR = -im/2

include("structs.jl")
include("atoms/atom_cube.jl")
include("atoms/atom_cylinder.jl")
include("atoms/atom_sphere.jl")

include("formulas.jl")
export ρ_of, b₀_of, single_sensor

include("IO.jl")
export Atom
export Square, Circle
export Cube, Sphere, Cylinder

export get_dimension
export Laser
export PlaneWave2D, PlaneWave3D
export Gaussian2D, Gaussian3D

export Physics, Linear
export LinearOptics, Scalar, Vectorial
export NonLinearOptics, MeanField, PairCorrelation

include("linear/linear_problems.jl")
include("nonlinear/nonlinear_problems.jl")
export turn_laser_off!, turn_off!

include("atoms/atoms.jl")
export cube_inputs, sphere_inputs, cylinder_inputs
export radius_of_exclusion
export get_pairwise_matrix

include("linear/kernels.jl")
export interaction_matrix
export green_scalar!

include("lasers.jl")
export laser_field, laser_intensity
export apply_laser_over_sensors, apply_laser_over_oneSensor
export apply_laser_over_oneAtom, laser_angles

include("scattered_fields.jl")
include("scattered_fields_implementations.jl")
include("laser_fields.jl")
include("scattering.jl")
include("linear/scalar_pump_scattering.jl")
include("linear/vectorial_pump_scattering.jl")
include("nonlinear/meanfield_pump_scattering.jl")
include("nonlinear/paircorrelation_pump_scattering.jl")
export scattered_electric_field, scattered_intensity, laser_and_scattered_intensity, laser_and_scattered_electric_field


include("integrate_over_angle.jl")
export get_intensity_over_an_angle

include("eigenanalysis/eigen_analysis.jl")
export get_spectrum
export get_IPRs, get_PRs
export get_localization_length
export get_spatial_profile_single_mode
export classify_modes

include("eigenanalysis/exponential_fit.jl")
export select_points


include("sensors.jl")
export get_sensors_ring
export get_sensors_sphereSurface

include("initial_conditions.jl")
export default_initial_condition

include("linear/time_evolution_linear.jl")
include("linear/steady_state_linear.jl")
include("nonlinear/time_evolution_nonlinear.jl")
include("nonlinear/steady_state_nonlinear.jl")
export steady_state
export time_evolution, get_evolution_function

include("transmission.jl")
export transmission, how_far_is_FarField
export _create_sphere_sensor, _create_plane_sensor

include("powers.jl")
export scattered_power

include("applications/cbs.jl")
export CBS_scalar


import SnoopPrecompile
SnoopPrecompile.@precompile_all_calls begin

    N = 40
    kL = 32.4
    w₀, s, Δ = 4π, 1e-5, 0.3
    tspan = (0, 15.0)

    atoms = Atom(Sphere(gaussian=true), N, kL; r_min=0.0)
    laser = Laser(Gaussian3D(w₀), s, Δ)
    sensors = get_sensors_ring(; num_pts = 10, kR = 300, θ = 5π / 12)

    simulation = LinearOptics(Scalar(), atoms, laser)
    βₙ = steady_state(simulation)
    scattered_intensity(simulation, βₙ, sensors; regime=:near_field)
    scattered_intensity(simulation, βₙ, sensors; regime=:far_field)
    G = interaction_matrix(simulation)
    ωₙ, Γₙ = get_spectrum(simulation)
    u₀ = default_initial_condition(simulation)
    βₜ = time_evolution(simulation, u₀, tspan) # using ODE solver
    βₜ = time_evolution(simulation, u₀, tspan; ode_solver=false) # using FORMAL solution
    P_total = scattered_power(simulation, βₜ.u[end])

    laser = Laser(Gaussian3D(w₀), s, Δ; polarization=[1,0,0])
    simulation = LinearOptics(Vectorial(), atoms, laser)
    G = interaction_matrix(simulation)
    βₙ = steady_state(simulation)
    scattered_intensity(simulation, βₙ, sensors; regime=:near_field)
    scattered_intensity(simulation, βₙ, sensors; regime=:far_field)

    simulation = NonLinearOptics(MeanField(), atoms, laser)
    u₀ = default_initial_condition(simulation)
    βₜ = time_evolution(simulation, u₀, tspan)
    βₛ = steady_state(simulation)
    scattered_intensity(simulation, βₛ, sensors; regime=:near_field)
    scattered_intensity(simulation, βₛ, sensors; regime=:far_field)
    P_total = scattered_power(simulation, βₛ)
end

end
