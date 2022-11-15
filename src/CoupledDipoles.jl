module CoupledDipoles

using LinearAlgebra
using Distances

# using DifferentialEquations
using OrdinaryDiffEq
using Logging: global_logger
using TerminalLoggers: TerminalLogger
global_logger(TerminalLogger())

using ThreadsX
using Printf
using LsqFit: curve_fit, coef
using Optim: minimizer, optimize, Options
using HCubature
# using MKL ## in some servers I can download this package
using ParallelStencil
@init_parallel_stencil(Threads, Float64, 2);

using Distributions: MvNormal
using Bessels

const k₀ = 1
const Γ = 1
const PROBABILITY_THRESHOLD = 0.9999
const R1_THRESHOLD = 0.5
const FARFIELD_FACTOR = 100

include("structs.jl")
include("atoms/atom_cube.jl")
include("atoms/atom_cylinder.jl")
include("atoms/atom_sphere.jl")

include("formulas.jl")
export ρ_of, b₀_of

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
export NonLinearOptics, MeanField, BBGKY

include("linear/linear_problems.jl")
include("nonlinear/nonlinear_problems.jl")
export turn_laser_off!, turn_off!

include("atoms/atoms.jl")
export cube_inputs, sphere_inputs, cylinder_inputs
export get_rₘᵢₙ
export get_pairwise_matrix

include("linear/kernels.jl")
export interaction_matrix
export green_scalar!

include("lasers.jl")
include("linear/scalar_pump_scattering.jl")
include("linear/vectorial_pump_scattering.jl")
export scattered_electric_field, scattered_intensity
export laser_field, apply_laser_over_oneAtom
export apply_laser_over_sensors, apply_laser_over_oneSensor

include("eigenanalysis/eigen_analysis.jl")
export get_spectrum
export get_IPRs, get_PRs
export get_localization_length
export get_spatial_profile_single_mode
export classify_modes

include("eigenanalysis/exponential_fit.jl")
export select_points

include("plot_modes.jl")
export get_spatial_profile_single_mode_and_fit
export get_decay_fit

include("scattering.jl")
include("linear/linear_scattering.jl")
include("nonlinear/nonlinear_scattering.jl")
export scattering_intensity, scattering_fuction
export get_intensity_over_an_angle
export get_intensity_over_an_angle_legacy

include("sensors.jl")
export get_sensors_ring
export get_sensors_sphereSurface

include("linear/linear_dynamics.jl")
include("nonlinear/nonlinear_dynamics.jl")
export steady_state
export time_evolution, default_initial_condition, get_evolution_function

include("transmission.jl")
export transmission, how_far_is_FarField
export _create_sphere_sensor, _create_plane_sensor

import SnoopPrecompile
SnoopPrecompile.@precompile_all_calls begin

    N = 40
    kL = 32.4
    w₀, s, Δ = 4π, 1e-5, 0.3
    tspan = (0, 15.0)

    atoms = Atom(Sphere(gaussian=true), N, kL; r_min=0.0)
    laser = Laser(Gaussian3D(w₀), s, Δ)

    simulation = LinearOptics(Scalar(), atoms, laser)
    βₙ = steady_state(simulation)
    G = interaction_matrix(simulation)
    ωₙ, Γₙ = get_spectrum(simulation)
    u₀ = default_initial_condition(simulation)
    βₜ = time_evolution(simulation, u₀, tspan)


    simulation = LinearOptics(Vectorial(), atoms, laser)
    G = interaction_matrix(simulation)
    βₙ = steady_state(simulation)

    simulation = NonLinearOptics(MeanField(), atoms, laser)
    u₀ = default_initial_condition(simulation)
    βₜ = time_evolution(simulation, u₀, tspan)
end

end
