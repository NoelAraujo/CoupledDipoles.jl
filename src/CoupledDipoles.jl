module CoupledDipoles

using LinearAlgebra
using Distances

using DifferentialEquations
using Logging: global_logger
using TerminalLoggers: TerminalLogger
global_logger(TerminalLogger())
# using Random
# using Statistics
using ThreadsX

using Printf
using LsqFit
using Optim: minimizer, optimize, Options
using SpecialFunctions
using Folds
using ThreadPools
using LazyArrays
# using Memoize
using SharedArrays
using HCubature
using MKL

const k₀ = 1
const Γ = 1
const probability_threshold = 0.9999
const R1_threshold = 0.5
const farField_factor = 100

include("structs.jl")
include("atom_cube.jl")
include("atom_sphere.jl")
include("atom_cylinder.jl")

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

include("linear_problems.jl")
include("nonlinear_problems.jl")


include("atoms.jl")
export cube_inputs, sphere_inputs, cylinder_inputs
export get_rₘᵢₙ
export get_pairwise_matrix

include("kernels.jl")
export interaction_matrix
export green_scalar!

include("lasers.jl")
export laser_field, apply_laser_over_oneAtom
export apply_laser_over_sensors, apply_laser_over_oneSensor

include("eigen_analysis.jl")
export get_spectrum
export get_IPRs, get_PRs
export get_localization_length
export get_spatial_profile_single_mode
export classify_modes

include("exponential_fit.jl")
export select_points

include("plot_modes.jl")
export get_spatial_profile_single_mode_and_fit
export get_decay_fit


include("scattering.jl")
export scattering_intensity, scattering_fuction 
export get_intensity_over_an_angle
export get_intensity_over_an_angle_legacy

include("sensors.jl")
export get_sensors_ring
export get_sensors_sphereSurface

include("dynamics.jl")
export steady_state
export time_evolution, default_evolution_initial_condition, get_evolution_function

include("transmission.jl")
export transmission, how_far_is_FarField
export _create_sphere_sensor, _create_plane_sensor

end
