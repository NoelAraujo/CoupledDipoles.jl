module CoupledDipoles

using LinearAlgebra
using Distances
# using Plots, LaTeXStrings

using DifferentialEquations
using Logging: global_logger
using TerminalLoggers: TerminalLogger
global_logger(TerminalLogger())
# using ProgressMeter
# using Random
# using Statistics

using Printf
using LsqFit
using Optim: minimizer, optimize, Options
using SpecialFunctions
using Folds
using ThreadPools
using LazyArrays
using Memoize
using SharedArrays

const k₀ = 1
const Γ = 1
const probability_threshold = 0.9999
const R1_threshold = 0.5

include("structs.jl")
include("atom_cube.jl")
include("atom_sphere.jl")
include("formulas.jl")
include("IO.jl")
export Atom
export Square, Circle
export Cube, Sphere

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
export cube_inputs, sphere_inputs
export get_rₘᵢₙ
export get_pairwise_matrix

include("kernels.jl")
export get_interaction_matrix
export green_scalar!

# include("lasers.jl")
# export laser_over_atoms, estimate_waist
include("pump_gaussians.jl")
include("pump_planewaves.jl")
export apply_laser_over_atoms, apply_laser_over_oneAtom

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
export get_intensities_over_sensors
export get_intensity_over_an_angle
export get_intensity_over_an_angle_legacy

include("sensors.jl")
export get_sensors_ring

include("dynamics.jl")
export get_steady_state
export time_evolution, default_evolution_initial_condition, get_evolution_function

#=
    One needs attention to use @memoize not forget to update this function
=#
# function clear_cache_for_long_simulations(problem)
#     empty!(memoize_cache(get_xyz_distances))
#     nothing
# end
# export clear_cache_for_long_simulations

end
