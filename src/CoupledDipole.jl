module CoupledDipole

using LinearAlgebra
using StaticArrays
using Distances
# using Plots, LaTeXStrings

# using DifferentialEquations
# using ProgressMeter
# using Random
# using Statistics

# using Crayons
using Printf
# using LsqFit
# using Optim: minimizer, optimize
# using SpecialFunctions: besseli
using Folds
using ThreadPools
using LazyArrays
using Memoize

const k₀ = 1
const Γ = 1
const probability_threshold = 0.9999
const R1_threshold = 0.5

include("structs.jl")
include("shape_cube.jl")
include("formulas.jl")
include("IO.jl")
export Shape 
export Square, Circle
export Cube, Sphere

export get_dimension
export Laser
export PlaneWave2D, PlaneWave3D
export Gaussian2D, Gaussian3D

export Physics, Linear
export LinearProblem, Scalar, Vectorial
export NonLinearProblem, MeanField, BBGKY

include("linear_problems.jl")
export myLinearProblem

include("atoms.jl")
# export cube_inputs, sphere_inputs
export get_rₘᵢₙ
export get_pairwise_matrix

include("kernels.jl")
export get_interaction_matrix
export green_scalar!

# include("lasers.jl")
# export laser_over_atoms, estimate_waist
include("pump_gaussians.jl")
export apply_laser_over_atoms, apply_laser_over_oneAtom

# include("eigen_analysis.jl")
# export get_interaction_matrix
# export get_spectrum
# export get_energy_shift_and_linewith
# export get_IPRs, get_PRs
# export get_all_ξ_and_R1
# export get_spatial_profile_single_mode
# export classify_modes

# include("exponential_fit.jl")

# include("plot_modes.jl")
# export plot_atoms_and_mode
# export get_spatial_profile_data_to_save
# export get_atoms_matrix
# export get_X_axes, get_Y_axes, get_Z_axes

include("scattering.jl")
export get_intensities_over_sensors
# export get_geometric_factor

include("sensors.jl")
export get_sensors_ring

include("dynamics.jl")
export get_steady_state
# export time_evolution


end
