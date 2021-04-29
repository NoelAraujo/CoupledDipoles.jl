module CoupledDipole

export functionHello
functionHello() = println("Huge Hello World")

using LinearAlgebra
using StaticArrays
using Distances
using Plots, LaTeXStrings

try
    using ArrayFire
catch
    println("no ArrayFire installed")
end
using DifferentialEquations
using ProgressMeter
using Random
using Statistics

using Crayons
using Printf
using LsqFit
using Optim: minimizer, optimize
using SpecialFunctions: besseli

export greet

const k₀ = 1
const Γ = 1
const probability_threshold = 0.999
const R1_threshold = 0.5

include("structs.jl")
export TwoD, ThreeD, Scalar
export Cube, Sphere
export SimulationScalar, SimulationMeanField
export PlaneWave_3D, Gaussian_3D
export ScalarProblem, MeanFieldProblem

include("atoms.jl")
export cube_inputs, sphere_inputs
export get_pairwise_matrix
export convert_matrix_to_StaticArray

include("kernels.jl")
export green_scalar!

include("lasers.jl")
export laser_over_atoms, estimate_waist

include("eigen_analysis.jl")
export get_interaction_matrix
export get_spectrum
export get_energy_shift_and_linewith
export get_IPRs, get_PRs
export get_all_ξ_and_R1
export get_spatial_profile_single_mode


include("exponential_fit.jl")

include("plot_modes.jl")
export plot_atoms_and_mode
export get_spatial_profile_data_to_save
export get_atoms_matrix
export get_X_axes, get_Y_axes, get_Z_axes

include("scattering.jl")
export get_scattered_intensity

include("sensors.jl")
export ring_on_space, ring_on_plane

include("dynamics.jl")
export get_steady_state
export time_evolution


end
