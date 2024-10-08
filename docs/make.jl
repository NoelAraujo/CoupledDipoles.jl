# Inside make.jl
push!(LOAD_PATH,"../src/")
using Documenter, CoupledDipoles
makedocs(;
    sitename = "CoupledDipoles.jl",
    authors = "Noel Araujo Moreira",
    pages=[
        "Home" => "index.md",
        "Atom" => "create_atoms/atoms.md",
        "Laser" => "lasers/lasers.md",
        "Problem" => "problems/problems.md",
        "Steady State" => "steady_state.md",
        "Scatttering" => "scattering/scattering.md",
        "Single Atom" => "dipole_example/single_atom_volume.md",
        "Intensity Statistics" => "variances_angles/rayleigh_variance.md",
        "Transmission" => "transmission/transmission.md",
        "Time Evolution" => "time_evolution/time_evolution.md",
        "Spectrum" => "spectrum/spectrum.md",
        "Spatial Profile" => "localization/spatial_modes.md",
        "Fitting" => "linear_regression/linear_regression_methods.md",
   ])


deploydocs(;
    repo="github.com/NoelAraujo/CoupledDipoles.jl.git",
)
