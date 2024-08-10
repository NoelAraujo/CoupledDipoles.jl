# Inside make.jl
push!(LOAD_PATH,"../src/")
using Documenter, CoupledDipoles
makedocs(;
    sitename = "CoupledDipoles.jl",
    authors = "Noel Araujo Moreira",
    pages=[
        "Home" => "index.md",
        "Steady State" => "steady_state.md",
        "Scatttering" => "scattering.md",
        "Single Atom" => "dipole_example/single_atom_volume.md",
        "Intensity Statistics" => "variances_angles/rayleigh_variance.md",
        "Fitting" => "linear_regression/linear_regression_methods.md",
        "Atom" => "create_atoms/atoms.md",
   ])


deploydocs(;
    repo="github.com/NoelAraujo/CoupledDipoles.jl.git",
)
