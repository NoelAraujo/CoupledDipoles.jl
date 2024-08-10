# Inside make.jl
push!(LOAD_PATH,"../src/")
using CoupledDipoles
using Documenter
makedocs(
         sitename = "CoupledDipoles.jl",
         modules  = [CoupledDipoles],
         authors = "Noel Araujo Moreira",
         pages=[
                "Home" => "index.md",
                "Atom" => "create_atoms/atoms.md",
                "Steady State" => "steady_state.md",
                "Scattering" => "scattering.md",
                "Single Atom" => "dipole_example/single_atom_volume.md",
                "Intensity Statistics" => "variances_angles/variaces_different_angles_and_shapes.md",
                "Fitting" => "linear_regression/linear_regression_methods.md",
               ])
deploydocs(;
    repo="github.com/NoelAraujo/CoupledDipoles.jl.git",
)