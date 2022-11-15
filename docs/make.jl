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
                "Steady State" => "steady_state.md",
               ])
deploydocs(;
    repo="github.com/NoelAraujo/CoupledDipoles.jl.git",
)