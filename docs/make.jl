# Inside make.jl
push!(LOAD_PATH,"../src/")
import Pkg; Pkg.add(url="https://github.com/JuliaMath/Bessels.jl")
using CoupledDipoles
using Documenter
makedocs(
         sitename = "CoupledDipoles.jl",
         modules  = [CoupledDipoles],
         authors = "Noel Araujo Moreira",
         pages=[
                "Home" => "index.md",
                "Steady State" => "src/steady_state.md",
               ])
deploydocs(;
    repo="github.com/NoelAraujo/CoupledDipoles.jl.git",
)