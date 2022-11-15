# Inside make.jl
push!(LOAD_PATH,"../src/")
import Pkg; Pkg.add(url="https://github.com/JuliaMath/Bessels.jl")
using CoupledDipoles
using Documenter
makedocs(
         sitename = "CoupledDipoles.jl",
         modules  = [CoupledDipoles],
         pages=[
                "Home" => "index.md"
               ])
deploydocs(;
    repo="github.com/NoelAraujo/CoupledDipoles.jl.git",
)