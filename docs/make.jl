# Inside make.jl
push!(LOAD_PATH,"../src/")
using CoupledDipoles
using Documenter
makedocs(
         sitename = "CoupledDipoles.jl",
         modules  = [CoupledDipoles],
         pages=[
                "Home" => "index.md"
               ])
deploydocs(;
    repo="github.com/NoelAraujo/CoupledDipoles.jl",
)