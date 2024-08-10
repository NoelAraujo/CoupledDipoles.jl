# CoupledDipoles.jl

## **Installation**
This is package is not registred in julia, you have to use its github URL

Due to its many dependencies, the installation process may take around 10 minutes and requires at least 16Gb of RAM on your device
```julia
import Pkg
# Pkg.add("MKL") # (if installation goes wrong, run this line, and try to install again)
Pkg.add(url="https://github.com/NoelAraujo/CoupledDipoles.jl")
```

To verify that everything is working correctly, run:

```julia
using CoupledDipoles
# Pkg.test("CoupledDipoles")
```

## First Example

Let's create a problem with `N=1000` atoms inside a `Cube` of size `kR=40`, pumped by a laser with saturation `s=1e-6` and on resonance `Δ=0`. We'll use the `Atom` and `Laser` constructors to define the `problem`, and then create a `Scalar Model`.

From the `problem` you get, for example, the [`steady_state`](@ref ss_page)

```julia
using CoupledDipoles
N, kR = 1000, 40
atoms = Atom(Cube(), N, kR)

s, Δ = 1e-6, 0.0
laser = Laser(PlaneWave3D(), s, Δ)
problem = LinearOptics(Scalar(), atoms, laser)
βₛₛ = steady_state(problem)
```




## Manual Outline
```@contents
```

