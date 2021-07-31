# CoupledDipoles

[![Build Status](https://travis-ci.com/NoelAraujo/CoupledDipole.jl.svg?branch=master)](https://travis-ci.com/NoelAraujo/CoupledDipole.jl)
[![Build Status](https://ci.appveyor.com/api/projects/status/github/NoelAraujo/CoupledDipole.jl?svg=true)](https://ci.appveyor.com/project/NoelAraujo/CoupledDipole-jl)
[![Coverage](https://codecov.io/gh/NoelAraujo/CoupledDipole.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/NoelAraujo/CoupledDipole.jl)
[![Coverage](https://coveralls.io/repos/github/NoelAraujo/CoupledDipole.jl/badge.svg?branch=master)](https://coveralls.io/github/NoelAraujo/CoupledDipole.jl?branch=master)


### Very experimental: Be aware that syntax may change
Based on publication *N.A. Moreira*. **Localization versus subradiance in three-dimensional scattering of light**. EPR. 2019
```julia
using CoupledDipoles
using Random; Random.seed!(1135)

#= Spherical Cloud and Plane Wave laser at z-direction =#
N, ρ = 2_000, 1.0
atoms = Atom(Sphere(), sphere_inputs(N, ρ)...)

s, Δ = 1e-5, 0.0
laser = Laser(PlaneWave3D([0,0,1]), s, Δ)
simulation = LinearOptics(Scalar(), atoms, laser)

ωₙ, Γₙ = get_spectrum(simulation)
modes = classify_modes(simulation)
```

After some customization (available at */images/Fig1.jl*), one may get the Spatial Profiles for each mode.

![GitHub Logo](/images/Fig1.png)