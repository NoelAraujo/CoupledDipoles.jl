# CoupledDipoles

[![Build Status](https://travis-ci.com/NoelAraujo/CoupledDipole.jl.svg?branch=master)](https://travis-ci.com/NoelAraujo/CoupledDipole.jl)
[![Build Status](https://ci.appveyor.com/api/projects/status/github/NoelAraujo/CoupledDipole.jl?svg=true)](https://ci.appveyor.com/project/NoelAraujo/CoupledDipole-jl)
[![Coverage](https://codecov.io/gh/NoelAraujo/CoupledDipole.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/NoelAraujo/CoupledDipole.jl)
[![Coverage](https://coveralls.io/repos/github/NoelAraujo/CoupledDipole.jl/badge.svg?branch=master)](https://coveralls.io/github/NoelAraujo/CoupledDipole.jl?branch=master)


### Still Experimental: Be aware that syntax may change

You create your simulation with 3 building blocks `Atom`, `Laser`, `Atom`, `LinearOpitcs`/`NonLinearOptics`.

# `Atom`

Example:

```
using CoupledDipoles
cubic_cloud = Atom(Cube(), 100, 5)

N, ρ = 100, 0.02
spherical_cloud = Atom(Sphere(), sphere_inputs(N, ρ)...)

#note: cube_inputs(N, ρ) function also exists
```

The first line creates `N=100` atoms inside a `Cube` of length of side `sizes=5`.  
The second example creates atoms inside a sphere, but we specify the atomic density to  an auxiliary functions to convert density into the correct size.


Your atoms will be stored inside a Julia `struct`, with 4 components `shape`, `r`, `N`, `sizes`.

- `shape` is the cloud distribution in space. Valid arguments are `Cube()`, `Sphere()`  - the distributions are homogeneous.
- `r` is a matrix with atomic positions. First **row** is x-component, second row the y-component, and if exists, third row is the z-component.
- `N` is the number of atoms.
- `sizes` is the size of the system, and is a variable context dependent. For example, the `size` of a sphere, is its radius.

To see the atoms:
```
using Plots
my_atoms = cubic_cloud.r
x,y,z = my_atoms[1,:], my_atoms[2,:], my_atoms[3,:]
scatter(x,y,z)
```


# `Laser`
```
s = 1e-5 # laser saturation
Δ = 0.0  # laser detuning (must be a Float64 number (Integers not allowed))

w₀ = size(spherical_cloud) # laser waist for (Gaussian lasers)
gaussian_laser = Laser( Gaussian3D(w₀), s, Δ)

k = [0,0,1] # direction of propagation
plane_wave_laser = Laser(PlaneWave3D(k), s, Δ)
```

- `pump` is the laser caracteristics, `Gaussian3D` or `PlaneWave3D`
- `s` is the laser saturation. For linear models, avoid values `s > 1e-2`.
- `Δ` is the laser detuning (difference of laser and atomic frequencies)

All of this laser features are **NOT** mutable. You need to create another laser with the new specs.

# `LinearOptics` and `NonLinearOptics`
These struct clues atoms and laser though out all the package. You pass them over the functions, and multiple-dispatch apply the correct formula.

```
linear_problem = LinearOptics(Scalar(), cubic_cloud, gaussian_laser)
nonlinear_problem = NonLinearOptics(MeanField(), cubic_cloud, gaussian_laser)
```

The syntax is simple `LinearOptics(physics, atoms, laser)` or `NonLinearOptics(physics, atoms, laser)`.

The `Linear Physics` can be `Scalar()` or `Vectorial()`.

The `Non Linear Physics` supports only `MeanField()`.

# `Functions`
Here a list of some commom functions:
- `Atom` related: `cube_inputs`, `sphere_inputs`, `get_rₘᵢₙ`, `get_pairwise_matrix`, `get_interaction_matrix`.
- `Laser` related: `apply_laser_over_atoms`, `apply_laser_over_oneAtom`, `apply_laser_over_sensors`, `apply_laser_over_oneSensor`.
- `Get Physics` related: `get_intensities_over_sensors`, `get_intensity_over_an_angle`, `get_steady_state`, `time_evolution`, `get_transmission`, `get_spectrum`, `get_IPRs`, `get_PRs`, `get_localization_length`, `get_spatial_profile_single_mode`, `classify_modes`.


# `Transmission Example`
You can calculate the transmission spectrum (currently only for `Scalar` case) with the following code:
```
N, ρ = 1500, 0.2
spherical_cloud = Atom(Sphere(), sphere_inputs(N, ρ)...)

s = 1e-5
w₀ = size(spherical_cloud)
gaussian_laser = Laser( Gaussian3D(w₀), s, Δ)

Δ_range = range(-50, 50, length = 50)
T = zeros(length(Δ_range))

for idx ∈ 1:length(Δ_range)
    Δ = Δ_range[idx]
    println("computing Δ = $(Δ)")
    _laser   = Laser(Gaussian3D(w₀), s, Δ)
    _problem = LinearOptics(Scalar(), spherical_cloud, _laser)
    _βₙ      = get_steady_state(_problem)

    T[idx] = get_transmission(_problem, _βₙ)
end


plot(
    Δ_range,
    T,
    label = "",
    ylims = (0, 1.2),
    size = (800, 400),
    lw = 3,
    legend = :bottomright,
)
hline!([1], linestyle = :dash, c = :black, label = "")
xlabel!("Δ")
ylabel!("Transmission")
```

![GitHub Logo](/images/FigTransmission.png)

# `Publication Example`
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