## Steady State

### `LinearOptics`

The time evolution of `LinearOptics` follows the equation:


$$\frac{d\vec{\beta}}{dt} = G\vec{\beta} + \vec{\Omega}$$

which has a well-defined formal solution:

$$\vec{\beta}_{s} = -G^{-1}\vec{\Omega}$$

#### `Scalar`
In the `Scalar` case, the soltuion is a `Vector` solution. 

```@example ss_scalar
using CoupledDipoles
r =[1 2 0;
    1 0 1.0]
atoms = Atom(Cube(), Array(transpose(r)), 10)
laser = Laser(PlaneWave3D(), 1e-6, 1.0)
problem = LinearOptics(Scalar(), atoms, laser)
βₛ = steady_state(problem) # N-array
```

#### `Vectorial`
In `Vectorial` case, one gets a `Matrix`:
- each row represents the `x,y,z`-components of the polarization
- each colums correspond to different atoms

```@example ss_scalar
# (...) same as Scalar case
problem = LinearOptics(Vectorial(), atoms, laser)
βₛ = steady_state(problem) # 3xN-matrix
```


### `NonLinearOptics`
`NonLinearOptics` does not have a formal solution, therefore `steady_state` have two approaches


- Use `NewtonRaphson` methods (default)
- If `ode_solver=true`, use the `time_evolution` function over the period `tspan = (0, 250)` and return the final state


#### `MeanField`
The solution is a `Vector` of size `2N`, corresponding to [$\langle \sigma^- \rangle$ $\langle \sigma^z \rangle$].

```@example ss_scalar
# (...) same as Scalar case
problem = NonLinearOptics(MeanField(), atoms, laser)
βₛ = steady_state(problem) # default with NewtonRaphson

βₛ = steady_state(problem; ode_solver=true) # bruteforce time evoltuion
```


#### `PairCorrelation`
The solution is a `Vector` of size `2N + 4N^2`, corresponding to [$\sigma^{-}, \sigma^{z}, \sigma^{z}\sigma^{-}, \sigma^{+}\sigma^{-}, \sigma^{-}\sigma^{-}, \sigma^{z}\sigma^{z}$].

```@example ss_scalar
# (...) same as Scalar case
problem = NonLinearOptics(PairCorrelation(), atoms, laser)
βₛ = steady_state(problem) # default with NewtonRaphson

βₛ = steady_state(problem; ode_solver=true) # bruteforce time evoltuion
```

---

```@docs
steady_state
```