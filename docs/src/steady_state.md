## Steady State

`LinearOptics` has a well defined analytical solution through a linear system of equations. In `Scalar` case, one get an `Array` solution, and in `Vectorial` case, one gets a `Matrix`, where each collum represents the `x,y,z`-components.


```@example ss_scalar
using CoupledDipoles
r =[1 2 0;
    1 0 1.0]
atoms = Atom(Cube(), Array(transpose(r)), 10)
laser = Laser(PlaneWave3D(), 1e-6, 1.0)
problem = LinearOptics(Scalar(), atoms, laser)
βₛₛ = steady_state(problem) # N-array
```


```@example ss_scalar
# (...) same as Scalar case
problem = LinearOptics(Vectorial(), atoms, laser)
βₛₛ = steady_state(problem) # 3xN-matrix
```

`NonLinearOptics` has not well defined solution, therefore `steady_state` apply `time_evolution` function over the period `tspan = (0, 500)` and return the final state.

The solution is an `Array` of size `2N`, corresponding to [$\langle \sigma^- \rangle$ $\langle \sigma^z \rangle$].

```@example ss_scalar
# (...) same as Scalar case
problem = NonLinearOptics(MeanField(), atoms, laser)
βₛₛ = steady_state(problem) # 2N-array
```


```@docs
steady_state
```