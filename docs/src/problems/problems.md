# Problems

Your simulation will only work if you use the right problem data type, provided by the constructors `LinearOptics` or `NonLinearOptics`.

For the `LinearOptics` you have the `Scalar` and `Vectorial` models, whereas the `NonLinearOptics` have `MeanField` and `PairCorrelation`. All of them deal with the 3 dimensions - different dimensions have different equations which are waiting for someone (not the author) to write them in the code engines.

Here all the possible configurations.

```julia
# settings
N, kR = 15, 32.4
w₀, s, Δ = 4π, 1e-5, 0.3   
atoms = Atom(CoupledDipoles.Sphere(gaussian=true), N, kL; r_min=0.0)
laser = Laser(Gaussian3D(w₀), s, Δ)

# LinearOptics options
prob_scalar = LinearOptics(Scalar(), atoms, laser)
prob_vectorial = LinearOptics(Vectorial(), atoms, laser)

# NonLinearOptics options
prob_meanfield = NonLinearOptics(MeanField(), atoms, laser)
prob_pair = NonLinearOptics(PairCorrelation(), atoms, laser)


# (do something with your problem)
# (...)
```

---

```@docs
LinearOptics
```

```@docs
NonLinearOptics
```