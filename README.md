# CoupledDipoles

Check Documentation

 [![][docs-master-img]][docs-master-url]

[docs-master-img]: https://img.shields.io/badge/docs-master-blue.svg
[docs-master-url]: https://noelaraujo.github.io/CoupledDipoles.jl/dev/

## What the software do?

`CoupledDipoles.jl` provides a set of fundamental tools to simulate the Coupled Dipoles model in cold atoms. While it includes some analysis tools, its primary focus is to provide a flexible package that can be extended to accommodate different equation models and atomic geometries.

### Theory and Background
For a detailed explanation of the theory please refer to the [Phd Thesis](https://doi.org/10.11606/T.76.2024.tde-26012024-114225). 

Our [documentation](https://noelaraujo.github.io/CoupledDipoles.jl/dev/) will focus on how to use the package, rather than justifying the physics. Note that the code implementation **is** the ground truth regargind the right set of equations, because it may fix errors not found at the time of the PhD defense.



## Key Features

Our tools have been extensively validated over the years and are highly optimized for performance:

- Minimized memory allocation to reduce computational overhead
- Utilization of the fastest packages and methods for core computations
- Internal parameters carefully selected to balance speed and accuracy


## Important Note

Please note that our package is distinct from [CoupledDipole.jl](https://nano-optics.ac.nz/CoupledDipole.jl/dev/theory/), which focuses on polarization-related properties of materials - it is a coincidence that the package lacking a simple `s` at the end of the name exists.

```
Our approach is specifically designed to study properties in the context of **Cold Atom** physics.
```
