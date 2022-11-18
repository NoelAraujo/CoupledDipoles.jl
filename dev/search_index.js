var documenterSearchIndex = {"docs":
[{"location":"steady_state/#Steady-State","page":"Steady State","title":"Steady State","text":"","category":"section"},{"location":"steady_state/","page":"Steady State","title":"Steady State","text":"LinearOptics has a defined analytica solution through a linear system of equations. In Scalar case, one get an Array solution, and in Vectorial case, one gets a Matrix, where each collum represents the x,y,z-components.","category":"page"},{"location":"steady_state/","page":"Steady State","title":"Steady State","text":"using CoupledDipoles\nr =[1 2 0;\n    1 0 1.0]\natoms = Atom(Cube(), Array(transpose(r)), 10)\nlaser = Laser(PlaneWave3D(), 1e-6, 1.0)\nproblem = LinearOptics(Scalar(), atoms, laser)\nβₛₛ = steady_state(problem) # N-array\n\nproblem = LinearOptics(Vectorial(), atoms, laser)\nβₛₛ = steady_state(problem) # 3xN-matrix","category":"page"},{"location":"steady_state/","page":"Steady State","title":"Steady State","text":"NonLinearOptics has not well defined solution, therefore steady_state apply time_evolution function over the period tspan = (0, 500) and return the final state.","category":"page"},{"location":"steady_state/","page":"Steady State","title":"Steady State","text":"The solution is an Array of size 2N, corresponding to [langle sigma^- rangle langle sigma^z rangle].","category":"page"},{"location":"steady_state/","page":"Steady State","title":"Steady State","text":"using CoupledDipoles\nr =[1 2 0;\n    1 0 1.0]\natoms = Atom(Cube(), Array(transpose(r)), 10)\nlaser = Laser(PlaneWave3D(), 1e-6, 1.0)\nproblem = NonLinearOptics(MeanField(), atoms, laser)\nβₛₛ = steady_state(problem) # 2N-array","category":"page"},{"location":"steady_state/","page":"Steady State","title":"Steady State","text":"steady_state","category":"page"},{"location":"steady_state/#CoupledDipoles.steady_state","page":"Steady State","title":"CoupledDipoles.steady_state","text":"steady_state(problem::LinearOptics{Scalar})\n\nSolve x=G\\Ω, with default interaction_matrix and laser_field.\n\n\n\n\n\nsteady_state(problem::LinearOptics{Vectorial})\n\nSolve x=G\\Ω, with default interaction_matrix and laser_field. The solution x is reshaped as a 3xN matrix.\n\n\n\n\n\nsteady_state(problem::NonLinearOptics{MeanField})\n\nFor Non Linear Optics makes a time evolution to infinity (t=500Γ by default) and returns the state.\n\n\n\n\n\n","category":"function"},{"location":"#CoupledDipoles.jl","page":"Home","title":"CoupledDipoles.jl","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"A Package In Development","category":"page"},{"location":"#**Installation**","page":"Home","title":"Installation","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"On the REPL type","category":"page"},{"location":"","page":"Home","title":"Home","text":"pkg> add https://github.com/NoelAraujo/CoupledDipoles.jl\n\njulia> using CoupledDipoles\n\npkg> test CoupledDipoles","category":"page"},{"location":"","page":"Home","title":"Home","text":"This package uses Bessel, because is faster than SpecialFunctions, but MAYBE you had to install it manually.","category":"page"},{"location":"","page":"Home","title":"Home","text":"import Pkg;\n\nPkg.add(url=\"https://github.com/JuliaMath/Bessels.jl\")\nPkg.add(url=\"https://github.com/NoelAraujo/CoupledDipoles.jl\")","category":"page"},{"location":"#Manual-Outline","page":"Home","title":"Manual Outline","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"","category":"page"}]
}
