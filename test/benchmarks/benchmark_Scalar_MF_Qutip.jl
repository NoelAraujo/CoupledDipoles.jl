using CoupledDipole, Revise, Plots, DelimitedFiles

### --- INITIAL CONDITIONS ARE CHANGED TO AGREE WITH QUTIP SAVED DATA ---
function CoupledDipole.get_initial_conditions(problem::SimulationScalar)
    return zeros(ComplexF64, problem.atoms.N)
end

function CoupledDipole.get_initial_conditions(problem::SimulationMeanField)   
    β₀ = zeros(ComplexF64, problem.atoms.N)
    z₀ = -ones(problem.atoms.N)
    u₀ = vcat(β₀, z₀)
    return u₀
end


### ------------ FIXED INPUTS ---------------------
r = [ -1.01015   -0.455099     0.543082;
      -2.96935   -0.0781967   -0.00764304;
      -0.714778  -0.248556     0.683987;
      -1.20512    0.376068    -1.85059;
      -1.89765    0.84797      0.248158;
      0.664911   0.00503915  -0.177228]

atoms = Cube(r, 10);

Δ = -2.0
s = 0.0022
laser = PlaneWave_3D(:z, s, Δ);

### ------------ CREATE PROBLEM AND EVOLVE ---------------------
simulationScalar = ScalarProblem(atoms, laser);
simulationMeanField = MeanFieldProblem(atoms, laser);

aS = time_evolution(simulationScalar, time_max=100);
aMF = time_evolution(simulationMeanField, time_max=100);


### ------------ GET INTENSITIES ---------------------
iS = zeros( length(aS.time_array) )
for t ∈ eachindex(iS)
    iS[t] = get_scattered_intensity(simulationScalar, aS.βₜ[t], deg2rad(35))
end


iMF = zeros( length(aMF.time_array) )
for t ∈ eachindex(iMF)
    iMF[t] = get_scattered_intensity(simulationMeanField, vcat(aMF.βₜ[t], aMF.zₜ[t]), deg2rad(35))
end


matrixData = readdlm(pwd()*"/test/benchmarks/Field_qutip.txt", '\t', Float64, '\n')
time_Qutip = matrixData[:,1]
iQutip = matrixData[:,2]

### ------------ STEADY STATES -----------------------
ss_S = get_steady_state(simulationScalar)
ss_MF = get_steady_state(simulationMeanField; time_max=1000)

i_ss_S = get_scattered_intensity(simulationScalar, ss_S, deg2rad(35))
i_ss_MF = get_scattered_intensity(simulationMeanField, ss_MF , deg2rad(35))

### ------------ PLOTS ---------------------
plot(aS.time_array, iS, label="Scalar", lw=5, ylabel="intensity", xlabel="time")
plot!(aMF.time_array, iMF, label="Mean Field", lw=4, linestyle=:dash, size=(800,600))
plot!( time_Qutip, iQutip, label="Qutip", lw=4, linestyle=:dot  )
plot!(guidefont=17, tickfont=15, legendfontsize=15, legend=:bottomright)

hline!([i_ss_S], label="i_ss_S", lw=5, linestyle=:dash)
hline!([i_ss_MF], label="i_ss_MF", lw=5, linestyle=:dot)