using CoupledDipoles, Revise, Plots, DelimitedFiles

### --- INITIAL CONDITIONS ARE CHANGED TO AGREE WITH QUTIP SAVED DATA ---
function get_initial_condition(problem::LinearOptics{Scalar})
    return zeros(ComplexF64, problem.atoms.N)
end

function get_initial_condition(problem::NonLinearOptics{MeanField})   
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

atoms = Atom(Cube(), Array(transpose(r)), 10);

Δ = -2.0
s = 0.0022
laser = Laser(PlaneWave3D([0,0,1]), s, Δ);

### ------------ CREATE PROBLEM AND EVOLVE ---------------------
simulationScalar =  LinearOptics(Scalar(), atoms, laser)
simulationMeanField = NonLinearOptics(MeanField(), atoms, laser)

u0_scalar    = get_initial_condition(simulationScalar)
u0_meanfield = get_initial_condition(simulationMeanField)
tspan = (0.0,100)


aS = time_evolution(simulationScalar, u0_scalar, tspan);
aMF = time_evolution(simulationMeanField, u0_meanfield, tspan);


### ------------ GET INTENSITIES ---------------------
iS = zeros( length(aS.t) );
for t ∈ eachindex(iS)
    iS[t] = get_intensity_over_an_angle(simulationScalar, aS.u[t], deg2rad(35))
end


iMF = zeros( length(aMF.t) );
for t ∈ eachindex(iMF)
    iMF[t] = get_intensity_over_an_angle(simulationMeanField, aMF.u[t], deg2rad(35))
end


matrixData = readdlm(pwd()*"/test/benchmarks/Field_qutip.txt", '\t', Float64, '\n')
time_Qutip = matrixData[:,1]
iQutip = matrixData[:,2]

### ------------ STEADY STATES -----------------------
ss_S = steady_state(simulationScalar)
ss_MF = steady_state(simulationMeanField)

i_ss_S = get_intensity_over_an_angle(simulationScalar, ss_S, deg2rad(35))
i_ss_MF =get_intensity_over_an_angle(simulationMeanField, ss_MF , deg2rad(35))

### ------------ PLOTS ---------------------
plot(aS.t, iS, label="Scalar", lw=5, ylabel="intensity", xlabel="time")
plot!(aMF.t, iMF, label="Mean Field", lw=4, linestyle=:dash, size=(800,600))
plot!( time_Qutip, iQutip, label="Qutip", lw=4, linestyle=:dot  )
plot!(guidefont=17, tickfont=15, legendfontsize=15, legend=:bottomright)

hline!([i_ss_S], label="i_ss_S", lw=5, linestyle=:dash)
hline!([i_ss_MF], label="i_ss_MF", lw=5, linestyle=:dot)