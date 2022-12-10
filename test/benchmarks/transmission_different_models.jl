using CoupledDipoles
using Plots, ProgressMeter
plotly()

# cloud settings
N = 350
ρ = 0.5

cloud = Atom(CoupledDipoles.Cylinder(), cylinder_inputs(N, ρ; h=5π)...)
R = cloud.sizes.R

# laser settings
w₀ = 4π
s = 1e0

# transmission settings
Δ_range = range(-20, 20; length=50)
T = zeros(length(Δ_range), 3)

plot(
    ylims = (0, 2),
    size=(800, 400),
    legend=:bottomright,
)
p = Progress(length(Δ_range); showspeed=true)
Threads.@threads for idx in 1:length(Δ_range)
    Δ = Δ_range[idx]

    _laser = Laser(Gaussian3D(w₀), s, Δ; polarization=[1,im,0]./√2)
    sca_problem = LinearOptics(Scalar(), cloud, _laser)
    vec_problem = LinearOptics(Vectorial(), cloud, _laser)
    mef_problem = NonLinearOptics(MeanField(), cloud, _laser)

    sca_βₙ = steady_state(sca_problem)
    vec_βₙ = steady_state(vec_problem)
    mef_βₙ = steady_state(mef_problem)


    T[idx, 1] = transmission(sca_problem, sca_βₙ; regime=:near_field)
    T[idx, 2] = transmission(vec_problem, vec_βₙ; regime=:near_field)
    T[idx, 3] = transmission(mef_problem, mef_βₙ; regime=:near_field)
    ProgressMeter.next!(p)
end

plot!(Δ_range, T[:, 1]; label="Scalar", lw=4)
plot!(Δ_range, T[:, 2]; label="Vectorial", lw=3, linestyle=:dash)
plot!(Δ_range, T[:, 3]; label="Mean Field", lw=3, c=:black, linestyle=:dot)

hline!([1]; linestyle=:dash, c=:black, label="")
xlabel!("Δ")
ylabel!("Transmission")
display(title!("Cylinder : N=$(N), ρ=$(round(ρ,digits=3)), R=$(round(R,digits=2)), w₀=$(round(w₀,digits=2)), λ=$(round(2π,digits=2))"))