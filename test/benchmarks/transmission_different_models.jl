@time using CoupledDipoles
using Plots, ProgressMeter, Random
plotly()

# cloud settings
N = 100
R = 2π
# laser settings
w₀ = 4π
polarization = [1, im, 0] ./ √2

# transmission settings
Δ_range = -20:1:20 # range(-20, 20; length = 40)
regime = :near_field
rtol = exp10(-7)
m = 100

for ρ = [0.01, 0.2], s =[1e-6, 1e-1]

	Random.seed!(999)
	cloud = Atom(CoupledDipoles.Cylinder(), cylinder_inputs(N, ρ; R = R)...)
	h = cloud.sizes.h

	T = zeros(length(Δ_range), 3)
	p = Progress(length(Δ_range); showspeed = true)
	Threads.@threads :dynamic for idx in eachindex(Δ_range)
		Δ = Δ_range[idx]

		_laser = Laser(Gaussian3D(w₀), s, Δ; polarization = polarization)
		sca_problem = LinearOptics(Scalar(), cloud, _laser)
		sca_βₙ = steady_state(sca_problem)
		T[idx, 1] = transmission(sca_problem, sca_βₙ; regime = regime, rtol = rtol)
		ProgressMeter.next!(p)
	end

	p = Progress(length(Δ_range); showspeed = true)
	Threads.@threads :dynamic for idx in eachindex(Δ_range)
		Δ = Δ_range[idx]

		_laser = Laser(Gaussian3D(w₀), s, Δ; polarization = polarization)
		vec_problem = LinearOptics(Vectorial(), cloud, _laser)
		vec_βₙ = steady_state(vec_problem)
		T[idx, 2] = transmission(vec_problem, vec_βₙ; regime = regime, rtol=rtol)
		ProgressMeter.next!(p)
	end

	p = Progress(length(Δ_range); showspeed = true)
	for idx in eachindex(Δ_range)
		Δ = Δ_range[idx]

		_laser = Laser(Gaussian3D(w₀), s, Δ; polarization = polarization)
		mef_problem = NonLinearOptics(MeanField(), cloud, _laser)
		mef_βₙ = steady_state(mef_problem; m=m)
		T[idx, 3] = transmission(mef_problem, mef_βₙ; regime = regime, rtol=rtol)
		ProgressMeter.next!(p)
	end


	plot(
		ylims = (0, 1.3),
		size = (800, 400),
		legend = :bottomright,
	)
	plot!(Δ_range, T[:, 1]; label = "Scalar", lw = 4)
	plot!(Δ_range, T[:, 2]; label = "Vectorial", lw = 3, linestyle = :dash)
	plot!(Δ_range, T[:, 3]; label = "Mean Field", lw = 3, c = :black, linestyle = :dot)

	hline!([1]; linestyle = :dash, c = :black, label = "")
	xlabel!("Δ")
	ylabel!("Transmission")
	display(title!("Cylinder : N=$(N), ρ=$(round(ρ,digits=3)), R=$(round(R,digits=2)),
		h=$(round(h,digits=2)), w₀=$(round(w₀,digits=2)), s=$(s)"))
end


