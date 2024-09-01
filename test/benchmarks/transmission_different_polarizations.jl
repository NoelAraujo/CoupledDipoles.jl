@time begin
	using CoupledDipoles
	using Plots, Random
	plotly()	
end

# cloud settings
R =  17.7183099768
L = 8.177681527782540982229875
k = 1
a,b,c = 0.1/k, 0.2/k, 0.3/k

X = [a,a,-a,-a] .* R
Y = [b,-b,b,-b] .* R
Z = [c,-c,-c,c] .* L 
r = vcat(X', Y', Z') |> Array # matriz dos atomos

# laser settings
w₀ = R/3
s = 1e-5
polarization_circular_positivo = [1, +im, 0] ./ √2 #circular polarization
polarization_circular_negativo = [1, -im, 0] ./ √2 #circular polarization
polarization_linear = [1, 0, 0]  # x-direction

# transmission settings
Δ_range = range(-0.75, 0.75; length = 50)

begin

	dummy_dimension = R
	cloud = Atom(Cube(), r, dummy_dimension)

	
	T = zeros(length(Δ_range), 3)
	
	laser_lin = Laser(Gaussian3D(w₀), s, Δ_range[1]; polarization = polarization_linear)
    problem_lin = LinearOptics(Vectorial(), cloud, laser_lin)
    @time T[:, 1] = transmission(problem_lin, Δ_range)

    laser_cir_p = Laser(Gaussian3D(w₀), s, Δ_range[1]; polarization = polarization_circular_positivo)
    problem_cir_p = LinearOptics(Vectorial(), cloud, laser_cir_p)
    @time T[:, 2] = transmission(problem_cir_p, Δ_range)

	laser_cir_n = Laser(Gaussian3D(w₀), s, Δ_range[1]; polarization = polarization_circular_negativo)
    problem_cir_n = LinearOptics(Vectorial(), cloud, laser_cir_n)
    @time T[:, 3] = transmission(problem_cir_n, Δ_range)
	
	plot(
		# ylims = (0, 1.3),
		size = (800, 400),
		legend = :bottomright,
	)
	plot!(Δ_range, log10.(T[:, 1]); label = "Linear em X", lw = 4)
	plot!(Δ_range, log10.(T[:, 2]); label = "Circular Positivo", lw = 3, linestyle = :dot)
	plot!(Δ_range, log10.(T[:, 3]); label = "Circular Negativo", lw = 3, linestyle = :dash)
	
	# hline!([1]; linestyle = :dash, c = :black, label = "")
	xlabel!("Δ")
	title!("log10( Transmission )") |> display
	# display(title!("Title"))
end