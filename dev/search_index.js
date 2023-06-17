var documenterSearchIndex = {"docs":
[{"location":"steady_state/#Steady-State","page":"Steady State","title":"Steady State","text":"","category":"section"},{"location":"steady_state/","page":"Steady State","title":"Steady State","text":"LinearOptics has a well defined analytical solution through a linear system of equations. In Scalar case, one get an Array solution, and in Vectorial case, one gets a Matrix, where each collum represents the x,y,z-components.","category":"page"},{"location":"steady_state/","page":"Steady State","title":"Steady State","text":"using CoupledDipoles\nr =[1 2 0;\n    1 0 1.0]\natoms = Atom(Cube(), Array(transpose(r)), 10)\nlaser = Laser(PlaneWave3D(), 1e-6, 1.0)\nproblem = LinearOptics(Scalar(), atoms, laser)\nβₛₛ = steady_state(problem) # N-array\n\nproblem = LinearOptics(Vectorial(), atoms, laser)\nβₛₛ = steady_state(problem) # 3xN-matrix","category":"page"},{"location":"steady_state/","page":"Steady State","title":"Steady State","text":"NonLinearOptics has not well defined solution, therefore steady_state apply time_evolution function over the period tspan = (0, 500) and return the final state.","category":"page"},{"location":"steady_state/","page":"Steady State","title":"Steady State","text":"The solution is an Array of size 2N, corresponding to [langle sigma^- rangle langle sigma^z rangle].","category":"page"},{"location":"steady_state/","page":"Steady State","title":"Steady State","text":"using CoupledDipoles\nr =[1 2 0;\n    1 0 1.0]\natoms = Atom(Cube(), Array(transpose(r)), 10)\nlaser = Laser(PlaneWave3D(), 1e-6, 1.0)\nproblem = NonLinearOptics(MeanField(), atoms, laser)\nβₛₛ = steady_state(problem) # 2N-array","category":"page"},{"location":"steady_state/","page":"Steady State","title":"Steady State","text":"steady_state","category":"page"},{"location":"steady_state/#CoupledDipoles.steady_state","page":"Steady State","title":"CoupledDipoles.steady_state","text":"steady_state(problem::LinearOptics{Scalar})\n\nSolve x=G\\Ω, with default interaction_matrix and laser_field.\n\n\n\n\n\nsteady_state(problem::LinearOptics{Vectorial})\n\nSolve x=-G\\Ω, with default interaction_matrix and laser_field. The solution x is reshaped as a 3xN matrix.\n\n\n\n\n\nsteady_state(problem::NonLinearOptics{MeanField})\n\n\n\n\n\n","category":"function"},{"location":"linear_regression/linear_regression_methods/#Linear-Regression-with-Noise-Data","page":"Fitting","title":"Linear Regression with Noise Data","text":"","category":"section"},{"location":"linear_regression/linear_regression_methods/","page":"Fitting","title":"Fitting","text":"The localization length comes from noise data, and one needs to pre process the data before extracting any analysis. In this section, we show how we use the package LinRegOutliers.jl for our application.","category":"page"},{"location":"linear_regression/linear_regression_methods/","page":"Fitting","title":"Fitting","text":"Localization length, xi, comes from the assumption that the data follows","category":"page"},{"location":"linear_regression/linear_regression_methods/","page":"Fitting","title":"Fitting","text":"y = A_0e^-xxi","category":"page"},{"location":"linear_regression/linear_regression_methods/","page":"Fitting","title":"Fitting","text":"However, the package LinRegOutliers.jl works best with linear functions, thus, we define","category":"page"},{"location":"linear_regression/linear_regression_methods/","page":"Fitting","title":"Fitting","text":"Y = ln(A_0) -xxi","category":"page"},{"location":"linear_regression/linear_regression_methods/","page":"Fitting","title":"Fitting","text":"or","category":"page"},{"location":"linear_regression/linear_regression_methods/","page":"Fitting","title":"Fitting","text":"Y = A + Bx","category":"page"},{"location":"linear_regression/linear_regression_methods/","page":"Fitting","title":"Fitting","text":"Our goal then becomes to identify xi = -1B.","category":"page"},{"location":"linear_regression/linear_regression_methods/#Exponential-with-Noise","page":"Fitting","title":"Exponential with Noise","text":"","category":"section"},{"location":"linear_regression/linear_regression_methods/","page":"Fitting","title":"Fitting","text":"We will create some data and check for the methods avaiable on LinRegOutliers.jl that leads to best results. Our test function will be","category":"page"},{"location":"linear_regression/linear_regression_methods/","page":"Fitting","title":"Fitting","text":"y = 2234e^-x1345","category":"page"},{"location":"linear_regression/linear_regression_methods/","page":"Fitting","title":"Fitting","text":"and we add white noise at each point of the domain x","category":"page"},{"location":"linear_regression/linear_regression_methods/","page":"Fitting","title":"Fitting","text":"using Random\nRandom.seed!(1234)\n\nx = range(0, 5, length = 40)\ny = log.(2.234exp.(-x/1.345) .+ rand(length(x)) ./ 10)","category":"page"},{"location":"linear_regression/linear_regression_methods/","page":"Fitting","title":"Fitting","text":"Now, we will create a functin that returns the coeficientes of the linear fit. For that, we specify our @formula, the Dataframe with our data, and use one of the available methods on LinRegOutliers.jldocs.","category":"page"},{"location":"linear_regression/linear_regression_methods/","page":"Fitting","title":"Fitting","text":"using LinRegOutliers, DataFrames\nfunction linear_fit_robust(x, y; regression_method = lta)\n\treg = createRegressionSetting(@formula(y ~ x), DataFrame([:x => x, :y => y]))\n\tresult = regression_method(reg)\n\n\tA, B = result[\"betas\"]\n\ty_fit = A .+ (B .* x)\n\n\treturn A, B, y_fit\nend\n## test\nA, B, y_fit = linear_fit_robust(x, y)\nξ = -1/B","category":"page"},{"location":"linear_regression/linear_regression_methods/","page":"Fitting","title":"Fitting","text":"Now, we can test and visualize the fitting result","category":"page"},{"location":"linear_regression/linear_regression_methods/","page":"Fitting","title":"Fitting","text":"using CairoMakie\nscatter(x, y, axis=(xlabel=\"x\", ylabel=\"y\"))\nlines!(x, y_fit, color = Cycled(2), linewidth = 5, label=\"ξ = $( round(ξ, digits=3))\")\naxislegend()\ncurrent_figure()","category":"page"},{"location":"linear_regression/linear_regression_methods/","page":"Fitting","title":"Fitting","text":"(Image: )","category":"page"},{"location":"linear_regression/linear_regression_methods/","page":"Fitting","title":"Fitting","text":"","category":"page"},{"location":"linear_regression/linear_regression_methods/#Methods-Comparison","page":"Fitting","title":"Methods Comparison","text":"","category":"section"},{"location":"linear_regression/linear_regression_methods/","page":"Fitting","title":"Fitting","text":"The exact metho depends on the level of noise on the data. On this specific example, the best regression method was the asm2000.","category":"page"},{"location":"linear_regression/linear_regression_methods/","page":"Fitting","title":"Fitting","text":"begin\n    fig = Figure()\n    ax = Axis(fig[1, 1], xlabel = \"x\", ylabel=\"y\")\n    scatter!(ax, x, y)\n\tmethods_list = [hs93, ks89, smr98, lms, lts, bch, py95, satman2013, satman2015,\n\t\tasm2000, lad, lta, galts, imon2005, ccf, cm97, quantileregression]\n    comparisons = []\n\tfor method in methods_list\n\t\ta, b, y_fit = linear_fit_robust(x, y; regression_method = method)\n        lines!(ax, x, y_fit, linewidth = 1, label = \"$(method)\")\n\n        comparison =  (-1/b)/1.345\n        push!(comparisons, comparison)\n\tend\n    value, idx  = findmin(comparisons)\n    println(\"$(methods_list[idx]) had the best fitting, with ξ = $(value*1.345)\")\n\n    Legend(fig[1, 2], ax)\n\tcurrent_figure()\nend","category":"page"},{"location":"linear_regression/linear_regression_methods/","page":"Fitting","title":"Fitting","text":"(Image: )","category":"page"},{"location":"scattering/#Scattering","page":"Scattering","title":"Scattering","text":"","category":"section"},{"location":"scattering/","page":"Scattering","title":"Scattering","text":"The electric field measured at sensor positioned at mathbfR = Rhatn results from the sum of all  dipoles located at mathbfr_j and their respective states beta_j. The exact formula depends on the model and regime.","category":"page"},{"location":"scattering/","page":"Scattering","title":"Scattering","text":"","category":"page"},{"location":"scattering/#Electric-Field","page":"Scattering","title":"Electric Field","text":"","category":"section"},{"location":"scattering/#Scalar","page":"Scattering","title":"Scalar","text":"","category":"section"},{"location":"scattering/","page":"Scattering","title":"Scattering","text":"regime = :near_field","category":"page"},{"location":"scattering/","page":"Scattering","title":"Scattering","text":"E_sc(mathbfR t) = +ifracGamma2sum_j frace^ ik_0mathbfR - mathbfr_j k_0mathbfR - mathbfr_jbeta_j(t)","category":"page"},{"location":"scattering/","page":"Scattering","title":"Scattering","text":"regime = :far_field","category":"page"},{"location":"scattering/","page":"Scattering","title":"Scattering","text":"E_sc(mathbfR t) approx +ifracGamma2 frace^ ik_0R k_0Rsum_j exp( -ik_0hatn cdot mathbfr_j )beta_j(t)","category":"page"},{"location":"scattering/#Vectorial","page":"Scattering","title":"Vectorial","text":"","category":"section"},{"location":"scattering/","page":"Scattering","title":"Scattering","text":"regime = :near_field","category":"page"},{"location":"scattering/","page":"Scattering","title":"Scattering","text":"E_vec(mathbfR t) = -ifracGamma2sum_jsum_etaG_mueta(mathbfR-mathbfr_j)beta_j^eta(t)","category":"page"},{"location":"scattering/","page":"Scattering","title":"Scattering","text":"regime = :far_field","category":"page"},{"location":"scattering/","page":"Scattering","title":"Scattering","text":"E^mu_vec(mathbfRt) approx -ifracGamma2cdotfrac32 frace^ik_0Rk_0Rsum_jsum_eta(delta_mu eta - hatn_muhatn_eta^*)exp(-ik_0hatmathbfncdotmathbfr_j)beta_j^eta(t)","category":"page"},{"location":"scattering/#Mean-Field","page":"Scattering","title":"Mean Field","text":"","category":"section"},{"location":"scattering/","page":"Scattering","title":"Scattering","text":"regime = :near_field|:far_field","category":"page"},{"location":"scattering/","page":"Scattering","title":"Scattering","text":"mathbfE_mf = mathbfE_sc","category":"page"},{"location":"scattering/","page":"Scattering","title":"Scattering","text":"","category":"page"},{"location":"scattering/#Intensity","page":"Scattering","title":"Intensity","text":"","category":"section"},{"location":"scattering/#Scalar-2","page":"Scattering","title":"Scalar","text":"","category":"section"},{"location":"scattering/","page":"Scattering","title":"Scattering","text":"regime = :near_field|:far_field","category":"page"},{"location":"scattering/","page":"Scattering","title":"Scattering","text":"I_sc(vecRt) = mathbfE_sc^2","category":"page"},{"location":"scattering/#Vectorial-2","page":"Scattering","title":"Vectorial","text":"","category":"section"},{"location":"scattering/","page":"Scattering","title":"Scattering","text":"regime = :near_field|:far_field","category":"page"},{"location":"scattering/","page":"Scattering","title":"Scattering","text":"I_vec(hatnt) = mathbfE_vec^2 = sum_muE^mu_vec^2","category":"page"},{"location":"scattering/#Mean-Field-2","page":"Scattering","title":"Mean Field","text":"","category":"section"},{"location":"scattering/","page":"Scattering","title":"Scattering","text":"regime = :near_field","category":"page"},{"location":"scattering/","page":"Scattering","title":"Scattering","text":"I_mf(mathbfRt) = I_sc(mathbfR t)\n+fracGamma^2(2k_0)^2 left   sum_j frac- beta_j^2    + frac1+langle sigma_j^z rangle 2R-r_j^2 right ","category":"page"},{"location":"scattering/","page":"Scattering","title":"Scattering","text":"regime = :far_field","category":"page"},{"location":"scattering/","page":"Scattering","title":"Scattering","text":"I_mf(mathbfRt) = I_sc(mathbfR t) + fracGamma^2(2k_0R)^2sum_j=1^N left ( -beta_j^2 + frac1 + langle sigma_j^z rangle 2right )","category":"page"},{"location":"scattering/","page":"Scattering","title":"Scattering","text":"","category":"page"},{"location":"scattering/#Functions","page":"Scattering","title":"Functions","text":"","category":"section"},{"location":"scattering/","page":"Scattering","title":"Scattering","text":"scattered_electric_field","category":"page"},{"location":"scattering/#CoupledDipoles.scattered_electric_field","page":"Scattering","title":"CoupledDipoles.scattered_electric_field","text":"scattered_electric_field(problem, atomic_states, sensor_positions; regime=:far_field)\n\nReturns a Matrix{ComplexF64} with value of the Eletric Laser + Electric Scattered from atoms\n\nproblem: LinearOptics or NonLinearOptics\natomic_states: β for Scalar/Vectorial Model, or [β,z] for Mean Field Model\nsensor_positions: matrix with measurement points\n\nNote:\n\nScalar problem returns a Matrix and not a Vector, to maintain consistency   with Vectorial problem that necessary returns a Matrix,   where each column has the [Ex, Ey, Ez] components of the field.\nAlso, even for single sensor, returns a Matrix of one element.\n\nExample\n\nusing CoupledDipoles, Random\nRandom.seed!(111)\nN = 5\nkR, kh = 1.0, 1.0\natoms = Atom(Cylinder(), N, kR, kh)\n\ns, Δ = 1e-5, 1.0\nlaser = Laser(PlaneWave3D(), s, Δ; polarization=[1,0,0])\n\nproblem_scalar = LinearOptics(Scalar(), atoms, laser)\nproblem_vectorial = LinearOptics(Vectorial(), atoms, laser)\n\natomic_states_scalar = steady_state(problem_scalar)\natomic_states_vectorial = steady_state(problem_vectorial)\n\n## 1 sensor\nRandom.seed!(222)\nnSensors = 1\nsensor = Matrix(rand(3, nSensors)) # '3' == sensor in position in 3D space\nscattered_electric_field(problem_scalar, atomic_states_scalar, sensor)\nscattered_electric_field(problem_vectorial, atomic_states_vectorial, sensor)\n\n\n## 10 sensors\nRandom.seed!(333)\nnSensors = 10\nsensor = rand(3, nSensors) # '3' == sensor in position in 3D space\nscattered_electric_field(problem_scalar, atomic_states_scalar, sensor)\nscattered_electric_field(problem_vectorial, atomic_states_vectorial, sensor)\n\n\n\n\n\n","category":"function"},{"location":"scattering/","page":"Scattering","title":"Scattering","text":"laser_and_scattered_intensity","category":"page"},{"location":"scattering/#CoupledDipoles.laser_and_scattered_intensity","page":"Scattering","title":"CoupledDipoles.laser_and_scattered_intensity","text":"laser_and_scattered_intensity(problem, atomic_states, sensor_positions; regime=:far_field)\n\nReturns a Vector{Float64} with value of the |Electric Laser + Electric Scattered|^2 from atoms\n\nproblem: LinearOptics or NonLinearOptics\natomic_states: β for Scalar/Vectorial Model, or [β,z] for Mean Field Model\nsensor_positions: matrix with measurement points\n\n\n\n\n\n","category":"function"},{"location":"scattering/","page":"Scattering","title":"Scattering","text":"scattered_intensity","category":"page"},{"location":"scattering/#CoupledDipoles.scattered_intensity","page":"Scattering","title":"CoupledDipoles.scattered_intensity","text":"scattered_intensity(problem, atomic_states, sensor_positions; regime=:far_field)\n\nReturns a Vector{Float64} with value of the |Electric Scattered|^2 from atoms\n\nproblem: LinearOptics or NonLinearOptics\natomic_states: β for Scalar/Vectorial Model, or [β,z] for Mean Field Model\nsensor_positions: matrix with measurement points\n\nExample\n\nusing CoupledDipoles, Random\nRandom.seed!(111)\nN = 5\nkR, kh = 1.0, 1.0\natoms = Atom(Cylinder(), N, kR, kh)\n\ns, Δ = 1e-5, 1.0\nlaser = Laser(PlaneWave3D(), s, Δ; polarization=[1,0,0])\n\nproblem_scalar = LinearOptics(Scalar(), atoms, laser)\nproblem_vectorial = LinearOptics(Vectorial(), atoms, laser)\n\natomic_states_scalar = steady_state(problem_scalar)\natomic_states_vectorial = steady_state(problem_vectorial)\n\n## 1 sensor\nRandom.seed!(222)\nnSensors = 1\nsensor = Matrix(rand(3, nSensors)) # '3' == sensor in position in 3D space\nscattered_intensity(problem_scalar, atomic_states_scalar, sensor)\nscattered_intensity(problem_vectorial, atomic_states_vectorial, sensor)\n\n\n## 10 sensors\nRandom.seed!(333)\nnSensors = 10\nsensor = rand(3, nSensors) # '3' == sensor in position in 3D space\nscattered_intensity(problem_scalar, atomic_states_scalar, sensor)\nscattered_intensity(problem_vectorial, atomic_states_vectorial, sensor)\n\n\n\n\n\n","category":"function"},{"location":"scattering/","page":"Scattering","title":"Scattering","text":"get_intensity_over_an_angle","category":"page"},{"location":"scattering/#CoupledDipoles.get_intensity_over_an_angle","page":"Scattering","title":"CoupledDipoles.get_intensity_over_an_angle","text":"get_intensity_over_an_angle(problem::LinearOptics{Scalar}, atoms_states::Vector{ComplexF64}, θ::Number; tol=exp10(-7.4))\n\nUsed for the single angle and single single state (most probably user case).\n\nExample:\n\nusing CoupledDipoles, Random\nRandom.seed!(111)\nN = 5\nkR, kh = 1.0, 1.0\natoms = Atom(Cylinder(), N, kR, kh)\n\ns, Δ = 1e-5, 1.0\nlaser = Laser(PlaneWave3D(), s, Δ; polarization=[1,0,0])\n\nproblem_scalar = LinearOptics(Scalar(), atoms, laser)\natomic_states_scalar = steady_state(problem_scalar)\n\nθ = deg2rad(48)\nget_intensity_over_an_angle(problem_scalar, atomic_states_scalar, θ)\n\n\n\n\n\nget_intensity_over_an_angle(problem::LinearOptics{Scalar}, atoms_states::Vector{Vector{ComplexF64}}, θ::Number; tol=exp10(-7.4), exact_solution=false)\n\nUsed for the single angle and different states (for example, the output of time_evolution).\n\nExample:\n\nusing CoupledDipoles, Random\nRandom.seed!(111)\nN = 5\nkR, kh = 1.0, 1.0\natoms = Atom(Cylinder(), N, kR, kh)\n\ns, Δ = 1e-5, 1.0\nlaser = Laser(PlaneWave3D(), s, Δ; polarization=[1,0,0])\n\nproblem_scalar = LinearOptics(Scalar(), atoms, laser)\nu0 = default_initial_condition(problem_scalar)\ntspan = (0.0, 10.0)\nsolutions = time_evolution(problem_scalar, u0, tspan)\nstates = solutions.u\n\nθ = deg2rad(48)\nget_intensity_over_an_angle(problem_scalar, states, θ)\n\n\n\n\n\n","category":"function"},{"location":"variances_angles/rayleigh_variance/#Deviations-from-Rayleigh's-law","page":"Intensity Statistics","title":"Deviations from Rayleigh's law","text":"","category":"section"},{"location":"variances_angles/rayleigh_variance/","page":"Intensity Statistics","title":"Intensity Statistics","text":"Our goal is to reproduce Figure 2 from the paper Cottier et all, 2019, where the authors studied the statistics of the scattered light, and found that the variance of the intensity distribution deviates from the expected Rayleigh's law.","category":"page"},{"location":"variances_angles/rayleigh_variance/","page":"Intensity Statistics","title":"Intensity Statistics","text":"With exception of the Step 1, the code is expected to run without any adjustments.","category":"page"},{"location":"variances_angles/rayleigh_variance/","page":"Intensity Statistics","title":"Intensity Statistics","text":"Step 1","category":"page"},{"location":"variances_angles/rayleigh_variance/","page":"Intensity Statistics","title":"Intensity Statistics","text":"Load the necessary packages. For me (the author), I prefer to execute all repetitions in parallel in my home made cluster. If you don't want parallel processing, just remove the process nodes.","category":"page"},{"location":"variances_angles/rayleigh_variance/","page":"Intensity Statistics","title":"Intensity Statistics","text":"@time begin\n    using CairoMakie, LinearAlgebra\n    using Statistics: mean, var\n    using StatsBase: fit, normalize, Histogram\nend\n\n@time begin\n    using Distributed\n\n    ## for local paralellism\n    addprocs(2; exeflags=`--project=$(Base.active_project()) --threads 4`, topology=:master_worker, enable_threaded_blas=true)\n\n    ## for remote machines\n    # ip_m2 = \"pc2@192.168.15.8\"\n    # nprocess_m2 = 2\n    # machine2 = [(ip_m2, nprocess_m2)]\n    # addprocs(\n    #     machine2;\n    #     topology = :master_worker,\n    #     exeflags = `--threads 6`,\n    #     enable_threaded_blas = true,\n    #     exename = \"/home/pc2/julia/julia-1.8.0-rc3/bin/julia\",\n    #     dir = \"/home/pc2\",\n    # )\n\n    # ip_m3 = \"pc3@192.168.15.7\"\n    # nprocess_m3 = 2\n    # machine3 = [(ip_m3, nprocess_m3)]\n    # addprocs(\n    #     machine3;\n    #     topology = :master_worker,\n    #     exeflags = `--threads 4`,\n    #     enable_threaded_blas = true,\n    #     exename = \"/home/pc3/julia_program/julia-1.8.0-rc3/bin/julia\",\n    #     dir = \"/home/pc3\",\n    # )\nend\n\n@time @everywhere begin\n    using CoupledDipoles\n    using ProgressMeter, Random\nend","category":"page"},{"location":"variances_angles/rayleigh_variance/","page":"Intensity Statistics","title":"Intensity Statistics","text":"Step 2 We use the exact configuration parameters from the paper. You will notice may Warning messages. This happen because a small laser waist leads to unreasonable results.","category":"page"},{"location":"variances_angles/rayleigh_variance/","page":"Intensity Statistics","title":"Intensity Statistics","text":"If you are studying Cottier's paper, note that the results from the paper are not accurate due to its small laser waist, even though they general paper's message is still correct.","category":"page"},{"location":"variances_angles/rayleigh_variance/","page":"Intensity Statistics","title":"Intensity Statistics","text":"### ------------ ATOMS SPECS ---------------------\nL = 32.4\nN = [684, 6066]\n\n\n### ------------ LASER SPECS ---------------------\nΔ = 1.0\ns = 1e-6\nw₀ = L / 4\n\n### ------------ SIMULATION SPECS ---------------------\nsensors = get_sensors_ring(; num_pts = 720, kR = 300, θ = 5π / 12)\nmaxRep = 15","category":"page"},{"location":"variances_angles/rayleigh_variance/","page":"Intensity Statistics","title":"Intensity Statistics","text":"Step 3","category":"page"},{"location":"variances_angles/rayleigh_variance/","page":"Intensity Statistics","title":"Intensity Statistics","text":"For each atom number N, create maxRep atomic configurations, compute their state states, and scattered light intensity. The normalization over the mean comes from the paper.","category":"page"},{"location":"variances_angles/rayleigh_variance/","page":"Intensity Statistics","title":"Intensity Statistics","text":"### -------- PRODUCE INTENSITIES -----------------\nall_intensities = map(N) do N\n    many_intensities = @showprogress pmap(1:maxRep) do rep\n        Random.seed!(1134 + rep)\n\n        atoms = Atom(Cube(), N, L)\n        laser = Laser(Gaussian3D(w₀), s, Δ)\n        simulation = LinearOptics(Scalar(), atoms, laser)\n\n        βₙ = steady_state(simulation)\n        intensities = scattered_intensity(simulation, βₙ, sensors; regime = :far_field)\n\n        intensities\n    end\n\n    many_intensities = reduce(vcat, many_intensities)\n    all_intensities_over_mean = many_intensities ./ mean(many_intensities)\n\n    all_intensities_over_mean\nend","category":"page"},{"location":"variances_angles/rayleigh_variance/","page":"Intensity Statistics","title":"Intensity Statistics","text":"Step 4","category":"page"},{"location":"variances_angles/rayleigh_variance/","page":"Intensity Statistics","title":"Intensity Statistics","text":"Instead of ploting histogram for each particle number, we are interested in the data from the histogram to display it in a scatter plot. Also, this is the moment to compute the variance of all intensities.","category":"page"},{"location":"variances_angles/rayleigh_variance/","page":"Intensity Statistics","title":"Intensity Statistics","text":"### ------------ CREATE HISTOGRAM ---------------------\nbins = 10.0 .^ range(log10(1e-6), log10(75); length = 30)\n\nxy_data = map(eachindex(N)) do n\n    h = fit(Histogram, all_intensities[n], bins)\n\n    h_norm = normalize(h; mode = :pdf)\n    bins_edges = collect(h_norm.edges[1])\n    bins_centers = [sqrt(bins_edges[i] * bins_edges[i+1]) for i = 1:(length(bins_edges)-1)]\n    variance = var(all_intensities[n])\n\n    # x_data_histogram, y_data_histogram, variance\n    (bins_centers, h_norm.weights, variance)\nend","category":"page"},{"location":"variances_angles/rayleigh_variance/","page":"Intensity Statistics","title":"Intensity Statistics","text":"Step 5","category":"page"},{"location":"variances_angles/rayleigh_variance/","page":"Intensity Statistics","title":"Intensity Statistics","text":"Overlay the Distribution Probability in a single figure.","category":"page"},{"location":"variances_angles/rayleigh_variance/","page":"Intensity Statistics","title":"Intensity Statistics","text":"The begin-end structure is just to facilitate the Figure development for the user. It is easier to just run the block at once at every little plot tweak, than select and run everything all the time.","category":"page"},{"location":"variances_angles/rayleigh_variance/","page":"Intensity Statistics","title":"Intensity Statistics","text":"### ------------ FINAL PLOT ---------------------\nbegin\n    fig = Figure(resolution = (800, 450))\n    ax = Axis(\n        fig[1, 1],\n        xlabel = \"Intensity\",\n        ylabel = \"Probability Distribution\",\n        title = \"\",\n        xlabelsize = 25,\n        ylabelsize = 25,\n        xticklabelsize = 20,\n        yticklabelsize = 20,\n        xscale = log10,\n        yscale = log10,\n    )\n\n    ## theoretical curve\n    x_ray = range(0.01, 50; step = 0.15)\n    y_ray = exp.(-x_ray)\n    lines!(ax, x_ray, y_ray, linestyle = :dash, label = \"Rayleigh\", color = :black, lw = 4)\n\n    for n = 1:2\n        x = xy_data[n][1]\n        y = xy_data[n][2]\n        v = xy_data[n][3] # variance\n        notNull = findall(y .> 0)\n        scatter!(\n            ax,\n            x[notNull],\n            y[notNull];\n            label = \"N=$(N[n]), Variance = $( round(v,digits=3 ))\",\n            markershape = :circle,\n            markersize = 20,\n        )\n    end\n    ylims!(1e-6, 10)\n    xlims!(1e-1, 100)\n    axislegend(position = :rt, labelsize = 20)\n    fig\nend","category":"page"},{"location":"variances_angles/rayleigh_variance/","page":"Intensity Statistics","title":"Intensity Statistics","text":"(Image: Rayleigh Deviation)","category":"page"},{"location":"create_atoms/atoms/#Creating-Atoms","page":"Atom","title":"Creating Atoms","text":"","category":"section"},{"location":"create_atoms/atoms/","page":"Atom","title":"Atom","text":"CoupledDipoles currently supports only 3D objects, that is, Sphere, Cube and Cylinder, even though it is possible to extend many functionalities to 1D and 2D in future developments. To create atoms, use the Atom constructor.","category":"page"},{"location":"create_atoms/atoms/","page":"Atom","title":"Atom","text":"using CoupledDipoles, CairoMakie, Random\nRandom.seed!(2354)\nnAtoms = 5000\n\n# CoupledDipoles and CairoMakie expoerts 'Sphere', therefore\n# one needs to be specific\nsphere_radius = 1.5\nsphere_cloud = Atom(CoupledDipoles.Sphere(), nAtoms, sphere_radius)\n\ncube_side = 1.0\ncube_cloud = Atom(Cube(), nAtoms, cube_side)\n\ncylinder_radius = 0.5\ncylinder_height = 2.0\ncylinder_cloud = Atom(Cylinder(), nAtoms, cylinder_radius, cylinder_height)","category":"page"},{"location":"create_atoms/atoms/","page":"Atom","title":"Atom","text":"The result should like something similar to","category":"page"},{"location":"create_atoms/atoms/","page":"Atom","title":"Atom","text":"fig = Figure(resolution = (800, 300))\nax_sphere = Axis3(fig[1:2, 1:2], aspect = (1, 1, 1))\nax_cube = Axis3(fig[1:2, 3:4], aspect = (1, 1, 1))\nax_cylinder = Axis3(fig[1:2, 5:6], aspect = (1, 1, 1))\n\nsx, sy, sz = sphere_cloud.r[1, :], sphere_cloud.r[2, :], sphere_cloud.r[3, :]\nscatter!(ax_sphere, sx, sy, sz, color = sz)\n\ncx, cy, cz = cube_cloud.r[1, :], cube_cloud.r[2, :], cube_cloud.r[3, :]\nscatter!(ax_cube, cx, cy, cz, color = cz)\n\ncyx, cyy, cyz = cylinder_cloud.r[1, :], cylinder_cloud.r[2, :], cylinder_cloud.r[3, :]\nscatter!(ax_cylinder, cyx, cyy, cyz, color = cyz)\n\nhidedecorations!(ax_sphere)\nhidedecorations!(ax_cube)\nhidedecorations!(ax_cylinder)\nfig\n\nsave(\"geometries.png\", fig)","category":"page"},{"location":"create_atoms/atoms/","page":"Atom","title":"Atom","text":"(Image: 3D examples)","category":"page"},{"location":"create_atoms/atoms/#User-defined-geometry","page":"Atom","title":"User defined geometry","text":"","category":"section"},{"location":"create_atoms/atoms/","page":"Atom","title":"Atom","text":"The Atom constructor holds the shape used for multiple-dispatch on the right physical equation, r is the matrix containg the atomic positions, N is the number of atom, and sizes is a generic radius for the system.","category":"page"},{"location":"create_atoms/atoms/","page":"Atom","title":"Atom","text":"struct Atom{T<:Dimension}\n    shape::T\n    r::Matrix{Float64}\n    N::Int64\n    sizes::Any\nend","category":"page"},{"location":"create_atoms/atoms/","page":"Atom","title":"Atom","text":"If a user wants to create their own atomic configuration, there are two constraints to consider:","category":"page"},{"location":"create_atoms/atoms/","page":"Atom","title":"Atom","text":"The shape field can be handled easily by choosing one of the available options: Sphere, Cube, or Cylinder.\nThe matrix containing the atom positions must have each Cartesian dimension represented by a row.","category":"page"},{"location":"create_atoms/atoms/","page":"Atom","title":"Atom","text":"Here's an example code snippet demonstrating the creation of a custom atomic configuration:","category":"page"},{"location":"create_atoms/atoms/","page":"Atom","title":"Atom","text":"# Define atom positions as separate arrays\natom_1 = [1, 1, 1]\natom_2 = [2, 2, 2]\natom_3 = [3, 3, 3]\natom_4 = [4, 4, 4]\n\n# Combine atom positions into a single matrix\nr = transpose(vcat(atom_1, atom_2, atom_3, atom_4))\n# Note: The `transpose` function is used to fulfill the matrix constraint.\n\n# Convert the transposed result to an actual matrix using `Array`\nr = Array(r)\n# Note: The `transpose` operation returns a non-matrix object, so we use `Array` to materialize it.\n\n# Create the `Atom` object with the chosen shape and custom positions\ndummy_dimension = 5\natoms = Atom(Cube(), r, dummy_dimension)","category":"page"},{"location":"create_atoms/atoms/","page":"Atom","title":"Atom","text":"Ensure that the atom positions are arranged correctly in the matrix, with each row representing a Cartesian dimension. The usage of transpose, followed by Array, is necessary to adhere to the matrix requirements.","category":"page"},{"location":"create_atoms/atoms/","page":"Atom","title":"Atom","text":"Atom","category":"page"},{"location":"create_atoms/atoms/#CoupledDipoles.Atom","page":"Atom","title":"CoupledDipoles.Atom","text":"Atom(geometry::Cube, N::Int64, kL::Union{Real,Integer}; r_min)\n\nArguments\n\ngeometry::Cube: The geometry of the atom object, which should be a Cube.\nN::Int64: The number of atoms.\nkL::Union{Real,Integer}: Cube's side length.\n\nKeyword Arguments\n\n:r_min: Optional keyword argument specifying the minimum distance between atoms.\n\nExample\n\natom = Atom(Cube(), 100, 5.0; r_min = 0.1)\n\n\n\n\n\nAtom(geometry::Cylinder, N::Int64, R::Union{Real,Integer}, h::Union{Real,Integer}; r_min)\n\nArguments\n\ngeometry::Cylinder: The geometry of the atom object, which should be a Cylinder.\nN::Int64: The number of atoms.\nR::Union{Real,Integer}: The radius of the cylinder.\nh::Union{Real,Integer}: The height of the cylinder.\n\nKeyword Arguments\n\n:r_min: Optional keyword argument specifying the minimum distance between atoms.\n\nExample\n\natom = Atom(Cylinder(), 100, 5.0, 10.0; r_min = 0.1)\n\n\n\n\n\nAtom(geometry::Sphere, N::Int64, kR::Union{Real,Integer}; kwargs...)\n\nArguments\n\ngeometry::Sphere: The geometry of the atom object, which should be a Sphere.\nfor Gaussian distribution, set gaussian=true\nN::Int64: The number of atoms.\nkR::Union{Real,Integer}: The radius for the sphere.\n\nKeyword Arguments\n\n:r_min: Optional keyword argument specifying the minimum distance between atoms.\n\nReturns\n\nAn Atom object with the specified geometry, atom positions, and sphere parameters.\n\nExample\n\natom_homogenous = Atom(Sphere(), 100, 5.0; r_min = 0.1)\natom_gaussian = Atom(Sphere(gaussian=true), 100, 5.0; r_min = 0.1)\n\n\n\n\n\n","category":"type"},{"location":"#CoupledDipoles.jl","page":"Home","title":"CoupledDipoles.jl","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"A Package In Development","category":"page"},{"location":"#**Installation**","page":"Home","title":"Installation","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"On the REPL type","category":"page"},{"location":"","page":"Home","title":"Home","text":"import Pkg\nPkg.add(url=\"https://github.com/NoelAraujo/CoupledDipoles.jl\")\n\n\njulia> using CoupledDipoles\n\n# pkg> test CoupledDipoles","category":"page"},{"location":"","page":"Home","title":"Home","text":"This package uses Bessel, because is faster than SpecialFunctions, but MAYBE you had to install it manually.","category":"page"},{"location":"","page":"Home","title":"Home","text":"import Pkg;\n\nPkg.add(url=\"https://github.com/JuliaMath/Bessels.jl\")\nPkg.add(url=\"https://github.com/NoelAraujo/CoupledDipoles.jl\")","category":"page"},{"location":"#Manual-Outline","page":"Home","title":"Manual Outline","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"dipole_example/single_atom_volume/#Single-Atom","page":"Single Atom","title":"Single Atom","text":"","category":"section"},{"location":"dipole_example/single_atom_volume/","page":"Single Atom","title":"Single Atom","text":"The special case of a single atom is the radiation pattern of a dipole, defined in electromagnetism books.","category":"page"},{"location":"dipole_example/single_atom_volume/","page":"Single Atom","title":"Single Atom","text":"To create a single atom system, set N=1 and any other function should function accordingly. Note that all functions will return a Matrix of 1 element. This was an intentional decision to make internal functions interoperate effectively - if all inputs were a matrix, there was no ambiguity about the number of particles.","category":"page"},{"location":"dipole_example/single_atom_volume/","page":"Single Atom","title":"Single Atom","text":"We need to use the Vectorial model to best visualize the radiation pattern. The code below is a minimal working example for checking the radiation pattern through a volume slice. The laser will be pointing in the negative x-direction, to make the visualization clearer.","category":"page"},{"location":"dipole_example/single_atom_volume/","page":"Single Atom","title":"Single Atom","text":"using CoupledDipoles, Random\n\n# cloud settings\nN = 1\nkR = 10\n\n# laser settings\nw₀ = 4π\ns = 1e-5\nΔ = 0.0\n\nRandom.seed!(2044)\nsingle_atom =Atom(CoupledDipoles.Cylinder(), N, kR, kR)\nlaser = Laser(Gaussian3D(w₀), s, Δ; direction=[-1,0,0], polarization=[0,0,1])\nproblem = LinearOptics(Vectorial(), single_atom, laser)\nβₙ = steady_state(problem)","category":"page"},{"location":"dipole_example/single_atom_volume/","page":"Single Atom","title":"Single Atom","text":"We had to apply :near_field regime manually to compute the intensity, because the default regime is :far_field. To evaulate the intensity in a certain spatial domain, we  use a list comprehension.","category":"page"},{"location":"dipole_example/single_atom_volume/","page":"Single Atom","title":"Single Atom","text":"x = LinRange(-100, 100, 100)\ny = LinRange(-100, 100, 100)\nz = LinRange(-100, 100, 100)\n\n## this plot only works for near field regime\n_vol = [laser_and_scattered_intensity(problem, βₙ, Matrix([X Y Z]');regime=:near_field)[1] for X ∈ x, Y ∈ y, Z ∈ z]\nlaserOn = log10.(_vol)\n\n_vol = [scattered_intensity(problem, βₙ, Matrix([X Y Z]');regime=:near_field)[1] for X ∈ x, Y ∈ y, Z ∈ z]\nlaserOff = log10.(_vol)","category":"page"},{"location":"dipole_example/single_atom_volume/","page":"Single Atom","title":"Single Atom","text":"The package Makie is not a dependence of CoupledDipoles. Please, install it.","category":"page"},{"location":"dipole_example/single_atom_volume/","page":"Single Atom","title":"Single Atom","text":"In the following, the figures represents the radiation in space for single atom and the color range was choosen ad hoc to higlight the expected dipole radiation pattern.","category":"page"},{"location":"dipole_example/single_atom_volume/","page":"Single Atom","title":"Single Atom","text":"using ColorSchemes\nusing WGLMakie\nWGLMakie.activate!()\n\n## Example 1\nlet\n    fig = Figure(resolution = (800, 600))\n    ax_on = Axis3(fig[1, 1], title = \"Dipole Radiation\", aspect=:data)\n\n    on_plt = volumeslices!(ax_on, x, y, z, laserOff,\n        colormap=cgrad( ColorSchemes.linear_kryw_0_100_c71_n256, rev=false),\n        colorrange=(-11, -9)\n        )\n    on_plt[:update_yz][](100)\n    on_plt[:update_xz][](50)\n    on_plt[:update_xy][](1)\n\n    # save(\"dipole_radiation.png\", fig, resolution = (800, 600))\n    fig\nend","category":"page"},{"location":"dipole_example/single_atom_volume/","page":"Single Atom","title":"Single Atom","text":"(Image: Dipole Radiation)","category":"page"},{"location":"dipole_example/single_atom_volume/","page":"Single Atom","title":"Single Atom","text":"## Example 2\nlet\n    fig = Figure(resolution = (900, 800), background_color=:transparent)\n    ax_on = Axis3(fig[1, 1], title = \"Laser On\", aspect=:data)\n    ax_off = Axis3(fig[1, 2], title = \"Laser Off\", aspect=:data)\n\n    on_plt = volumeslices!(ax_on, x, y, z, laserOn,\n        colormap=cgrad( ColorSchemes.linear_kryw_0_100_c71_n256, rev=false),\n        colorrange=(-11, -9)\n        )\n    on_plt[:update_yz][](100)\n    on_plt[:update_xz][](50)\n    on_plt[:update_xy][](1)\n\n    off_plt = volumeslices!(ax_off, x, y, z, laserOff,\n        colormap=cgrad( ColorSchemes.linear_kryw_0_100_c71_n256, rev=false),\n        colorrange=(-11, -9)\n        )\n    off_plt[:update_yz][](100)\n    off_plt[:update_xz][](50)\n    off_plt[:update_xy][](1)\n\n    cbar = Colorbar(fig, off_plt; label=\"log10( Intensity )\", flipaxis=false,  vertical = false, width = Relative(4/5),ticks=WilkinsonTicks(3))\n    fig[2, :] = cbar\n\n    # save(\"on_off_radiation.png\", fig, resolution = (800, 600))\n    fig\nend","category":"page"},{"location":"dipole_example/single_atom_volume/","page":"Single Atom","title":"Single Atom","text":"(Image: Laser On and Off)","category":"page"}]
}
