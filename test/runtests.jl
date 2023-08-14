using CoupledDipoles
using Test
using Tullio

@testset begin
    @testset "Green Matrix 1 Atom" begin
        singleAtom = Matrix([0 -1 0]')
        atoms = Atom(Cube(), singleAtom, 1)

        s, Δ = 1e-6, 1.5
        laser = Laser(PlaneWave3D(), s, Δ; direction=[0,0,1])

        ## Scalar case
        problem = LinearOptics(Scalar(), atoms, laser)
        G = interaction_matrix(problem)
        G_expected = [im*1.5 - 1/2]
        @test all(G .≈ G_expected)

        ## Vectorial case
        problem = LinearOptics(Vectorial(), atoms, laser)
        G = interaction_matrix(problem)
        G_expected = fill(im*1.5 - 1/2, 3, 3)
        @test all(G .≈ G_expected)


        r = [1, 0.5, -1]
        Δ = 0
        Γ = 1

        s, Δ = 1e-4, 0.0
        laser = Laser(PlaneWave3D(), s, Δ)

        Ω₀ = CoupledDipoles.raby_frequency(laser)
        β_expected = -(0.5*im*Ω₀*cis(r[3]))/(Γ/2)

        singleAtom = Atom(Cube(), Array(Matrix(r')'), 1)
        problem = LinearOptics(Scalar(), singleAtom, laser)
        β_scalar = steady_state(problem)

        @test β_expected ≈ β_scalar[1][1]
    end
    @testset "Vectorial Kernel Matrix 3D" begin
        r_jm = [-3, √3, -2]
        G = CoupledDipoles._vectorial_3D_green_kernel(r_jm)
        k₀ = CoupledDipoles.k₀
        G_expected = (3cis(4k₀)/(8k₀))*(
            (1 + im/(4k₀) - 1/(16k₀^2)).*[1 0 0; 0 1 0; 0 0 1] +
            (-1 - 3im/(4k₀) + 3/(16k₀^2)).*[9 -3√3 6; -3√3 3 -2√3; 6 -2√3 4]./16
        )
        @test all(G .≈ G_expected)
    end

    @testset "Green Matrix 2 Atoms" begin
        r1 = [0,-1,0]
        r2 = [3,1,0]
        r = hcat(r1,r2)
        atoms = Atom(Cube(), Array(r), 1)
        Γ = CoupledDipoles.Γ
        s, Δ = 1e-6, 1.0
        laser = Laser(PlaneWave3D(), s, Δ; direction=[0,0,1])

        ## Scalar case
        problem = LinearOptics(Scalar(), atoms, laser)
        G = interaction_matrix(problem)

        G_expected = zeros(ComplexF64, 2, 2)
        G_expected[1,1] = im*1.0 - 1/2
        G_expected[1,2] = +im*(Γ / 2) *(cis(√13))/√13
        G_expected[2,1] = +im*(Γ / 2) *(cis(√13))/√13
        G_expected[2,2] = im*1.0 - 1/2

        @test all(G .≈ G_expected)

        ## Vectorial case
        problem = LinearOptics(Vectorial(), atoms, laser)
        G = interaction_matrix(problem)

        G_expected = zeros(ComplexF64, 6, 6)
        # r_12
        G_expected[1,1] = im*1.0 - 1/2
        G_expected[1,3] = 0.0
        G_expected[1,5] = 0.0
        G_expected[3,1] = 0.0
        G_expected[3,3] = im*1.0 - 1/2
        G_expected[3,5] = 0.0
        G_expected[5,1] = 0.0
        G_expected[5,3] = 0.0
        G_expected[5,5] = im*1.0 - 1/2

        # r_12
        G_expected[1,2] = -(Γ / 2) *(3cis(√13))/(2im√13)*( 4/13 - (14/13)*(im/√13 - 1/13))
        G_expected[1,4] = -(Γ / 2) *(3cis(√13))/(2im√13)*(-6/13 - (18/13)*(im/√13 - 1/13))
        G_expected[1,6] = 0.0
        G_expected[3,2] = -(Γ / 2) *(3cis(√13))/(2im√13)*(-6/13 - (18/13)*(im/√13 - 1/13))
        G_expected[3,4] = -(Γ / 2) *(3cis(√13))/(2im√13)*( 9/13 +  (1/13)*(im/√13 - 1/13))
        G_expected[3,6] = 0.0
        G_expected[5,2] = 0.0
        G_expected[5,4] = 0.0
        G_expected[5,6] = -(Γ / 2) *(3cis(√13))/(2im√13)*(    1 +        1*(im/√13 - 1/13))

        # r_21
        G_expected[2,1] = -(Γ / 2) *(3cis(√13))/(2im√13)*( 4/13 - (14/13)*(im/√13 - 1/13))
        G_expected[2,3] = -(Γ / 2) *(3cis(√13))/(2im√13)*(-6/13 - (18/13)*(im/√13 - 1/13))
        G_expected[2,5] = 0.0
        G_expected[4,1] = -(Γ / 2) *(3cis(√13))/(2im√13)*(-6/13 - (18/13)*(im/√13 - 1/13))
        G_expected[4,3] = -(Γ / 2) *(3cis(√13))/(2im√13)*( 9/13 +  (1/13)*(im/√13 - 1/13))
        G_expected[4,5] = 0.0
        G_expected[6,1] = 0.0
        G_expected[6,3] = 0.0
        G_expected[6,5] = -(Γ / 2) *(3cis(√13))/(2im√13)*(    1 +        1*(im/√13 - 1/13))

        # r_22
        G_expected[2,2] = im*1.0 - 1/2
        G_expected[2,4] = 0.0
        G_expected[2,6] = 0.0
        G_expected[4,2] = 0.0
        G_expected[4,4] = im*1.0 - 1/2
        G_expected[4,6] = 0.0
        G_expected[6,2] = 0.0
        G_expected[6,4] = 0.0
        G_expected[6,6] = im*1.0 - 1/2

        @test all(G .≈ G_expected)
    end


    @testset "3D scalar: 2 atoms" begin
        r =[1 2 0;
            1 0 1.0]
        atoms = Atom(Cube(), Array(transpose(r)), 10)
        laser = Laser(PlaneWave3D(), 1e-6, 1.0)
        problem = LinearOptics(Scalar(), atoms, laser)

        βₛₛ = steady_state(problem)
        G = interaction_matrix(problem)
        βₛₛ_expected = -CoupledDipoles.inverseMatrix2x2(G)*laser_field(laser, Array(transpose(r)))
        @test all(βₛₛ_expected .≈ βₛₛ)
    end
    @testset "3D scalar: 3 atoms" begin
        r =[ 0 0 0;
                1 1 1;
            -1 1 0.5]
        atoms = Atom(Cube(), Array(transpose(r)), 10)
        laser = Laser(PlaneWave3D(), 1e-6, 0.0; direction=[0,0,1])
        problem = LinearOptics(Scalar(), atoms, laser)

        R_nm = get_pairwise_matrix(atoms.r)
        R_nm_ans = [0 1.7320508075688 1.5;
                    1.73205080756888 0 2.06155281280883;
                    1.5 2.06155281280883 0]
        @test all(R_nm_ans .≈ R_nm)

        G = interaction_matrix(problem)
        G_expected = zeros(ComplexF64, 3, 3)
        G_expected[1,1] = -0.5
        G_expected[1,2] = +im*0.5 *(cis(√3))/√3
        G_expected[1,3] = +im*0.5 *(cis(√2.25))/√2.25

        G_expected[2,1] = +im*0.5 *(cis(√3))/√3
        G_expected[2,2] = -0.5
        G_expected[2,3] = +im*0.5 *(cis(√4.25))/√4.25

        G_expected[3,1] = +im*0.5 *(cis(√2.25))/√2.25
        G_expected[3,2] = +im*0.5 *(cis(√4.25))/√4.25
        G_expected[3,3] = -0.5

        @test all(G .≈ G_expected)

        βₛₛ = steady_state(problem)
        βₛₛ_expected = -CoupledDipoles.inverseMatrix3x3(G)*laser_field(laser, Array(transpose(r)))
        @test all(βₛₛ_expected .≈ βₛₛ)
    end

    @testset "Vectorial Polarization" begin
        w₀, s, Δ = 4π, 1e-5, 0.3
        laser = Laser(Gaussian3D(w₀), s, Δ; direction = [0,1,1], polarization=[1,0,0]) # OK
        @test laser.pump.w₀ == w₀
        @test laser.s == s
        @test laser.Δ == Δ
        @test laser.direction == [0,1,1]
        @test laser.polarization == [1,0,0]
    end
    @testset "Laser Over Points" begin
        N = 3; r = rand(3, N)
        w₀, s, Δ = 4π, 1e-5, 0.3

        atoms = Atom(Cube(), r, 1.5)
        laser = Laser(Gaussian3D(w₀), s, Δ)
        simulation = LinearOptics(    Scalar(),    atoms, laser)
        @test laser_field(laser, atoms) ≈ laser_field(laser, copy(atoms.r))

        sensor = [-2, 4, 6]

        laser_expected = -0.5im*CoupledDipoles.Γ * √(0.5s)* CoupledDipoles._scalar_laser_field(laser, sensor)
        @test laser_field(laser, sensor)[1][1] ≈ laser_expected
        @test laser_field(simulation, sensor)[1][1] ≈ laser_expected

        sensors = [ 9  3  2   3  10;
                    7  5  4  10   4;
                    2  8  1   7   4]
        laser_expected = map(eachcol(sensors)) do oneSensor
            -0.5im * CoupledDipoles.Γ * √(0.5s)* CoupledDipoles._scalar_laser_field(laser, oneSensor)
        end
        @test all(laser_field(laser, sensors) .≈ laser_expected)
        @test all(laser_field(laser, sensors)[1] .≈ laser_field(simulation, sensors)[1])

    end
    @testset "Vectorial Scattering - Single Atom" begin
        using LinearAlgebra
        atoms = Atom(Cube(), Matrix([1.0 1.0 1.0]'), 1.0)
        sensor = Matrix([-1000 -1000 -500]')
        β = [3, 4im, 5.0]
        E_μ = CoupledDipoles._vectorial_scattering_far_field(atoms, β, sensor)

        R = norm(sensor)
        n = sensor./R
        nx, ny, nz = n[1], n[2], n[3]
        C = cis(-nx - ny - nz)
        E_x_expected = +(im/2)*(3/2)*(cis(R)/R)*C*((1 - nx^2)*3 + -nx*ny*4im - nx*nz*5)
        E_y_expected = +(im/2)*(3/2)*(cis(R)/R)*C*( -ny*nx*3 + (1 -ny^2)*4im - ny*nz*5)
        E_z_expected = +(im/2)*(3/2)*(cis(R)/R)*C*( -nz*nx*3 - nz*ny*4im + (1 - nz^2)*5)
        @test all(E_μ .≈ [E_x_expected, E_y_expected, E_z_expected])




        sensor = Matrix([100 -1000 0]')
        β = [1, -2im, 15.0]
        E_μ = CoupledDipoles._vectorial_scattering_far_field(atoms, β, sensor)

        R = norm(sensor)
        n = sensor./R
        nx, ny, nz = n[1], n[2], n[3]
        C = cis(-nx - ny)
        E_x_expected = +(im/2)*(3/2)*(cis(R)/R)*C*((1 - nx^2) + nx*ny*2im - nx*nz*15)
        E_y_expected = +(im/2)*(3/2)*(cis(R)/R)*C*( -ny*nx - (1 -ny^2)*2im - ny*nz*15)
        E_z_expected = +(im/2)*(3/2)*(cis(R)/R)*C*( -nz*nx + nz*ny*2im + (1 - nz^2)*15)

        @test all(E_μ .≈ [E_x_expected, E_y_expected, E_z_expected])
    end

    @testset "angles directions" begin
        @test all(rad2deg.(laser_angles([1,1,0])) .≈ (90, 45))
        @test all(rad2deg.(laser_angles([0,0,+1])) .≈ (0, 0))
        @test all(rad2deg.(laser_angles([1,0,1])) .≈ (45, 0))
    end

    function _test_sum_all_off_diag(beta, r)
        result = zero(eltype(beta))
        for j in eachindex(beta), m in eachindex(beta)
            if j ≠ m
                rjm = sqrt((r[j, 1] - r[m, 1])^2 + (r[j, 2] - r[m, 2])^2 + (r[j, 3] - r[m, 3])^2)
                result += (sin(rjm) / rjm) * beta[j] * conj(beta[m])
            end
        end
        return real(result)
    end
    @testset "Scattered Power" begin
        w₀, s, Δ = 4π, 1e-5, 1.0

        N = 1
        atoms = Atom(Cylinder(), cylinder_inputs(N, 0.3)...)
        laser = Laser(Gaussian3D(w₀), s, Δ)
        simulation = NonLinearOptics(MeanField(), atoms, laser)

        single_atom_state = steady_state(simulation)

        P_tot = scattered_power(simulation, single_atom_state; part=:total)
        P_def = scattered_power(simulation, single_atom_state)
        @test P_tot ≈ P_def


        P_coh = scattered_power(simulation, single_atom_state; part=:coherent)
        P_inc = scattered_power(simulation, single_atom_state; part=:incoherent)

        σ⁻ = single_atom_state[1]
        R = CoupledDipoles.how_far_is_farField(simulation)
        expected_coh = 4π*(1^2)/(2*R)*abs2(σ⁻)
        @test expected_coh ≈ P_coh

        σᶻ = real(single_atom_state[2])
        expected_inc = 4π*(1^2)/(2*R)*(-abs2(σ⁻) + 0.5*(1 + σᶻ))
        @test expected_inc ≈ P_inc



        N = 200
        atoms = Atom(Cylinder(), cylinder_inputs(N, 0.3)...)
        laser = Laser(Gaussian3D(w₀), s, Δ)
        simulation = NonLinearOptics(MeanField(), atoms, laser)

        atomic_state = steady_state(simulation)

        result_1 = CoupledDipoles.sum_upper_diagonal(atomic_state[1:N], transpose(atoms.r));
        result_2 = _test_sum_all_off_diag(atomic_state[1:N], transpose(atoms.r));

        @test result_1 ≈ result_2


        N = 10
        atoms = Atom(Cylinder(), cylinder_inputs(N, 0.3)...)
        laser = Laser(Gaussian3D(w₀), s, Δ)

        simulation_scalar = LinearOptics(Scalar(), atoms, laser)
        ss_scalar = steady_state(simulation_scalar)

        simulation_meanfield = NonLinearOptics(MeanField(), atoms, laser)
        # !!! TRICKING the solutions to be equal !!!
        ss_meanfield = vcat(vec(ss_scalar[1:N]), rand(N))

        P_scalar = scattered_power(simulation_scalar, ss_scalar; part=:coherent)
        P_meanfield = scattered_power(simulation_meanfield, ss_meanfield; part=:coherent)

        @test P_scalar ≈ P_meanfield
    end
end