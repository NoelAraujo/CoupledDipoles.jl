using CoupledDipoles
using Test

@testset begin
    @testset "Green Matrix 1 Atom" begin
        singleAtom = Matrix([0 -1 0]')
        atoms = Atom(Cube(), singleAtom, 1)

        s, Δ = 1e-6, 1.5
        laser = Laser(PlaneWave3D([0,0,1]), s, Δ)

        ## Scalar case
        problem = LinearOptics(Scalar(), atoms, laser)
        G = interaction_matrix(problem)
        G_expected = im*1.5 - 1/2
        @test all(G .≈ G_expected)

        ## Vectorial case
        problem = LinearOptics(Vectorial(), atoms, laser)
        G = interaction_matrix(problem)
        G_expected = fill(im*1.5 - 1/3, 3, 3)
        @test all(G .≈ G_expected)
    end


    @testset "Green Matrix 2 Atoms" begin
        r1 = [0,-1,0]
        r2 = [3,1,0]
        r = hcat(r1,r2)
        atoms = Atom(Cube(), Array(r), 1)
        Γ = CoupledDipoles.Γ
        s, Δ = 1e-6, 1.0
        laser = Laser(PlaneWave3D([0,0,1]), s, Δ)

        ## Scalar case
        problem = LinearOptics(Scalar(), atoms, laser)
        G = interaction_matrix(problem)

        G_expected = zeros(ComplexF64, 2, 2)
        G_expected[1,1] = im*1.0 - 1/2
        G_expected[1,2] = -(Γ / 2) *(cis(√13))/(im√13)
        G_expected[2,1] = -(Γ / 2) *(cis(√13))/(im√13)
        G_expected[2,2] = im*1.0 - 1/2

        @test all(G .≈ G_expected)

        ## Vectorial case
        problem = LinearOptics(Vectorial(), atoms, laser)
        G = interaction_matrix(problem)

        G_expected = zeros(ComplexF64, 6, 6)
        # r_11
        G_expected[1,1] = im*1.0 - 1/3
        G_expected[1,3] = im*1.0 - 1/3
        G_expected[1,5] = im*1.0 - 1/3
        G_expected[3,1] = im*1.0 - 1/3
        G_expected[3,3] = im*1.0 - 1/3
        G_expected[3,5] = im*1.0 - 1/3
        G_expected[5,1] = im*1.0 - 1/3
        G_expected[5,3] = im*1.0 - 1/3
        G_expected[5,5] = im*1.0 - 1/3

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
        G_expected[2,2] = im*1.0 - 1/3
        G_expected[2,4] = im*1.0 - 1/3
        G_expected[2,6] = im*1.0 - 1/3
        G_expected[4,2] = im*1.0 - 1/3
        G_expected[4,4] = im*1.0 - 1/3
        G_expected[4,6] = im*1.0 - 1/3
        G_expected[6,2] = im*1.0 - 1/3
        G_expected[6,4] = im*1.0 - 1/3
        G_expected[6,6] = im*1.0 - 1/3

        @test all(G .≈ G_expected)
    end


    @testset "3D scalar: 2 atoms" begin
        r =[1 2 0;
            1 0 1.0]
        atoms = Atom(Cube(), Array(transpose(r)), 10)
        laser = Laser(PlaneWave3D([0,0,1]), 1e-6, 1.0)
        problem = LinearOptics(Scalar(), atoms, laser)
        ωₙ, Γₙ = get_spectrum(problem)

        Γₙ_ans = [0.675922453937850, 0.324077546062151]
        ωₙ_ans = [0.861973588757495, 1.13802641124251]
        @test all(sort(ωₙ_ans) .≈ sort(ωₙ))
        @test all(sort(Γₙ_ans) .≈ sort(Γₙ))

        βₛₛ = steady_state(problem)
        βₛₛ_ans = [0.000695478188875652-0.000416450295287958*im,
                    0.000554066244758147+0.000208376732488666*im ]
        @test all(βₛₛ_ans .≈ βₛₛ)
    end
    @testset "3D scalar: 3 atoms" begin
        r =[ 0 0 0;
                1 1 1;
            -1 1 0.5]
        atoms = Atom(Cube(), Array(transpose(r)), 10)
        laser = Laser(PlaneWave3D([0,0,1]), 1e-6, 0.0)
        problem = LinearOptics(Scalar(), atoms, laser)
        ωₙ, Γₙ = get_spectrum(problem)

        Γₙ_ans = [0.293248907595325, 0.15447403278592, 1.05227705967608]
        ωₙ_ans = [-0.0853652197846037, -0.0371827945713394,0.12254801435943]
        @test all(sort(ωₙ_ans) .≈ sort(ωₙ))
        @test all(sort(Γₙ_ans) .≈ sort(Γₙ))

        R_nm = get_pairwise_matrix(atoms.r)
        R_nm_ans = [0 1.7320508075688 1.5;
                    1.73205080756888 0 2.06155281280883;
                    1.5 2.06155281280883 0]
        @test all(R_nm_ans .≈ R_nm)

        G = interaction_matrix(problem)
        G_ans = [-0.5 -0.284930049591257-0.04634868038121*im -0.332498328868018+0.0235790672225676*im;
                    -0.284930049591257-0.046348680381261*im -0.5 -0.213910742134572-0.11430539712927*im;
                    -0.332498328868018+0.0235790672225676*im -0.213910742134572-0.11430539712197*im -0.5]
        @test all(G .≈ G_ans)
    end

end