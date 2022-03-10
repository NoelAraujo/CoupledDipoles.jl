using CoupledDipoles
using Test

@testset "CoupledDipoles.jl" begin
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
        βₛₛ_ans = -[0.000695478188875652-0.000416450295287958*im,
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