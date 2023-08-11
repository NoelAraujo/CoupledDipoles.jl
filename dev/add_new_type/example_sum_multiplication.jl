using Random, LinearAlgebra, Tullio

function truncation(A, B, C, AB, AC, BC, N)
	f = (-2.0 .* reshape(A * transpose(B), (N, N, 1)) .* reshape(C, (1, 1, N))
		 + reshape(BC, (1, N, N)) .* reshape(A, (N, 1, 1))
		 + reshape(AB, (N, N, 1)) .* reshape(C, (1, 1, N))
		 + reshape(AC, (N, 1, N)) .* reshape(B, (1, N, 1)))
	return f
end

function getExclusionMatrices(N)
	# creating matrix of "Bool" uses less memory than "Integer"
	exclusionDiagonal = ones(Bool, N, N)
	exclusionDiagonal[diagind(exclusionDiagonal)] .= false

	exclusion3D = ones(Bool, N, N)
	exclusion3D[diagind(exclusion3D)] .= false
	exclusion3D = (exclusion3D .* reshape(exclusion3D, (N, 1, N))) .* reshape(exclusion3D, (1, N, N))

	return exclusionDiagonal, exclusion3D
end
function check_matrices(a, b)
    if all(a .≈ b)
        println("all elements are equal")
    else
        @error "some elements are different "
    end
end

N = 50

Random.seed!(4555)
G = rand(N, N)

Random.seed!(52)
A, B, C, AB, AC, BC = rand(N), rand(N), rand(N), rand(N, N), rand(N, N), rand(N, N)
σσσ = truncation(A, B, C, AB, AC, BC, N)

# ----------- no constraints -----------
@time matrix_1 = let
	sum(reshape(G, (N, 1, N)) .* σσσ, dims = 3)[:, :]
end;
@time matrix_2 = @tullio s[j,m] := begin
    G[j, k]*σσσ[j, m, k]
end;
check_matrices(matrix_1, matrix_2)

# ----------- with constraints -----------
exclusionDiagonal, exclusion3D = getExclusionMatrices(N)

@time matrix_3 = let
	sum(reshape(G, (N, 1, N)) .* (exclusion3D.*σσσ), dims = 3)[:, :]
end;
@time matrix_4 = @tullio s[j,m] := begin
    if ((k ≠ j) && (k ≠ m))
        G[j, k]*σσσ[j, m, k]
    else
        0.0
    end
end;
check_matrices(matrix_3, matrix_4)

## the diagonals are different
matrix_3 .≈ matrix_4

# --> matrix_3 has ZERO at diagonal <--




# ---------------------------------------------------
# ------------------ testing tullio -----------------
# ---------------------------------------------------


σσσ_v1 = @tullio σ[j,m,k] := begin
    A[j]*B[m]*C[k]
end
t_1 = reshape(A * transpose(B), (N, N, 1)) .* reshape(C, (1, 1, N))
all(σσσ_v1 .≈ t_1)

σσσ_v2 = @tullio σ[j,m,k] := begin
    A[j]*BC[m,k]
end
t_2 = reshape(BC, (1, N, N)) .* reshape(A, (N, 1, 1))
all(σσσ_v2 .≈ t_2)

σσσ_v3 = @tullio σ[j,m,k] := begin
    AB[j,m]*C[k]
end
t_3 = reshape(AB, (N, N, 1)) .* reshape(C, (1, 1, N))
all(σσσ_v3 .≈ t_3)

σσσ_v4 = @tullio σ[j,m,k] := begin
    AC[j,k]*B[m]
end
t_4 = reshape(AC, (N, 1, N)) .* reshape(B, (1, N, 1))
all(σσσ_v4 .≈ t_4)



function truncation_v2(A, B, C, AB, AC, BC, N)
    f =  @tullio σ[j,m,k] := begin
        -2A[j]*B[m]*C[k] + A[j]*BC[m,k] + AB[j,m]*C[k] + AC[j,k]*B[m]
    end
	return f
end

all(truncation(A, B, C, AB, AC, BC, N) .≈ truncation_v2(A, B, C, AB, AC, BC, N))

@time truncation(A, B, C, AB, AC, BC, N);
@time truncation_v2(A, B, C, AB, AC, BC, N);