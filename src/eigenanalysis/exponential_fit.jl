"""
'satman2015' has better precision
'lta' is faster
"""
function get_single_ξ_and_R2(x, y; regression_method = satman2015, probability_threshold=0.999)
    x_use, y_use = select_points(x, y; probability_threshold=probability_threshold)
    A, B, y_fit = linear_fit_robust(x_use, log10.(y_use); regression_method = regression_method)

    ξ = abs(-1/B) # 'abs' to hangle extended modes
    R2 = my_r2score(log10.(y_use), y_fit)

    return ξ, R2
end
function linear_fit_robust(x, y; regression_method = satman2015)
	reg = createRegressionSetting(@formula(y ~ x), DataFrame([:x => x, :y => y]))
	result = regression_method(reg)

	A, B = result["betas"]
	y_fit = A .+ (B .* x)

	return A, B, y_fit
end

function my_r2score(y_true::AbstractVector{T}, y_pred::AbstractVector{T}) where T<:Real
    n = length(y_true)
    y_mean = mean(y_true)
    ss_total = sum((y .- y_mean)^2 for y in y_true)
    ss_residual = sum((y_true[i] - y_pred[i])^2 for i in 1:n)
    r2 = max(1 - ss_residual / ss_total, exp10(-16))
    return r2
end

"""
    x_use, y_use = select_points(DCM::Vector, ψ²::Vector)

computes how many points since the center of mass are important based upon the `probability_threshold` (constant defined as 0.9999)
"""
function select_points(DCM::Vector, ψ²::Vector; probability_threshold=0.999)
    if length(ψ²) < 10
        error("Very few data. Need at least 10 points")
    end

    cumulative_probability = cumsum(ψ²)
    threshold = max(sum(cumulative_probability .< probability_threshold), 10)

    x_use = DCM[1:threshold]
    y_use = ψ²[1:threshold]

    return x_use, y_use
end