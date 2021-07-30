"""
For localized modes, the L₁-norm is more reliable, but, if the result is extended,
the fit procedure is unstable, and it more meaninfull to apply fit with L₂-norm instead.
"""
function get_single_ξ_and_R1(x, y)
    x_use, y_use = select_points(x, y)
    A_L₁, ξ_L₁, R1, y_fit_L₁ = get_coeffs_L1(x_use, y_use)

    if is_R1_localized(R1)
        return ξ_L₁, R1
    else
        A_L₂, ξ_L₂, R2, y_fit_L₂ = get_coeffs_L2(x_use, y_use)
        return ξ_L₂, R2
    end
end
function get_coeffs_L1(x_use, y_use)
    A_L₁, ξ_L₁ = get_exp_fit_coeffs_with_L1(x_use, y_use)
    y_fit_L₁ = model_exp(x_use, [A_L₁, ξ_L₁])
    R1 = compute_R1(y_use, y_fit_L₁)
    return A_L₁, ξ_L₁, R1, y_fit_L₁
end
function get_coeffs_L2(x_use, y_use)
    A_L₂, ξ_L₂ = get_exp_fit_coeffs_with_L2(x_use, y_use)
    y_fit_L₂ = model_exp(x_use, [A_L₂, ξ_L₂])
    R2 = compute_R2(y_use, y_fit_L₂)
    return A_L₂, ξ_L₂, R2, y_fit_L₂
end
function is_R1_localized(R1)
    # is an arbitrary decision to say that a mode is  "localized" if R1 > "somthing"
    if R1 > R1_threshold
        return true
    else
        return false
    end
end

"""
    x_use, y_use = select_points(DCM, ψ²; probability_threshold = 0.999)

computes how many points since the center of mass are important based upon the `probability_threshold`
and return them in variables `x_use`, `y_use`
"""
function select_points(DCM, ψ²)
    if length(ψ²) < 10
        error("Very few data. Need at least 10 points")
    end

    cumulative_probability = cumsum(ψ²)
    threshold = sum(cumulative_probability .< probability_threshold)

    x_use = DCM[1:threshold]
    y_use = ψ²[1:threshold]

    return x_use, y_use
end

######### Now, we actually compute the coefficients
"""
    A,ξ = get_exp_fit_coeffs_with_L1(x, y)

Assumes that (x,y) have a exponential decay of form `y(x) = A*exp(-x/ξ)`
"""
function get_exp_fit_coeffs_with_L1(x, y)
    fit_result = optimize(xx -> model_exp_minimize(xx, x, y), rand(2), Options(iterations = 30))
    A = minimizer(fit_result)[1]
    ξ = minimizer(fit_result)[2]
    return A, ξ
end
model_exp_minimize(x0, x, y) = sum(abs.(x0[1] * exp.(-x ./ x0[2]) .- y))
@. model_exp(x, p) = p[1] * exp.(-x ./ p[2])

"""
    A,ξ = get_exp_fit_coeffs_with_L2(x, y)

Assumes that (x,y) have a exponential decay in linear form, that is, `log(y(x)) = log(A) - x/ξ`
"""
function get_exp_fit_coeffs_with_L2(x, y)
    fit = curve_fit(model_linear, x, log.(y), [maximum(log.(y)), rand()])
    A = exp(coef(fit)[1])
    ξ = coef(fit)[2]
    return A, ξ
end
@. model_linear(x, p) = p[1] - x / p[2]

######### Quality of the Fit

"""
    compute_R2(y, y_predict)

returns the R2 (adjusted) for a data set "y", and its modeled values "y_predict"

both inputs should have the same size

example:
actual = [2, 4, 6, 7] 
predict = [2.601,3.83,5.059,7.517]

compute_R2(actual, predict)
0.8430978644067797

Ref: 
Coefficient of Determination, R-squared, ASK Academic Skills Kit, Newcastle University
Coefficient of Determination (R-Squared), Matlab Documentation
"""
compute_R2(y, y_predict) = generic_R(log.(y), log.(y_predict), SSE_R2, SST_R2)
compute_R1(y, y_predict) = generic_R(y, y_predict, SSE_R1, SST_R1)

function generic_R(y, y_predict, SSE_fnc::Function, SST_fnc::Function)
    if size(y) ≠ size(y_predict)
        @error("Inputs does not have similar dimension")
    end

    average_y = sum(y) / length(y)

    SSE = SSE_fnc(y, y_predict)# the sum of squared error
    SST = SST_fnc(y, average_y) # sum of squared total

    n = length(y)
    p = 2 # number of regression coefficients, for linear, p=2
    R2_adjusted = 1 - ((n - 1) / (n - p)) * (SSE / SST)

    if R2_adjusted < 0
        R2_adjusted = 1e-16
    end
    return R2_adjusted
end

SSE_R1(y, y_predict) = sum(abs(y[i] - y_predict[i]) for i in eachindex(y)) # the sum of squared error
SST_R1(y, average_y) = sum(abs(y[i] - average_y) for i in eachindex(y)) # sum of squared total

SSE_R2(y, y_predict) = sum((y[i] - y_predict[i])^2 for i in eachindex(y)) # the sum of squared error
SST_R2(y, average_y) = sum((y[i] - average_y)^2 for i in eachindex(y)) # sum of squared total
