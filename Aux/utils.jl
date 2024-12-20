#########################################
###       Auxiliary Functions       #####
######################################### 

function Tstudent(x::Real, η::Real)
    a = gamma((1 + η) / (2 * η))/(sqrt(π * (1 - 2 * η)/η) * gamma(1 / (2 * η)));
    b = (1 + η * x^2/(1 - 2 * η))^((η + 1)/(2 * η));
    return a/b;
end

function probability_regime_given_time_n(p::Real, q::Real, σ::Vector{Float64}, r::Real, Pt::Real)
    numA = (1 - q) * pdf(Normal(0, σ[2]), r) * (1 - Pt);
    numB = p * pdf(Normal(0, σ[1]), r) * Pt;
    deno = pdf(Normal(0, σ[1]), r) * Pt + pdf(Normal(0, σ[2]), r) * (1 - Pt);
    l = numA/deno + numB/deno;
    return l;
end

function probability_regime_given_time_it(p::Real, q::Real, σ::Vector{Float64}, r::Real, Pt::Real, η::Real)
    numA = (1 - q) * 1/σ[2] * Tstudent(r / σ[2], η) *(1 - Pt);
    numB = p * 1/σ[1] * Tstudent(r / σ[1], η) * Pt;
    deno = 1/σ[1] * Tstudent(r / σ[1], η) * Pt +  1/σ[2] * Tstudent(r / σ[2], η) *(1 - Pt);
    return numA/deno + numB/deno;
end

function probability_regime_given_time_t(p::Real, q::Real, σ::Vector{Float64}, r::Real, Pt::Real, ν::Real)
    numA = (1 - q) * sqrt(ν/(ν-2))/σ[2] * pdf(TDist(ν), r*sqrt(ν/(ν-2))/σ[2]) * (1 - Pt);
    numB = p * sqrt(ν/(ν-2))/σ[1] * pdf(TDist(ν), r*sqrt(ν/(ν-2))/σ[1]) * Pt;
    deno = sqrt(ν/(ν-2))/σ[1] * pdf(TDist(ν), r*sqrt(ν/(ν-2))/σ[1]) * Pt + sqrt(ν/(ν-2))/σ[2] * pdf(TDist(ν), r*sqrt(ν/(ν-2))/σ[2]) * (1 - Pt);
    return numA/deno + numB/deno;
end

function param_transform(t_param)
    param = similar(t_param);
    param[1:2] .= 10 ./ (1 .+ exp.(-t_param[1:2]));
    param[3:4] .= exp.(-t_param[3:4]) ./ (1 .+ exp.(-t_param[3:4]) .+ exp.(-t_param[5:6]));
    param[5:6] .= exp.(-t_param[5:6]) ./ (1 .+ exp.(-t_param[3:4]) .+ exp.(-t_param[5:6]));
    param[7:8] .= 1 ./(1 .+ exp.(-t_param[7:8])); 
    return param;
end

function min_var_n(α::Float64, p1::Float64, p2::Float64, σ₁::Float64, σ₂::Float64, x)
    return (α - p1 * cdf(Normal(0, σ₁), -x[1]) - p2 * cdf(Normal(0, σ₂), -x[1]))^2;
end

function min_var_t(α::Float64, p1::Float64, p2::Float64, σ₁::Float64, σ₂::Float64, ν::Float64, x)
    return (α - p1 * cdf(TDist(ν), -x[1]*sqrt(ν/(ν-2))/σ₁) - p2 * cdf(TDist(ν), -x[1]*sqrt(ν/(ν-2))/σ₂))^2;
end
