##################################################
###            RSGARCH: Forecasts             ####
##################################################
function fore_gray(r::Vector{Float64}, k::Int64, par, distri::String)
    # Regime 1: Low Vol
    # Regime 2: High Vol
    n = length(r);
    h = Matrix{Float64}(undef, n + 1, k + 1);
    s = Vector{Int32}(undef, n + 1);
    Pt = Vector{Float64}(undef, n + 1);

    ω = par[1:k];
    α = par[k + 1 : 2 * k];
    β = par[2 * k + 1 : 3 * k];
    p = par[3 * k + 1];
    q = par[4 * k];
    P = [p 1-q; 1-p q];

    Pt[1] = (1 - q) / (2 - p - q);              
    h[1, 1:k] .= var(r);                        
    h[1, k + 1] = Pt[1] * h[1, 1] + (1 - Pt[1]) * h[1, 2];
    s[1] = wsample([1, 2], [Pt[1], 1 - Pt[1]])[1];
    if distri == "norm"
        @inbounds for i = 2:n+1
            h[i, 1:k] = ω .+ α .* r[i - 1]^2 + β .* h[i - 1, k + 1];
            Pt[i] = probability_regime_given_time_n(p, q, sqrt.(h[i - 1, :]), r[i - 1], Pt[i - 1]);
            h[i, k + 1] = Pt[i] * h[i, 1] + (1 - Pt[i]) * h[i, 2];
            s[i] = wsample([1, 2], P[:, s[i-1]])[1]; 
        end
    else
        ν = par[4 * k + 1];
        @inbounds for i = 2:n+1
            h[i, 1:k] = ω .+ α .* r[i - 1]^2 + β .* h[i - 1, k + 1];
            Pt[i] = probability_regime_given_time_t(p, q, sqrt.(h[i- 1, :]), r[i - 1], Pt[i - 1], ν);
            h[i, k + 1] = Pt[i] * h[i, 1] + (1 - Pt[i]) * h[i, 2];     
            s[i] = wsample([1, 2], P[:, s[i-1]])[1]; 
        end
    end
    return h[end,:], Pt, s;
end
##################################################
function fore_haas(r::Vector{Float64}, k::Int64, par, distri::String)
    # Regime 1: Low Vol
    # Regime 2: High Vol
    n = length(r);
    h = Matrix{Float64}(undef, n + 1, k + 1);
    s = Vector{Int32}(undef, n + 1);
    Pt = Vector{Float64}(undef, n + 1);
    M = Matrix{Float64}(undef, 4, 4);
    I4 = [1.0 0.0 0.0 0.0; 0.0 1.0 0.0 0.0; 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 1.0];

    ω = par[1:k];
    α = par[k + 1 : 2 * k];
    β = par[2 * k + 1 : 3 * k];
    p = par[3 * k + 1];
    q = par[4 * k];
    P = [p 1-q; 1-p q];
    M[1, 1] = P[1, 1] * (α[1] + β[1]);
    M[1, 2] = 0.0;
    M[1, 3] = P[1, 2] * (α[1] + β[1]);
    M[1, 4] = 0.0;
    M[2, 1] = P[1, 1] * α[2];
    M[2, 2] = P[1, 1] * β[2];
    M[2, 3] = P[1, 2] * α[2];
    M[2, 4] = P[1, 2] * β[2];
    M[3, 1] = P[2, 1] * β[1];
    M[3, 2] = P[2, 1] * α[1];
    M[3, 3] = P[2, 2] * β[1];
    M[3, 4] = P[2, 2] * α[1];
    M[4, 1] = 0.0;
    M[4, 2] = P[2, 1] * (α[2] + β[2]);
    M[4, 3] = 0.0;
    M[4, 4] = P[2, 2] * (α[2] + β[2]);

    Pt[1] = (1 - q) / (2 - p - q);       
    π∞ = [Pt[1]; 1 - Pt[1]];      
    h[1, 1:k] .= [1.0 0.0 1.0 0.0; 0.0 1.0 0.0 1.0] * inv(I4 - M) * kronecker(π∞, ω);
    h[1, k + 1] = Pt[1] * h[1, 1] + (1 - Pt[1]) * h[1, 2];
    s[1] = wsample([1, 2], [Pt[1], 1 - Pt[1]])[1];
    if distri == "norm"
        @inbounds for i = 2:n+1
            h[i, 1:k] = ω .+ α .* r[i - 1]^2 + β .* h[i - 1, 1:k];
            Pt[i] = probability_regime_given_time_n(p, q, sqrt.(h[i - 1, 1:k]), r[i - 1], Pt[i - 1]);
            h[i, k + 1] = Pt[i] * h[i, 1] + (1 - Pt[i]) * h[i, 2];
            s[i] = wsample([1, 2], P[:, s[i-1]])[1]; 
        end
    else
        ν = par[4 * k + 1];
        @inbounds for i = 2:n+1
            h[i, 1:k] = ω .+ α .* r[i - 1]^2 + β .* h[i - 1, 1:k];
            Pt[i] = probability_regime_given_time_t(p, q, sqrt.(h[i- 1, 1:k]), r[i - 1], Pt[i - 1], ν);
            h[i, k + 1] = Pt[i] * h[i, 1] + (1 - Pt[i]) * h[i, 2];     
            s[i] = wsample([1, 2], P[:, s[i-1]])[1]; 
        end
    end
    return h[end,:], Pt, s;

end
##################################################
function fore_klaassen(r::Vector{Float64}, k::Int64, par, distri::String)
    # Regime 1: Low Vol
    # Regime 2: High Vol
    n = length(r);
    h = Matrix{Float64}(undef, n + 1, k + 1);
    s = Vector{Int32}(undef, n + 1);
    Pt = Vector{Float64}(undef, n + 1);
    A = Matrix{Float64}(undef, 2, 2);
    I2 = [1.0 0.0; 0.0 1.0];

    ω = par[1:k];
    α = par[k + 1 : 2 * k];
    β = par[2 * k + 1 : 3 * k];
    p = par[3 * k + 1];
    q = par[4 * k];
    P = [p 1-q; 1-p q];
    Pt[1] = (1 - q) / (2 - p - q);              
    A[1, 1] = p * (α[1] + β[1]);
    A[1, 2] = (1 - p) * (α[1] + β[1]);
    A[2, 1] = (1 - q) * (α[2] + β[2]);
    A[2, 2] = q * (α[2] + β[2]);
    
    h[1, 1:k] .= inv(I2 - A) * ω;                                         
    h[1, k + 1] = Pt[1] * h[1, 1] + (1 - Pt[1]) * h[1, 2];
    s[1] = wsample([1, 2], [Pt[1], 1 - Pt[1]])[1];
    if distri == "norm"
        @inbounds for i = 2:n+1
            Pt[i] = probability_regime_given_time_n(p, q, sqrt.(h[i - 1, :]), r[i - 1], Pt[i - 1]);
            h[i, 1] = ω[1] + α[1] * r[i - 1]^2 + β[1] * (P[1,1] * pdf(Normal(0, sqrt(h[i - 1, 1])), r[i - 1]) * Pt[i - 1] * h[i - 1, 1] + P[2,1] * pdf(Normal(0, sqrt(h[i - 1, 2])), r[i - 1]) * (1 - Pt[i - 1]) * h[i - 1, 2])/(Pt[i] * (pdf(Normal(0, sqrt(h[i - 1, 1])), r[i - 1]) * Pt[i - 1] + pdf(Normal(0, sqrt(h[i - 1, 2])), r[i - 1]) * (1 - Pt[i - 1])));
            h[i, 2] = ω[2] + α[2] * r[i - 1]^2 + β[2] * (P[1,2] * pdf(Normal(0, sqrt(h[i - 1, 1])), r[i - 1]) * Pt[i - 1] * h[i - 1, 1] + P[2,2] * pdf(Normal(0, sqrt(h[i - 1, 2])), r[i - 1]) * (1 - Pt[i - 1]) * h[i - 1, 2])/((1 - Pt[i]) * (pdf(Normal(0, sqrt(h[i - 1, 1])), r[i - 1]) * Pt[i - 1] + pdf(Normal(0, sqrt(h[i - 1, 2])), r[i - 1]) * (1 - Pt[i - 1])));
            h[i, k + 1] = Pt[i] * h[i, 1] + (1 - Pt[i]) * h[i, 2];
            s[i] = wsample([1, 2], P[:, s[i-1]])[1]; 
        end
    else
        η = 1/par[4 * k + 1];
        @inbounds for i = 2:n+1
            Pt[i] = probability_regime_given_time_it(p, q, sqrt.(h[i - 1, :]), r[i - 1], Pt[i - 1], η);
            h[i, 1] = ω[1] + α[1] * r[i - 1]^2 + β[1] * (P[1,1] * (1/ sqrt(h[i - 1, 1]) * Tstudent(r[i - 1] / sqrt(h[i - 1, 1]), η)) * Pt[i - 1] * h[i - 1, 1] + P[2,1] * (1/ sqrt(h[i - 1, 2]) * Tstudent(r[i - 1] / sqrt(h[i - 1, 2]), η)) * (1 - Pt[i - 1]) * h[i - 1, 2])/(Pt[i] * ((1/ sqrt(h[i - 1, 1]) * Tstudent(r[i - 1] / sqrt(h[i - 1, 1]), η)) * Pt[i - 1] + (1/ sqrt(h[i - 1, 2]) * Tstudent(r[i - 1] / sqrt(h[i - 1, 2]), η)) * (1 - Pt[i - 1])));
            h[i, 2] = ω[2] + α[2] * r[i - 1]^2 + β[2] * (P[1,2] * (1/ sqrt(h[i - 1, 1]) * Tstudent(r[i - 1] / sqrt(h[i - 1, 1]), η)) * Pt[i - 1] * h[i - 1, 1] + P[2,2] * (1/ sqrt(h[i - 1, 2]) * Tstudent(r[i - 1] / sqrt(h[i - 1, 2]), η)) * (1 - Pt[i - 1]) * h[i - 1, 2])/((1 - Pt[i]) * ((1/ sqrt(h[i - 1, 1]) * Tstudent(r[i - 1] / sqrt(h[i - 1, 1]), η)) * Pt[i - 1] + (1/ sqrt(h[i - 1, 2]) * Tstudent(r[i - 1] / sqrt(h[i - 1, 2]), η)) * (1 - Pt[i - 1])));
            h[i, k + 1] = Pt[i] * h[i, 1] + (1 - Pt[i]) * h[i, 2];
            s[i] = wsample([1, 2], P[:, s[i-1]])[1]; 
        end
    end
    return h[end,:], Pt, s;
end


##################################################
### Forecast VaR and ES
##################################################
function var_es_rsgarch(α, p1, p2, σ₁, σ₂, distri, ν = nothing)
    if distri == "norm"
        opt = optimize(par -> min_var_n(α, p1, p2, σ₁, σ₂, par), [0.0], iterations = 10000).minimizer;
        f(x) = (x * p1 * pdf(Normal(0, σ₁), x) + x * p2 * pdf(Normal(0, σ₂), x))/α
        a = -Inf;
        b = -opt[1];
        ES, E = quadgk(f, a, b, rtol=1e-8);
    elseif distri == "student"
        opt = optimize(par -> min_var_t(α, p1, p2, σ₁, σ₂, ν, par), [0.0], iterations = 10000).minimizer;
        g(x) = (x * p1 * sqrt(ν/(ν-2))/σ₁ * pdf(TDist(ν), x*sqrt(ν/(ν-2))/σ₁) + x * p2 * + sqrt(ν/(ν-2))/σ₂ * pdf(TDist(ν), x*sqrt(ν/(ν-2))/σ₂))/α;
        a = -Inf;
        b = -opt[1];
        ES, E = quadgk(g, a, b, rtol=1e-8);
    else
        opt = NaN64
        ES = NaN64
        println("Error: Distribution should be Normal ('norm') or Student-T ('student').")
    end
    return [-opt; ES];
end


