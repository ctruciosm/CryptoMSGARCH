#################################################
###  RSGARCH GPD: Simulate RSGARCH Models    ####
#################################################


#################################################
function simulate_gray(n, distri, ω, α, β, P, burnin)
    # P = [p11, 1 - p22; 1 - p11, p22]
    ntot = n + burnin + 1;
    k = length(ω);
    h = Matrix{Float64}(undef, ntot, k + 1);
    r = Vector{Float64}(undef, ntot);
    s = Vector{Int32}(undef, ntot);
    Pt = Vector{Float64}(undef, ntot);

    e = ifelse(distri == "norm", rand(Normal(), ntot), sqrt(5/7).* rand(TDist(7), ntot));
    h[1, 1:k] .= ω ./ (1 .- (α .+ β));
    p = P[1, 1];
    q = P[2, 2];
    Pt[1] = (1 - q) / (2 - p - q);  
    h[1, k + 1] = Pt[1] * h[1, 1] + (1 - Pt[1]) * h[1, 2];
    s[1] = wsample([1, 2], [Pt[1], 1 - Pt[1]])[1];
    r[1] = e[1] * sqrt(h[1, s[1]]);
    if distri == "norm"
        for i = 2:ntot
            h[i, 1:k] = ω .+ α.* r[i - 1]^2 + β.* h[i - 1, k + 1];
            Pt[i] = probability_regime_given_time_n(p, q, sqrt.(h[i - 1, :]), r[i - 1], Pt[i - 1]);
            h[i, k + 1] = Pt[i] * h[i, 1] + (1 - Pt[i]) * h[i, 2];
            s[i] = wsample([1, 2], P[:, s[i-1]])[1]; 
            r[i] = e[i] * sqrt(h[i, s[i]]);
        end
    else
        for i = 2:ntot
            h[i, 1:k] = ω .+ α.* r[i - 1]^2 + β.* h[i - 1, k + 1];
            Pt[i] = probability_regime_given_time_t(p, q, sqrt.(h[i - 1, :]), r[i- 1], Pt[i - 1], 7);
            h[i, k + 1] = Pt[i] * h[i, 1] + (1 - Pt[i]) * h[i, 2];
            s[i] = wsample([1, 2], P[:, s[i-1]])[1];
            r[i] = e[i] * sqrt(h[i, s[i]]);
        end
    end
    return r[burnin + 1: end - 1], h[burnin + 1: end, :], Pt[burnin + 1: end], s[burnin + 1: end];
end
#################################################
function simulate_haas(n, distri, ω, α, β, P, burnin)
    # P = [p11, 1 - p22; 1 - p11, p22]
    ntot = n + burnin + 1;
    s = Vector{Int32}(undef, ntot);
    k = length(ω);
    h = Matrix{Float64}(undef, ntot, k + 1);
    r = Vector{Float64}(undef, ntot);
    Pt = Vector{Float64}(undef, ntot);
    M = Matrix{Float64}(undef, 4, 4);
    I4 = [1.0 0.0 0.0 0.0; 0.0 1.0 0.0 0.0; 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 1.0];

    e = ifelse(distri == "norm", rand(Normal(), ntot), sqrt(5/7).* rand(TDist(7), ntot));
    p = P[1, 1];
    q = P[2, 2];
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
    s[1] =  1;                                  
    r[1] = e[1] * sqrt(h[1, s[1]]);
    if distri == "norm"
        for i = 2:ntot
            h[i, 1:k] .= ω .+ α.* r[i - 1]^2 + β.* h[i - 1, 1:k];
            Pt[i] = probability_regime_given_time_n(p, q, sqrt.(h[i - 1, 1:k]), r[i - 1], Pt[i - 1]);
            h[i, k + 1] = Pt[i] * h[i, 1] + (1 - Pt[i]) * h[i, 2];
            s[i] = wsample([1, 2], P[:, s[i-1]])[1];
            r[i] = e[i] * sqrt(h[i, s[i]]);
        end
    else
        for i = 2:ntot
            h[i, 1:k] = ω .+ α.* r[i - 1]^2 + β.* h[i - 1, 1:k];
            Pt[i] = probability_regime_given_time_t(p, q, sqrt.(h[i - 1, 1:k]), r[i- 1], Pt[i - 1], 7);
            h[i, k + 1] = Pt[i] * h[i, 1] + (1 - Pt[i]) * h[i, 2];
            s[i] = wsample([1, 2], P[:, s[i-1]])[1];
            r[i] = e[i] * sqrt(h[i, s[i]]);
        end
    end
    return r[burnin + 1: end - 1], h[burnin + 1: end, :], Pt[burnin + 1: end], s[burnin + 1: end];
end
#################################################
function simulate_klaassen(n, distri, ω, α, β, P, burnin)
    # P = [p11, 1 - p22; 1 - p11, p22]
    ntot = n + burnin + 1;
    k = length(ω);
    h = Matrix{Float64}(undef, ntot, k + 1);
    r = Vector{Float64}(undef, ntot);
    s = Vector{Int32}(undef, ntot);
    Pt = Vector{Float64}(undef, ntot);
    A = Matrix{Float64}(undef, 2, 2);
    I2 = [1.0 0.0; 0.0 1.0];

    e = ifelse(distri == "norm", rand(Normal(), ntot), sqrt(5/7).* rand(TDist(7), ntot));
    h[1, 1:k] .= ω ./ (1 .- (α .+ β));
    p = P[1, 1];
    q = P[2, 2];
    Pt[1] = (1 - q) / (2 - p - q);  
    A[1, 1] = p * (α[1] + β[1]);
    A[1, 2] = (1 - p) * (α[1] + β[1]);
    A[2, 1] = (1 - q) * (α[2] + β[2]);
    A[2, 2] = q * (α[2] + β[2]);
    h[1, 1:k] .= inv(I2 - A) * ω;                    
    h[1, k + 1] = Pt[1] * h[1, 1] + (1 - Pt[1]) * h[1, 2];
    s[1] = wsample([1, 2], [Pt[1], 1 - Pt[1]])[1];
    r[1] = e[1] * sqrt(h[1, s[1]]);
    if distri == "norm"
        for i = 2:ntot
            Pt[i] = probability_regime_given_time_n(p, q, sqrt.(h[i - 1, :]), r[i - 1], Pt[i - 1]);
            h[i, 1] = ω[1] + α[1] * r[i - 1]^2 + β[1] * (P[1,1] * pdf(Normal(0, sqrt(h[i - 1, 1])), r[i - 1]) * Pt[i - 1] * h[i - 1, 1] + P[2,1] * pdf(Normal(0, sqrt(h[i - 1, 2])), r[i - 1]) * (1 - Pt[i - 1]) * h[i - 1, 2])/(Pt[i] * (pdf(Normal(0, sqrt(h[i - 1, 1])), r[i - 1]) * Pt[i - 1] + pdf(Normal(0, sqrt(h[i - 1, 2])), r[i - 1]) * (1 - Pt[i - 1])));
            h[i, 2] = ω[2] + α[2] * r[i - 1]^2 + β[2] * (P[1,2] * pdf(Normal(0, sqrt(h[i - 1, 1])), r[i - 1]) * Pt[i - 1] * h[i - 1, 1] + P[2,2] * pdf(Normal(0, sqrt(h[i - 1, 2])), r[i - 1]) * (1 - Pt[i - 1]) * h[i - 1, 2])/((1 - Pt[i]) * (pdf(Normal(0, sqrt(h[i - 1, 1])), r[i - 1]) * Pt[i - 1] + pdf(Normal(0, sqrt(h[i - 1, 2])), r[i - 1]) * (1 - Pt[i - 1])));
            h[i, k + 1] = Pt[i] * h[i, 1] + (1 - Pt[i]) * h[i, 2];
            s[i] = wsample([1, 2], P[:, s[i-1]])[1]; 
            r[i] = e[i] * sqrt(h[i, s[i]]);
        end
    else
        for i = 2:ntot
            Pt[i] = probability_regime_given_time_t(p, q, sqrt.(h[i - 1, :]), r[i- 1], Pt[i - 1], 7);
            h[i, 1] = ω[1] + α[1] * r[i - 1]^2 + β[1] * (P[1,1] * (sqrt(7/5)/sqrt(h[i - 1, 1]) * pdf(TDist(7), r[i - 1]*sqrt(7/5)/sqrt(h[i - 1, 1]))) * Pt[i - 1] * h[i - 1, 1] + P[2,1] * (sqrt(7/5)/sqrt(h[i - 1, 2]) * pdf(TDist(7), r[i - 1]*sqrt(7/5)/sqrt(h[i - 1, 2]))) * (1 - Pt[i - 1]) * h[i - 1, 2])/(Pt[i] * ((sqrt(7/5)/sqrt(h[i - 1, 1]) * pdf(TDist(7), r[i - 1]*sqrt(7/5)/sqrt(h[i - 1, 1]))) * Pt[i - 1] + (sqrt(7/5)/sqrt(h[i - 1, 2]) * pdf(TDist(7), r[i - 1]*sqrt(7/5)/sqrt(h[i - 1, 2]))) * (1 - Pt[i - 1])));
            h[i, 2] = ω[2] + α[2] * r[i - 1]^2 + β[2] * (P[1,2] * (sqrt(7/5)/sqrt(h[i - 1, 1]) * pdf(TDist(7), r[i - 1]*sqrt(7/5)/sqrt(h[i - 1, 1]))) * Pt[i - 1] * h[i - 1, 1] + P[2,2] * (sqrt(7/5)/sqrt(h[i - 1, 2]) * pdf(TDist(7), r[i - 1]*sqrt(7/5)/sqrt(h[i - 1, 2]))) * (1 - Pt[i - 1]) * h[i - 1, 2])/((1 - Pt[i]) * ((sqrt(7/5)/sqrt(h[i - 1, 1]) * pdf(TDist(7), r[i - 1]*sqrt(7/5)/sqrt(h[i - 1, 1]))) * Pt[i - 1] + (sqrt(7/5)/sqrt(h[i - 1, 2]) * pdf(TDist(7), r[i - 1]*sqrt(7/5)/sqrt(h[i - 1, 2]))) * (1 - Pt[i - 1])));
            h[i, k + 1] = Pt[i] * h[i, 1] + (1 - Pt[i]) * h[i, 2];
            s[i] = wsample([1, 2], P[:, s[i-1]])[1];
            r[i] = e[i] * sqrt(h[i, s[i]]);
        end
    end
    return r[burnin + 1: end - 1], h[burnin + 1: end, :], Pt[burnin + 1: end], s[burnin + 1: end];
end


