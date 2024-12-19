##################################################
###  RSGARCH: Maximum Likelihood Estimation   ####
##################################################

##################################################
### Likelihoods: par .= [ω, α, β, p11, p22, 1/gl]
##################################################
function gray_likelihood(r::Vector{Float64}, k::Int64, distri::String, par)
    n = length(r);
    h = Matrix{Float64}(undef, n, k + 1);
    Pt = Vector{Float64}(undef, n);
    log_lik = Vector{Float64}(undef, n - 1);
    # Transformations
    ω = 10 ./ (1 .+ exp.(-par[1:2]));
    α = exp.(-par[3:4]) ./ (1 .+ exp.(-par[3:4]) .+ exp.(-par[5:6]));
    β = exp.(-par[5:6]) ./ (1 .+ exp.(-par[3:4]) .+ exp.(-par[5:6]));
    p = 1 ./(1 .+ exp(-par[7]));
    q = 1 ./(1 .+ exp(-par[8]));
    # Likelihood
    Pt[1] = (1 - q) / (2 - p - q);              
    h[1, 1:k] .= var(r);                       
    h[1, k + 1] = Pt[1] * h[1, 1] + (1 - Pt[1]) * h[1, 2];
    if distri == "norm"
        @inbounds for i = 2:n
            h[i, 1:k] .= ω .+ α .* r[i - 1]^2 + β .* h[i - 1, k + 1];
            Pt[i] = probability_regime_given_time_n(p, q, sqrt.(h[i - 1, :]), r[i - 1], Pt[i - 1]);
            h[i, k + 1] = Pt[i] * h[i, 1] + (1 - Pt[i]) * h[i, 2];
            log_lik[i - 1] = log(pdf(Normal(0, sqrt(h[i, 1])), r[i]) * Pt[i] + pdf(Normal(0, sqrt(h[i, 2])), r[i]) * (1 - Pt[i]));
        end
    else
        η = 1 / (2 + exp(-par[9]));
        @inbounds for i = 2:n
            h[i, 1:k] .= ω .+ α .* r[i - 1]^2 + β .* h[i - 1, k + 1];
            Pt[i] = probability_regime_given_time_it(p, q, sqrt.(h[i- 1, :]), r[i - 1], Pt[i - 1], η);
            h[i, k + 1] = Pt[i] * h[i, 1] + (1 - Pt[i]) * h[i, 2];
            log_lik[i - 1] = log(1/ sqrt(h[i, 1]) * Tstudent(r[i] / sqrt(h[i, 1]), η)* Pt[i]  + 1 / sqrt(h[i, 2]) * Tstudent(r[i] / sqrt(h[i, 2]), η) * (1 - Pt[i]));
        end
    end
    return -sum(log_lik)/2;
end
##################################################
function haas_likelihood(r::Vector{Float64}, k::Int64, distri::String, par)
    n = length(r);
    h = Matrix{Float64}(undef, n, k);
    Pt = Vector{Float64}(undef, n);
    log_lik = Vector{Float64}(undef, n - 1);
    M = Matrix{Float64}(undef, 4, 4);
    I4 = [1.0 0.0 0.0 0.0; 0.0 1.0 0.0 0.0; 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 1.0];
    # Transformations
    ω = 10 ./ (1 .+ exp.(-par[1:2]));
    α = exp.(-par[3:4]) ./ (1 .+ exp.(-par[3:4]) .+ exp.(-par[5:6]));
    β = exp.(-par[5:6]) ./ (1 .+ exp.(-par[3:4]) .+ exp.(-par[5:6]));
    p = 1 ./(1 .+ exp(-par[7]));
    q = 1 ./(1 .+ exp(-par[8]));
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
    # Likelihood
    Pt[1] = (1 - q) / (2 - p - q);       
    π∞ = [Pt[1]; 1 - Pt[1]];      
    h[1, :] .= [1.0 0.0 1.0 0.0; 0.0 1.0 0.0 1.0] * inv(I4 - M) * kronecker(π∞, ω);
    if (distri == "norm")
        @inbounds for i = 2:n
            h[i, :] .= ω .+ α .* r[i - 1]^2 + β .* h[i - 1, :];
            Pt[i] = probability_regime_given_time_n(p, q, sqrt.(h[i - 1, :]), r[i - 1], Pt[i - 1]);
            log_lik[i - 1] = log(pdf(Normal(0, sqrt(h[i, 1])), r[i]) * Pt[i] + pdf(Normal(0, sqrt(h[i, 2])), r[i]) * (1 - Pt[i]));
        end
    else
        η = 1 / (2 + exp(-par[9]));
        @inbounds for i = 2:n
            h[i, :] .= ω .+ α .* r[i - 1]^2 + β .* h[i - 1, :];
            Pt[i] = probability_regime_given_time_it(p, q, sqrt.(h[i- 1, :]), r[i - 1], Pt[i - 1], η);
            log_lik[i - 1] = log(1/ sqrt(h[i, 1]) * Tstudent(r[i] / sqrt(h[i, 1]), η)* Pt[i]  + 1 / sqrt(h[i, 2]) * Tstudent(r[i] / sqrt(h[i, 2]), η) * (1 - Pt[i])); 
        end
    end
    return -sum(log_lik)/2;
end
##################################################
function klaassen_likelihood(r::Vector{Float64}, k::Int64, distri::String, par)
    n = length(r);
    h = Matrix{Float64}(undef, n, k + 1);
    Pt = Vector{Float64}(undef, n);
    log_lik = Vector{Float64}(undef, n - 1);
    A = Matrix{Float64}(undef, 2, 2);
    I2 = [1.0 0.0; 0.0 1.0];
    # Transformations
    ω = 10 ./ (1 .+ exp.(-par[1:2]));
    α = exp.(-par[3:4]) ./ (1 .+ exp.(-par[3:4]) .+ exp.(-par[5:6]));
    β = exp.(-par[5:6]) ./ (1 .+ exp.(-par[3:4]) .+ exp.(-par[5:6]));
    p = 1 ./(1 .+ exp(-par[7]));
    q = 1 ./(1 .+ exp(-par[8]));
    P = [p 1-q; 1-p q];
    # Likelihood
    Pt[1] = (1 - q) / (2 - p - q);    
    A[1, 1] = p * (α[1] + β[1]);
    A[1, 2] = (1 - p) * (α[1] + β[1]);
    A[2, 1] = (1 - q) * (α[2] + β[2]);
    A[2, 2] = q * (α[2] + β[2]);
    h[1, 1:k] .= inv(I2 - A) * ω;                    
    h[1, k + 1] = Pt[1] * h[1, 1] + (1 - Pt[1]) * h[1, 2];
    if distri == "norm"
        @inbounds for i = 2:n
            Pt[i] = probability_regime_given_time_n(p, q, sqrt.(h[i - 1, :]), r[i - 1], Pt[i - 1]);
            h[i, 1] = ω[1] + α[1] * r[i - 1]^2 + β[1] * (P[1,1] * pdf(Normal(0, sqrt(h[i - 1, 1])), r[i - 1]) * Pt[i - 1] * h[i - 1, 1] + P[2,1] * pdf(Normal(0, sqrt(h[i - 1, 2])), r[i - 1]) * (1 - Pt[i - 1]) * h[i - 1, 2])/(Pt[i] * (pdf(Normal(0, sqrt(h[i - 1, 1])), r[i - 1]) * Pt[i - 1] + pdf(Normal(0, sqrt(h[i - 1, 2])), r[i - 1]) * (1 - Pt[i - 1])));
            h[i, 2] = ω[2] + α[2] * r[i - 1]^2 + β[2] * (P[1,2] * pdf(Normal(0, sqrt(h[i - 1, 1])), r[i - 1]) * Pt[i - 1] * h[i - 1, 1] + P[2,2] * pdf(Normal(0, sqrt(h[i - 1, 2])), r[i - 1]) * (1 - Pt[i - 1]) * h[i - 1, 2])/((1 - Pt[i]) * (pdf(Normal(0, sqrt(h[i - 1, 1])), r[i - 1]) * Pt[i - 1] + pdf(Normal(0, sqrt(h[i - 1, 2])), r[i - 1]) * (1 - Pt[i - 1])));
            h[i, k + 1] = Pt[i] * h[i, 1] + (1 - Pt[i]) * h[i, 2];
            log_lik[i - 1] = log(pdf(Normal(0, sqrt(h[i, 1])), r[i]) * Pt[i] + pdf(Normal(0, sqrt(h[i, 2])), r[i]) * (1 - Pt[i]));
        end
    else
        η = 1 / (2 + exp(-par[9]));
        @inbounds for i = 2:n  
            Pt[i] = probability_regime_given_time_it(p, q, sqrt.(h[i - 1, :]), r[i - 1], Pt[i - 1], η);
            h[i, 1] = ω[1] + α[1] * r[i - 1]^2 + β[1] * (P[1,1] * (1/ sqrt(h[i - 1, 1]) * Tstudent(r[i - 1] / sqrt(h[i - 1, 1]), η)) * Pt[i - 1] * h[i - 1, 1] + P[2,1] * (1/ sqrt(h[i - 1, 2]) * Tstudent(r[i - 1] / sqrt(h[i - 1, 2]), η)) * (1 - Pt[i - 1]) * h[i - 1, 2])/(Pt[i] * ((1/ sqrt(h[i - 1, 1]) * Tstudent(r[i - 1] / sqrt(h[i - 1, 1]), η)) * Pt[i - 1] + (1/ sqrt(h[i - 1, 2]) * Tstudent(r[i - 1] / sqrt(h[i - 1, 2]), η)) * (1 - Pt[i - 1])));
            h[i, 2] = ω[2] + α[2] * r[i - 1]^2 + β[2] * (P[1,2] * (1/ sqrt(h[i - 1, 1]) * Tstudent(r[i - 1] / sqrt(h[i - 1, 1]), η)) * Pt[i - 1] * h[i - 1, 1] + P[2,2] * (1/ sqrt(h[i - 1, 2]) * Tstudent(r[i - 1] / sqrt(h[i - 1, 2]), η)) * (1 - Pt[i - 1]) * h[i - 1, 2])/((1 - Pt[i]) * ((1/ sqrt(h[i - 1, 1]) * Tstudent(r[i - 1] / sqrt(h[i - 1, 1]), η)) * Pt[i - 1] + (1/ sqrt(h[i - 1, 2]) * Tstudent(r[i - 1] / sqrt(h[i - 1, 2]), η)) * (1 - Pt[i - 1])));
            h[i, k + 1] = Pt[i] * h[i, 1] + (1 - Pt[i]) * h[i, 2];
            log_lik[i - 1] = log(1/ sqrt(h[i, 1]) * Tstudent(r[i] / sqrt(h[i, 1]), η)* Pt[i]  + 1 / sqrt(h[i, 2]) * Tstudent(r[i] / sqrt(h[i, 2]), η) * (1 - Pt[i]));
        end
    end
    return -sum(log_lik)/2;
end


##################################################
### Fit models
##################################################
function fit_gray(r::Vector{Float64}, k::Int64, par_ini, distri::String)
    sd = std(r);
    r = r/sd;
    ll = Vector{Float64}(undef, 5000);
    if distri == "norm"
        if isnothing(par_ini)
            par_ini = Matrix{Float64}(undef, 5000, 8);
            @try begin
                opt = optimize(par -> gray_likelihood(r, k, distri, par), [5, 2.5, 0.3, 0.1, 0.6, 4, 4, 3], iterations = 10000).minimizer;
                if isnan(gray_likelihood(r, k, distri, opt))
                    error("NaN")
                end
            @catch
                opt = NaN64;
                m = 0
                while typeof(opt) != Vector{Float64}
                    m = m + 1;
                    Random.seed!(1234 + m);
                    ω = rand(Uniform(0, 10), k);
                    α = rand(Uniform(0, 1), k);
                    β = [rand(Uniform(0, 1 - α[1])), rand(Uniform(0, 1 - α[2]))];
                    a = (α .+ β) ./ (1 .- α .- β);
                    tω = -log.(10 ./ ω .- 1);
                    tα = -log.((1 .+ a) .* α);
                    tβ = -log.((1 .+ a) .* β);
                    tp = rand(Uniform(1, 5), k);
                    par_random = [tω; tα; tβ; tp];
                    @try begin
                        opt = optimize(par -> gray_likelihood(r, k, distri, par), par_random, iterations = 10000).minimizer;
                    @catch
                        opt = NaN64;
                    end
                end
            end
            for i in 1:5000
                Random.seed!(i);
                ω = rand(Uniform(0, 10), k);
                α = rand(Uniform(0, 1), k);
                β = [rand(Uniform(0, 1 - α[1])), rand(Uniform(0, 1 - α[2]))];
                a = (α .+ β) ./ (1 .- α .- β);
                tω = -log.(10 ./ ω .- 1);
                tα = -log.((1 .+ a) .* α);
                tβ = -log.((1 .+ a) .* β);
                tp = rand(Uniform(1, 5), k);
                par_random = [tω; tα; tβ; tp];
                @try begin
                    gray_likelihood(r, k, distri, par_random);
                @catch
                    par_random = [-tω; tα; tβ; tp];
                    ll[i] = gray_likelihood(r, k, distri, par_random);
                    par_ini[i, :] .= par_random;
                @else 
                    ll[i] = gray_likelihood(r, k, distri, par_random);
                    par_ini[i, :] .= par_random;
                end
            end
            par_ini = par_ini[sortperm(ll),:];
            j = 1;
            l = 1;
            while l < 4
                @try begin
                    aux = optimize(par -> gray_likelihood(r, k, distri, par), par_ini[j,:], iterations = 10000).minimizer;
                    l = l + 1;
                    if gray_likelihood(r, k, distri, aux) < gray_likelihood(r, k, distri, opt)
                        opt = aux;
                    end
                @catch
                    l = l;
                end
                j = j + 1;
            end
            mle = param_transform(opt);
            mle[1:2] .= sd^2 .* mle[1:2];
        else
            opt = optimize(par -> gray_likelihood(r, k, distri, par), par_ini, iterations = 10000).minimizer;
            mle = param_transform(opt);
            mle[1:2] .= sd^2 .* mle[1:2];
        end
    elseif distri == "student"
        if isnothing(par_ini)
            par_ini = Matrix{Float64}(undef, 5000, 9);
            @try begin
                opt = optimize(par -> gray_likelihood(r, k, distri, par), [5, 2.5, 0.3, 0.1, 0.6, 4, 4, 3, 1.5], iterations = 10000).minimizer;
                if isnan(gray_likelihood(r, k, distri, opt))
                    error("NaN")
                end
            @catch
                opt = NaN64;
                m = 0
                while typeof(opt) != Vector{Float64}
                    m = m + 1;
                    Random.seed!(1234 + m);
                    ω = rand(Uniform(0, 10), k);
                    α = rand(Uniform(0, 1), k);
                    β = [rand(Uniform(0, 1 - α[1])), rand(Uniform(0, 1 - α[2]))];
                    a = (α .+ β) ./ (1 .- α .- β);
                    tω = -log.(10 ./ ω .- 1);
                    tα = -log.((1 .+ a) .* α);
                    tβ = -log.((1 .+ a) .* β);
                    tp = rand(Uniform(1, 5), k);
                    tν = rand(Uniform(-5, 3), 1);
                    par_random = [tω; tα; tβ; tp; tν];
                    @try begin
                        opt = optimize(par -> gray_likelihood(r, k, distri, par), par_random, iterations = 10000).minimizer;
                    @catch
                        opt = NaN64;
                    end
                end
            end
            for i in 1:5000
                Random.seed!(i);
                ω = rand(Uniform(0, 10), k);
                α = rand(Uniform(0, 1), k);
                β = [rand(Uniform(0, 1 - α[1])), rand(Uniform(0, 1 - α[2]))];
                a = (α .+ β) ./ (1 .- α .- β);
                tω = -log.(10 ./ ω .- 1);
                tα = -log.((1 .+ a) .* α);
                tβ = -log.((1 .+ a) .* β);
                tp = rand(Uniform(1, 5), k);
                tν = rand(Uniform(-5, 3), 1);
                par_random = [tω; tα; tβ; tp; tν];
                @try begin
                    gray_likelihood(r, k, distri, par_random);
                @catch
                    par_random = [-tω; tα; tβ; tp; tν];
                    ll[i] = gray_likelihood(r, k, distri, par_random);
                    par_ini[i, :] .= par_random;
                @else 
                    ll[i] = gray_likelihood(r, k, distri, par_random);
                    par_ini[i, :] .= par_random;
                end
            end
            par_ini = par_ini[sortperm(ll),:];
            j = 1;
            l = 1;
            while l < 4
                @try begin
                    aux = optimize(par -> gray_likelihood(r, k, distri, par), par_ini[j,:], iterations = 10000).minimizer;
                    l = l + 1;
                    if gray_likelihood(r, k, distri, aux) < gray_likelihood(r, k, distri, opt)
                        opt = aux;
                    end
                @catch
                    l = l;
                end
                j = j + 1;
            end
            mle = [param_transform(opt[1:8]); 2 + exp(-opt[9])];
            mle[1:2] .= sd^2 .* mle[1:2];
        else
            opt = optimize(par -> gray_likelihood(r, k, distri, par), par_ini, iterations = 10000).minimizer;
            mle = [param_transform(opt[1:8]); 2 + exp(-opt[9])];
            mle[1:2] .= sd^2 .* mle[1:2];
        end
    else
        mle = NaN64;
        println("Only Normal ('norm') and Student-T ('student') distributions are available")
    end
    return mle;
end
##################################################
function fit_haas(r::Vector{Float64}, k::Int64, par_ini, distri::String)
    sd = std(r);
    r = r/sd;
    ll = Vector{Float64}(undef, 5000);
    if distri == "norm"
        if isnothing(par_ini)
            par_ini = Matrix{Float64}(undef, 5000, 8);
            @try begin
                opt = optimize(par -> haas_likelihood(r, k, distri, par), [5, 2.5, 0.3, 0.1, 0.6, 4, 4, 3], iterations = 10000).minimizer;
                if isnan(haas_likelihood(r, k, distri, opt))
                    error("NaN")
                end
            @catch
                opt = NaN64;
                m = 0
                while typeof(opt) != Vector{Float64}
                    m = m + 1;
                    Random.seed!(1234 + m);
                    ω = rand(Uniform(0, 10), k);
                    α = rand(Uniform(0, 1), k);
                    β = [rand(Uniform(0, 1 - α[1])), rand(Uniform(0, 1 - α[2]))];
                    a = (α .+ β) ./ (1 .- α .- β);
                    tω = -log.(10 ./ ω .- 1);
                    tα = -log.((1 .+ a) .* α);
                    tβ = -log.((1 .+ a) .* β);
                    tp = rand(Uniform(1, 5), k);
                    par_random = [tω; tα; tβ; tp];
                    @try begin
                        opt = optimize(par -> haas_likelihood(r, k, distri, par), par_random, iterations = 10000).minimizer;
                    @catch
                        opt = NaN64;
                    end
                end
            end
            for i in 1:5000
                Random.seed!(i);
                ω = rand(Uniform(0, 10), k);
                α = rand(Uniform(0, 1), k);
                β = [rand(Uniform(0, 1 - α[1])), rand(Uniform(0, 1 - α[2]))];
                a = (α .+ β) ./ (1 .- α .- β);
                tω = -log.(10 ./ ω .- 1);
                tα = -log.((1 .+ a) .* α);
                tβ = -log.((1 .+ a) .* β);
                tp = rand(Uniform(1, 5), k);
                par_random = [tω; tα; tβ; tp];
                @try begin
                    haas_likelihood(r, k, distri, par_random);
                @catch
                    par_random = [-tω; tα; tβ; tp];
                    ll[i] = haas_likelihood(r, k, distri, par_random);
                    par_ini[i, :] .= par_random;
                @else 
                    ll[i] = haas_likelihood(r, k, distri, par_random);
                    par_ini[i, :] .= par_random;
                end
            end
            par_ini = par_ini[sortperm(ll),:];
            j = 1;
            l = 1;
            while l < 4
                @try begin
                    aux = optimize(par -> haas_likelihood(r, k, distri, par), par_ini[j,:], iterations = 10000).minimizer;
                    l = l + 1;
                    if haas_likelihood(r, k, distri, aux) < haas_likelihood(r, k, distri, opt)
                        opt = aux;
                    end
                @catch
                    l = l;
                end
                j = j + 1;
            end
            mle = param_transform(opt);
            mle[1:2] .= sd^2 .* mle[1:2];
        else
            opt = optimize(par -> haas_likelihood(r, k, distri, par), par_ini, iterations = 10000).minimizer;
            mle = param_transform(opt);
            mle[1:2] .= sd^2 .* mle[1:2];
        end
    elseif distri == "student"
        if isnothing(par_ini)
            par_ini = Matrix{Float64}(undef, 5000, 9);
            @try begin
                opt = optimize(par -> haas_likelihood(r, k, distri, par), [5, 2.5, 0.3, 0.1, 0.6, 4, 4, 3, 1.5], iterations = 10000).minimizer;
                if isnan(haas_likelihood(r, k, distri, opt))
                    error("NaN")
                end
            @catch
                opt = NaN64;
                m = 0
                while typeof(opt) != Vector{Float64}
                    m = m + 1;
                    Random.seed!(1234 + m);
                    ω = rand(Uniform(0, 10), k);
                    α = rand(Uniform(0, 1), k);
                    β = [rand(Uniform(0, 1 - α[1])), rand(Uniform(0, 1 - α[2]))];
                    a = (α .+ β) ./ (1 .- α .- β);
                    tω = -log.(10 ./ ω .- 1);
                    tα = -log.((1 .+ a) .* α);
                    tβ = -log.((1 .+ a) .* β);
                    tp = rand(Uniform(1, 5), k);
                    tν = rand(Uniform(-5, 3), 1);
                    par_random = [tω; tα; tβ; tp; tν];
                    @try begin
                        opt = optimize(par -> haas_likelihood(r, k, distri, par), par_random, iterations = 10000).minimizer;
                    @catch
                        opt = NaN64;
                    end
                end
            end
            for i in 1:5000
                Random.seed!(i);
                ω = rand(Uniform(0, 10), k);
                α = rand(Uniform(0, 1), k);
                β = [rand(Uniform(0, 1 - α[1])), rand(Uniform(0, 1 - α[2]))];
                a = (α .+ β) ./ (1 .- α .- β);
                tω = -log.(10 ./ ω .- 1);
                tα = -log.((1 .+ a) .* α);
                tβ = -log.((1 .+ a) .* β);
                tp = rand(Uniform(1, 5), k);
                tν = rand(Uniform(-5, 3), 1);
                par_random = [tω; tα; tβ; tp; tν];
                @try begin
                    haas_likelihood(r, k, distri, par_random);
                @catch
                    par_random = [-tω; tα; tβ; tp; tν];
                    ll[i] = haas_likelihood(r, k, distri, par_random);
                    par_ini[i, :] .= par_random;
                @else 
                    ll[i] = haas_likelihood(r, k, distri, par_random);
                    par_ini[i, :] .= par_random;
                end
            end
            par_ini = par_ini[sortperm(ll),:];
            j = 1;
            l = 1;
            while l < 4
                @try begin
                    aux = optimize(par -> haas_likelihood(r, k, distri, par), par_ini[j,:], iterations = 10000).minimizer;
                    l = l + 1;
                    if haas_likelihood(r, k, distri, aux) < haas_likelihood(r, k, distri, opt)
                        opt = aux;
                    end
                @catch
                    l = l;
                end
                j = j + 1;
            end
            mle = [param_transform(opt[1:8]); 2 + exp(-opt[9])];
            mle[1:2] .= sd^2 .* mle[1:2];
        else
            opt = optimize(par -> haas_likelihood(r, k, distri, par), par_ini, iterations = 10000).minimizer;
            mle = [param_transform(opt[1:8]); 2 + exp(-opt[9])];
            mle[1:2] .= sd^2 .* mle[1:2];
        end
    else
        mle = NaN64;
        println("Only Normal ('norm') and Student-T ('student') distributions are available")
    end
    return mle;
end
##################################################
function fit_klaassen(r::Vector{Float64}, k::Int64, par_ini, distri::String)
    sd = std(r);
    r = r/sd;
    ll = Vector{Float64}(undef, 5000);
    if distri == "norm"
        if isnothing(par_ini)
            par_ini = Matrix{Float64}(undef, 5000, 8);
            @try begin
                opt = optimize(par -> klaassen_likelihood(r, k, distri, par), [5, 2.5, 0.3, 0.1, 0.6, 4, 4, 3], iterations = 10000).minimizer;
                if isnan(klaassen_likelihood(r, k, distri, opt))
                    error("NaN")
                end
            @catch
                opt = NaN64;
                m = 0
                while typeof(opt) != Vector{Float64}
                    m = m + 1;
                    Random.seed!(1234 + m);
                    ω = rand(Uniform(0, 10), k);
                    α = rand(Uniform(0, 1), k);
                    β = [rand(Uniform(0, 1 - α[1])), rand(Uniform(0, 1 - α[2]))];
                    a = (α .+ β) ./ (1 .- α .- β);
                    tω = -log.(10 ./ ω .- 1);
                    tα = -log.((1 .+ a) .* α);
                    tβ = -log.((1 .+ a) .* β);
                    tp = rand(Uniform(1, 5), k);
                    par_random = [tω; tα; tβ; tp];
                    @try begin
                        opt = optimize(par -> klaassen_likelihood(r, k, distri, par), par_random, iterations = 10000).minimizer;
                    @catch
                        opt = NaN64;
                    end
                end
            end
            for i in 1:5000
                Random.seed!(i);
                ω = rand(Uniform(0, 10), k);
                α = rand(Uniform(0, 1), k);
                β = [rand(Uniform(0, 1 - α[1])), rand(Uniform(0, 1 - α[2]))];
                a = (α .+ β) ./ (1 .- α .- β);
                tω = -log.(10 ./ ω .- 1);
                tα = -log.((1 .+ a) .* α);
                tβ = -log.((1 .+ a) .* β);
                tp = rand(Uniform(1, 5), k);
                par_random = [tω; tα; tβ; tp];
                @try begin
                    klaassen_likelihood(r, k, distri, par_random);
                @catch
                    par_random = [-tω; tα; tβ; tp];
                    ll[i] = klaassen_likelihood(r, k, distri, par_random);
                    par_ini[i, :] .= par_random;
                @else 
                    ll[i] = klaassen_likelihood(r, k, distri, par_random);
                    par_ini[i, :] .= par_random;
                end
            end
            par_ini = par_ini[sortperm(ll),:];
            j = 1;
            l = 1;
            while l < 4
                @try begin
                    aux = optimize(par -> klaassen_likelihood(r, k, distri, par), par_ini[j,:], iterations = 10000).minimizer;
                    l = l + 1;
                    if klaassen_likelihood(r, k, distri, aux) < klaassen_likelihood(r, k, distri, opt)
                        opt = aux;
                    end
                @catch
                    l = l;
                end
                j = j + 1;
            end
            mle = param_transform(opt);
            mle[1:2] .= sd^2 .* mle[1:2];
        else
            opt = optimize(par -> klaassen_likelihood(r, k, distri, par), par_ini, iterations = 10000).minimizer;
            mle = param_transform(opt);
            mle[1:2] .= sd^2 .* mle[1:2];
        end
    elseif distri == "student"
        if isnothing(par_ini)
            par_ini = Matrix{Float64}(undef, 5000, 9);
            @try begin
                opt = optimize(par -> klaassen_likelihood(r, k, distri, par), [5, 2.5, 0.3, 0.1, 0.6, 4, 4, 3, 1.5], iterations = 10000).minimizer;
                if isnan(klaassen_likelihood(r, k, distri, opt))
                    error("NaN")
                end
            @catch
                opt = NaN64;
                m = 0
                while typeof(opt) != Vector{Float64}
                    m = m + 1;
                    Random.seed!(1234 + m);
                    ω = rand(Uniform(0, 10), k);
                    α = rand(Uniform(0, 1), k);
                    β = [rand(Uniform(0, 1 - α[1])), rand(Uniform(0, 1 - α[2]))];
                    a = (α .+ β) ./ (1 .- α .- β);
                    tω = -log.(10 ./ ω .- 1);
                    tα = -log.((1 .+ a) .* α);
                    tβ = -log.((1 .+ a) .* β);
                    tp = rand(Uniform(1, 5), k);
                    tν = rand(Uniform(-5, 3), 1);
                    par_random = [tω; tα; tβ; tp; tν];
                    @try begin
                        opt = optimize(par -> klaassen_likelihood(r, k, distri, par), par_random, iterations = 10000).minimizer;
                    @catch
                        opt = NaN64;
                    end
                end
            end
            for i in 1:5000
                Random.seed!(i);
                ω = rand(Uniform(0, 10), k);
                α = rand(Uniform(0, 1), k);
                β = [rand(Uniform(0, 1 - α[1])), rand(Uniform(0, 1 - α[2]))];
                a = (α .+ β) ./ (1 .- α .- β);
                tω = -log.(10 ./ ω .- 1);
                tα = -log.((1 .+ a) .* α);
                tβ = -log.((1 .+ a) .* β);
                tp = rand(Uniform(1, 5), k);
                tν = rand(Uniform(-5, 3), 1);
                par_random = [tω; tα; tβ; tp; tν];
                @try begin
                    klaassen_likelihood(r, k, distri, par_random);
                @catch
                    par_random = [-tω; tα; tβ; tp; tν];
                    ll[i] = klaassen_likelihood(r, k, distri, par_random);
                    par_ini[i, :] .= par_random;
                @else 
                    ll[i] = klaassen_likelihood(r, k, distri, par_random);
                    par_ini[i, :] .= par_random;
                end
            end
            par_ini = par_ini[sortperm(ll),:];
            j = 1;
            l = 1;
            while l < 4
                @try begin
                    aux = optimize(par -> klaassen_likelihood(r, k, distri, par), par_ini[j,:], iterations = 10000).minimizer;
                    l = l + 1;
                    if klaassen_likelihood(r, k, distri, aux) < klaassen_likelihood(r, k, distri, opt)
                        opt = aux;
                    end
                @catch
                    l = l;
                end
                j = j + 1;
            end
            mle = [param_transform(opt[1:8]); 2 + exp(-opt[9])];
            mle[1:2] .= sd^2 .* mle[1:2];
        else
            opt = optimize(par -> klaassen_likelihood(r, k, distri, par), par_ini, iterations = 10000).minimizer;
            mle = [param_transform(opt[1:8]); 2 + exp(-opt[9])];
            mle[1:2] .= sd^2 .* mle[1:2];
        end
    else
        mle = NaN64;
        println("Only Normal ('norm') and Student-T ('student') distributions are available")
    end
    return mle;
end

