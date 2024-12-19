##################################################
###            RSGARCH: Forecasts             ####
##################################################
using Distributions, Optim, Statistics, StatsFuns, Random, SpecialFunctions, TryCatch, DelimitedFiles, StatsBase, QuadGK, CSV, DataFrames, LinearAlgebra, Kronecker
include("utils.jl")
include("DGP.jl")
include("MaximumLikelihood.jl")
include("Forecast.jl")


# Import Data
prices = DataFrame(CSV.File("./Data/ETHUSDT_1d.csv", header = 1, delim=","));
select!(prices, [:OpenTime, :Close]);
rename!(prices, Symbol.(["Date","Price"]));
dropmissing!(prices);
prices.returns = [missing; 100*(log.(prices.Price[2:end]) -log.(prices.Price[1:end-1]))];
dropmissing!(prices);



# Settings
InS = 1000;
Tot = size(prices)[1];
OoS = Tot - InS;
VaR_1 = Matrix{Float64}(undef, OoS, 6);
VaR_2 = Matrix{Float64}(undef, OoS, 6);
ES_1  = Matrix{Float64}(undef, OoS, 6);
ES_2  = Matrix{Float64}(undef, OoS, 6);
r_oos  = Vector{Float64}(undef, OoS);



# Out-of-sample VaR and ES
k = 2;
for i in 1:OoS
    println(i)
    r = prices.returns[i:(InS + i - 1)];
    r_oos[i] = prices.returns[InS + i];
    μ = mean(r);
    r = r .- μ;
    Random.seed!(1234 + i);
    θ = fit_gray(r, k, nothing, "norm");
    (h, Pt, s)  = fore_gray(r, k, θ, "norm");
    (VaR_1[i,1], ES_1[i,1]) = μ .+ var_es_rsgarch(0.01, Pt[end], 1 - Pt[end], sqrt(h[1]), sqrt(h[2]), "norm");
    (VaR_2[i,1], ES_2[i,1]) = μ .+ var_es_rsgarch(0.025, Pt[end], 1 - Pt[end], sqrt(h[1]), sqrt(h[2]), "norm");

    θ = fit_gray(r, k, nothing, "student"); 
    (h, Pt, s)  = fore_gray(r, k, θ, "student");
    (VaR_1[i,2], ES_1[i,2]) = μ .+ var_es_rsgarch(0.01, Pt[end], 1 - Pt[end], sqrt(h[1]), sqrt(h[2]), "student", θ[9]);
    (VaR_2[i,2], ES_2[i,2]) = μ .+ var_es_rsgarch(0.025, Pt[end], 1 - Pt[end], sqrt(h[1]), sqrt(h[2]), "student", θ[9]);

    θ = fit_klaassen(r, k, nothing, "norm");
    (h, Pt, s)  = fore_klaassen(r, k, θ, "norm");
    (VaR_1[i,3], ES_1[i,3]) = μ .+ var_es_rsgarch(0.01, Pt[end], 1 - Pt[end], sqrt(h[1]), sqrt(h[2]), "norm");
    (VaR_2[i,3], ES_2[i,3]) = μ .+ var_es_rsgarch(0.025, Pt[end], 1 - Pt[end], sqrt(h[1]), sqrt(h[2]), "norm");

    θ = fit_klaassen(r, k, nothing, "student");
    (h, Pt, s)  = fore_klaassen(r, k, θ, "student");
    (VaR_1[i,4], ES_1[i,4]) = μ .+ var_es_rsgarch(0.01, Pt[end], 1 - Pt[end], sqrt(h[1]), sqrt(h[2]), "student", θ[9]);
    (VaR_2[i,4], ES_2[i,4]) = μ .+ var_es_rsgarch(0.025, Pt[end], 1 - Pt[end], sqrt(h[1]), sqrt(h[2]), "student", θ[9]);

    θ = fit_haas(r, k, nothing, "norm");
    (h, Pt, s)  = fore_haas(r, k, θ, "norm");
    (VaR_1[i,5], ES_1[i,5]) = μ .+ var_es_rsgarch(0.01, Pt[end], 1 - Pt[end], sqrt(h[1]), sqrt(h[2]), "norm");
    (VaR_2[i,5], ES_2[i,5]) = μ .+ var_es_rsgarch(0.025, Pt[end], 1 - Pt[end], sqrt(h[1]), sqrt(h[2]), "norm");

    θ = fit_haas(r, k, nothing, "student");
    (h, Pt, s)  = fore_haas(r, k, θ, "student");
    (VaR_1[i,6], ES_1[i,6]) = μ .+ var_es_rsgarch(0.01, Pt[end], 1 - Pt[end], sqrt(h[1]), sqrt(h[2]), "student", θ[9]);
    (VaR_2[i,6], ES_2[i,6]) = μ .+ var_es_rsgarch(0.025, Pt[end], 1 - Pt[end], sqrt(h[1]), sqrt(h[2]), "student", θ[9]);
end

writedlm("VaR1_ETH_1000.csv",  VaR_1, ',')
writedlm("VaR2_ETH_1000.csv",  VaR_2, ',')
writedlm("ES1_ETH_1000.csv",  ES_1, ',')
writedlm("ES2_ETH_1000.csv",  ES_2, ',')
writedlm("r_oos_ETH_1000.csv",  r_oos, ',')