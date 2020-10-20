# necessary packages #

#using Pkg
#Pkg.add("Distances")
using Distributions
using Random
using Distances
using LinearAlgebra
using SparseArrays
using IterativeSolvers
using ProgressMeter
using JLD2
using MCMCDiagnostics

include("../../../util2.j")

@load "../data/sim3data.jld";

# priors #
p = size(X_ord)[2]; q = size(Y_ord)[2]; K = 1;
μβ = spzeros(p, q); inv_Vr = spzeros(p, p);
μΛ = spzeros(K, q); inv_VΛ = spzeros(K, K)
aΣ = 2; bΣ = fill(1.0, q);
inv_Lr = spzeros(p, p); inv_LΛ = spzeros(K, K);
aϕ = fill(2.0, K); bϕ = fill(33 / (2 * sqrt(2)), K);

# Some data preparations #

dim_invD = sum(index_S_M);
invD_yind = 1:dim_invD; invD_xind = 1:dim_invD;
Xtilde_indy_up = vcat([S .+ (ind - 1) * N for ind in 1:K]...);
nsam = length(perm_ind) + (K * n);

# preallocation #
μ_m = [Array{Float64, 2}(undef, length(M_ind[i]), q) for i in 1:q];
nIndx = length(NN.nnIndx);
A = [Array{Float64}(undef, nIndx) for i in 1:K];
D = [Array{Float64}(undef, n) for i in 1:K];
I_A = [spzeros(n, n) for i in 1:K];
A_new = [Array{Float64}(undef, nIndx) for i in 1:K];
D_new = [Array{Float64}(undef, n) for i in 1:K];
I_A_new = [spzeros(n, n) for i in 1:K];
Ystar = vcat(Y_ord[S, :], inv_Lr * μβ, inv_LΛ * μΛ);             # will be updated after imputing missing response
Xstar = vcat([X_ord[S, :] spzeros(n, K)], [inv_Lr spzeros(p, K)], 
    [spzeros(K, p) inv_LΛ]);      
bstar = fill(0.0, q); astar = aΣ + 0.5 * n;
μγstar = vcat(μβ, μΛ); invVγstar = fill(0.0, p + K, p + K);
Y_Xm = spzeros(n + p + K, q);
nsam = length(perm_ind) + (K * n);
Ytilde =  Array{Float64}(undef, nsam);
Xtilde = SparseMatrixCSC{Float64,Int64};
lll = fill(1.0, (n, 1)); 

# Preallocation for MCMC samples and Initalization #
N_sam = 10000;
N_pre_burn = Integer(trunc(0.75 * N_sam));
N_pre_adapt = Integer(trunc(0.25 * N_sam));
N_after_burn = N_sam - N_pre_burn;
ω_incp_sam = Array{Float32, 2}(undef, n, q);

ω_incp_sam_mean = fill(0.0, n, q);
ω_incp_sam_var = fill(0.0, n, q);
Y_m_sam_mean = [fill(0.0, length(M_ind[i])) for i in 1:length(M_ind)];
Y_m_sam_var = [fill(0.0, length(M_ind[i])) for i in 1:length(M_ind)];

F_sam = Array{Float64, 2}(undef, n, K);
Y_m_sam =  [Array{Float64, 1}(undef, length(M_ind[i])) for i in 1:q];
A_sam = Array{Float64, 2}(undef, N_sam, K); # acceptance rate
lh_old = 1; lh_new = 1;     # record the likelihood for updating ranges

ϕ_sam = Array{Float64, 2}(undef, K, N_sam + 1);

γ_sam = vcat(fill(0.0, p, q), fill(0.001, K, q));
γ_sam[(p + 1):(p + K), 1:K] = γ_sam[(p + 1):(p + K), 1:K] + I;
Σ_sam = fill(0.2, q);
ω_cov_sam = fill(0.0, q, q);
ϕ_sam[:, 1] = fill(6.0, K);

#precond_D = Array{Float64, 1}(undef, K * n);
inv_sqrt_Σ_diag = Array{Float64, 1}(undef, q);
RWM_scale = fill(0.1, K);

using DelimitedFiles
writedlm("../results/K1/γ_sam.csv", vcat(fill(0.0, 1, q), γ_sam), ", ");
writedlm("../results/K1/Σ_sam.csv", vcat(0.0, Σ_sam), ", ");
writedlm("../results/K1/ω_cov_sam.csv", [fill(0.0, q)], ", ");
writedlm("../results/K1/ϕ_sam.csv", vcat(0.0, ϕ_sam[:, 1]), ", ");
writedlm("../results/K1/A_sam.csv", 0.0, ", ");

# for loop for MCMC chain #
Random.seed!(123);
prog = Progress(N_sam, 1, "Computing initial pass...", 50)
for l in 1:N_sam
    # Build the matrix D_Sigma_o^{1/2} #
    inv_sqrt_Σ_diag = 1 ./ (sqrt.(Σ_sam));
    invD_ele = [x for ind in 1:q for x in fill(inv_sqrt_Σ_diag[ind], length(S_ind[ind]))];
    invD = sparse(invD_xind, invD_yind, invD_ele);
                
    if l == 1
        for k in 1:K
            getAD(coords_ord[S, :], NN.nnIndx, NN.nnDist, NN.nnIndxLU, ϕ_sam[k, l], 0.5, 
                            A[k], D[k]);
            I_A[k] = sparse(nnIndx_row, nnIndx_col, vcat(-A[k], ones(n)));
        end
    end
    Ytilde = vcat(invD * vcat([Y_ord[S_ind[ind], ind] - X_ord[S_ind[ind], :] * 
                            γ_sam[1:p, ind] for ind in 1:q]...), zeros(K * n));
    Xtilde = vcat(invD * kron(sparse(transpose(γ_sam[(p + 1):(p + K), :])), 
                            sparse(1:N, 1:N, ones(N)))[obs_ind, Xtilde_indy_up],
             blockdiag([Diagonal(1 ./ sqrt.(D[ind])) * I_A[ind] for ind in 1:K]...));
                  
    # use LSMR to generate sample of F #                    
    F_sam = reshape(lsmr(Xtilde, collect(Ytilde) + rand(Normal(), nsam)), :, K);
                
    # impute missing response  over S#
    Xstar[1:n, (p + 1):(p + K)] = F_sam;        # update matrix Xstar with F
    if(l > N_pre_burn) # only save ω_incp_sam after burn-in
        ω_incp_sam = F_sam * γ_sam[(p + 1):(p + K), :] + lll * transpose(γ_sam[1, :]); 
        ω_incp_sam_mean = ω_incp_sam_mean + (ω_incp_sam ./ N_after_burn);
        ω_incp_sam_var = ω_incp_sam_var + ((ω_incp_sam.^2) ./ N_after_burn);  
        ω_cov_sam = cov(ω_incp_sam);
    else
        ω_cov_sam = cov(F_sam * γ_sam[(p + 1):(p + K), :]);
    end                    
    io1 = open("../results/K1/ω_cov_sam.csv", "a" ); # covariance of latent process
    writedlm(io1, ω_cov_sam, ", ");
    close(io1);            
           
    # impute missing response  over S#            
    for ind in 1:q
        Y_m_sam[ind] = Xstar[M_Sind[ind], :] * γ_sam[:, ind]+   
            rand(Normal(0.0, sqrt(Σ_sam[ind])), length(M_ind[ind]));
    end
         
    if (l > N_pre_burn)  # only save imputed Y after burn-in
        for ind in 1:q
            Y_m_sam_mean[ind] = Y_m_sam_mean[ind] + (Y_m_sam[ind] ./ N_after_burn);
            Y_m_sam_var[ind] = Y_m_sam_var[ind] + ((Y_m_sam[ind].^2) ./ N_after_burn);
        end
    end
                
                        
    # use MNIW to sample γ Σ #
    for ind in 1:q
        Ystar[M_Sind[ind], ind] = Y_m_sam[ind]   # update Ystar with imputed response
    end
    
    invVγstar = cholesky(Xstar'Xstar);
    mul!(μγstar, transpose(Xstar), Ystar); μγstar = invVγstar.U \ (invVγstar.L \ μγstar);
    Y_Xm = Ystar - Xstar * μγstar;      # maybe improve?
    bstar = [bΣ[ind] + 0.5 * (norm(Y_Xm[:, ind])^2) for ind in 1:q];
    Σ_sam = [rand(InverseGamma(astar, bstar[ind]), 1)[1] for ind in 1:q];          # sample Σ
    γ_sam = (invVγstar.U \ reshape(rand(Normal(), (p + K) * q), (p + K), q)) * 
                    Diagonal(sqrt.(Σ_sam)) + μγstar;          # sample γ    
    io2 = open("../results/K1/Σ_sam.csv", "a" );
    writedlm(io2, Σ_sam, ", ");
    close(io2); 
    io3 = open("../results/K1/γ_sam.csv", "a" );
    writedlm(io3, γ_sam, ", ");
    close(io3);
                
    # use adaptive metropolis-hasting to update range
    if l > 3 && l < N_pre_adapt
        RWM_scale = [sqrt(2.38^2 * var(ϕ_sam[i, Integer(floor(l / 2)):l], 
                                corrected=true) * 0.95^2 + 0.05^2 * 0.1^2) for i in 1:K];
    end
    # use metropolis-hasting to update range
    ϕ_sam[:, l + 1] = [ϕ_sam[i, l] + RWM_scale[i] * rand(Normal(), 1)[1] for i in 1:K]; # propose next sample point
  
    for i in 1:K
        if ϕ_sam[i, l + 1] > 0.0
            lh_old = -0.5 * (sum(log.(D[i])) + norm((I_A[i] * F_sam[:, i]) ./ sqrt.(D[i]))^2) + 
               loglikelihood(Gamma(aϕ[i], bϕ[i]), [ϕ_sam[i, l]]);
            getAD(coords_ord[S, :], NN.nnIndx, NN.nnDist, NN.nnIndxLU, ϕ_sam[i, l + 1], 
                0.5, A_new[i], D_new[i]);
            I_A_new[i] = sparse(nnIndx_row, nnIndx_col, vcat(-A_new[i], ones(n)));
            lh_new = -0.5 * (sum(log.(D_new[i]))  + norm((I_A_new[i] * F_sam[:, i]) ./ sqrt.(D_new[i]))^2) +
                   loglikelihood(Gamma(aϕ[i], bϕ[i]), [ϕ_sam[i, l + 1]]);     
            A_sam[l, i] = min(exp(lh_new - lh_old), 1.0);
            if rand(1)[1] < A_sam[l, i]
                I_A[i] = copy(I_A_new[i]); D[i] = copy(D_new[i]);        # update and update the corresponding I_A D
            else 
                ϕ_sam[i, l + 1] = ϕ_sam[i, l];
            end
        else 
            A_sam[l, i] = 0.0;
            ϕ_sam[:, l + 1] = ϕ_sam[:, l];   
        end
    end                   
    
    io4 = open("../results/K1/ϕ_sam.csv", "a" );
    writedlm(io4, ϕ_sam[:, l + 1], ", ");
    close(io4); 
    io5 = open("../results/K1/A_sam.csv", "a" );
    writedlm(io5, A_sam[l, :], ", ");
    close(io5);
                
    next!(prog); # monitor the progress
end
ω_incp_sam_var = (ω_incp_sam_var - ω_incp_sam_mean.^2) * (N_after_burn / (N_after_burn - 1));
Y_m_sam_var = [(Y_m_sam_var[ind] - Y_m_sam_mean[ind].^2) * 
               (N_after_burn / (N_after_burn - 1)) for ind in 1:q];

#load data
using CSV
γ_sam = convert(Matrix{Float64}, CSV.read("../results/K1/γ_sam.csv"));
ind_γ_sam = 1: (p + K) :((p + K) * N_sam + 1);
Σ_sam = convert(Matrix{Float64}, CSV.read("../results/K1/Σ_sam.csv"));
ind_Σ_sam = 1: q :(q * N_sam + 1);
ω_cov_sam = convert(Matrix{Float64}, CSV.read("../results/K1/ω_cov_sam.csv"));
ind_ω_cov_sam = 1: q :(q * (N_sam - 1) + 1);
            
covω = cov(ω_ord[S, :]);
corω = cor(ω_ord[S, :]);
ω_cov_pos_sam_mean = [mean(ω_cov_sam[ind_ω_cov_sam .+ (i - 1), j][(N_pre_burn + 1):N_sam]) 
    for i in 1:q, j in 1:q];
ω_corr_sam = [(Diagonal([1 / sqrt(ω_cov_sam[ind_ω_cov_sam[l] .+ (i - 1), i]) for i in 1:q]) * 
    ω_cov_sam[ind_ω_cov_sam[l] .+ (1:q) .- 1, 1:q] * 
        Diagonal([1 / sqrt(ω_cov_sam[ind_ω_cov_sam[l] .+ (i - 1), i]) for i in 1:q])) for l in 1:N_sam];
ω_corr_sam_mean = [mean([ω_corr_sam[i][j , k] for i in (N_pre_burn + 1):N_sam]) 
    for j in 1:q, k in 1:q];
# CVG-slope #
count_slope = 0.0
for i in 2:p
    for j in 1:q
        if ((quantile(γ_sam[ind_γ_sam .+ (i - 1), j][(N_pre_burn + 1):(N_sam + 1)], [0.025])[1] < 
                β[i, j] ) && (quantile(
                        γ_sam[ind_γ_sam .+ (i - 1), j][(N_pre_burn + 1):(N_sam + 1)], [0.975])[1] > 
                β[i, j] ))
            count_slope = count_slope + 1.0;
        end
        
    end
end
count_slope
# ESS-slope #
ESS_slope = fill(0.0, (p - 1), q)
MCSE_slope = fill(0.0, (p - 1), q)
for i in 2:p
    for j in 1:q
        ESS_slope[i - 1, j] = 
            effective_sample_size(γ_sam[ind_γ_sam .+ (i - 1), j][(N_pre_burn + 1):(N_sam + 1)]);
        MCSE_slope[i - 1, j] = std(γ_sam[ind_γ_sam .+ (i - 1), j][(N_pre_burn + 1):(N_sam + 1)]) / 
            sqrt(ESS_slope[i - 1, j]);
    end
end
print("slope ESS&MCSE:", [minimum(ESS_slope), maximum(MCSE_slope)], "\n")
ESS_ω_cov = fill(0.0, q, q)
MCSE_ω_cov = fill(0.0, q, q)
for j in 1:q
    for i in 1:q
        ESS_ω_cov[i, j] = 
            effective_sample_size(ω_cov_sam[ind_ω_cov_sam .+ (i - 1), j][(N_pre_burn + 1):N_sam]);
        MCSE_ω_cov[i, j] = std(ω_cov_sam[ind_ω_cov_sam .+ (i - 1), j][(N_pre_burn + 1):N_sam]) / 
            sqrt(ESS_ω_cov[i, j]);
    end
end
print("ω_cov ESS&MCSE:", [minimum(ESS_ω_cov), maximum(MCSE_ω_cov)], "\n")
ESS_ω_corr = fill(0.0, q, q)
MCSE_ω_corr = fill(0.0, q, q)
for j in 1:q
    for i in 1:q
        ESS_ω_corr[i, j] = 
            effective_sample_size([ω_corr_sam[k][i, j] for k in (N_pre_burn + 1):N_sam]);
        MCSE_ω_corr[i, j] = std([ω_corr_sam[k][i, j] for k in (N_pre_burn + 1):N_sam]) / 
            sqrt(ESS_ω_corr[i, j]);
    end
end
print("ω_corr ESS&MCSE:",[minimum(ESS_ω_corr), maximum(MCSE_ω_corr)], "\n")
ESS_Σ = fill(0.0, q)
MCSE_Σ = fill(0.0, q)
for i in 1:q
    ESS_Σ[i] = 
        effective_sample_size(Σ_sam[ind_Σ_sam .+ (i - 1)][(N_pre_burn + 1):N_sam]);
    MCSE_Σ[i] = std(Σ_sam[ind_Σ_sam .+ (i - 1)][(N_pre_burn + 1):N_sam]) / sqrt(ESS_Σ[i]);
end
print("Σ ESS&MCSE:", [minimum(ESS_Σ), maximum(MCSE_Σ)], "\n")
            
# CVL #
count_ω_incp = fill(0.0, q);
for j in 1:q
    for i in 1:n
        count_ω_incp[j] = count_ω_incp[j] + 
        (((ω_incp_sam_mean[i, j] - 1.96 * sqrt(ω_incp_sam_var[i, j])) < ω_incp_obs[S[i], j]) && 
            ((ω_incp_sam_mean[i, j] + 1.96 * sqrt(ω_incp_sam_var[i, j])) > ω_incp_obs[S[i], j]))
    end
end
count_ω_incp;
print("count_ω_incp", round.(count_ω_incp ./ n, digits = 5), sum(count_ω_incp) / (q*n), "\n")    
            
# CVG #
count_Y_M = fill(0.0, q);
Y_m_pos_qt = [Array{Float64, 2}(undef, length(M_ind[ind]), 3) for ind in 1:q];

for i in 1:q
    for j in 1:length(M_ind[i])
        count_Y_M[i] = count_Y_M[i] + (((Y_m_sam_mean[i][j] - 
                1.96 * sqrt(Y_m_sam_var[i][j])) < Y_ord[M_ind[i][j], i]) && 
         ((Y_m_sam_mean[i][j] + 1.96 * sqrt(Y_m_sam_var[i][j])) > 
                Y_ord[M_ind[i][j], i]))
    end
end
print("count_Y_M: ", round.(count_Y_M ./ 200, digits = 4), 
                round.(mean(count_Y_M ./ 200), digits = 4), "\n");
# RMSPE #
MSPE = (sum([sum((Y_m_sam_mean[i] - Y_ord[S[M_Sind[i]], i]).^2) for i in 1:q])) / 
    (sum([length(M_Sind[i]) for i in 1:q]))
RMSPE = sqrt(MSPE); 
print("RMSPE:", round.([sqrt((sum((Y_m_sam_mean[i] - Y_ord[S[M_Sind[i]], i]).^2)) / 
        length(M_Sind[i])) for i in 1:q], digits = 4), round(RMSPE, digits = 4), "\n");
@save "../results/K1/Factor_mean_var_K1.jld" ESS_slope MCSE_slope ESS_ω_corr MCSE_ω_corr ESS_ω_cov MCSE_ω_cov ESS_Σ MCSE_Σ count_slope count_ω_incp count_Y_M RMSPE ω_incp_sam_mean ω_incp_sam_var Y_m_sam_mean Y_m_sam_var Y_ord S_ind M_Sind K p q N_sam N_pre_burn
