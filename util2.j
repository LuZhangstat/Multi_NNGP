# build nearest neighbor for observed location set S #
function BuildNN(coords, m, nThreads = 1.0, brute = 0.0)
    # coords: n by 2 array
    # m:      number of neighbors
    n = size(coords)[1]
    # coords_fit = [coords[1, :] coords[2, :]]
    nIndx = trunc(Cint,(1 + m) / 2 * m + (n - m - 1) * m);
    nnIndx = zeros(Cint, nIndx);
    nnDist = zeros(Float64, nIndx);
    ##first of nnIndxLU column holds the nnIndx index for the i-th location and the second columns holds the number of neighbors the i-th location has (the second column is a bit of a waste but simplifies my life in the spNNGP)
    nnIndxLU = zeros(Cint, 2 * n);
    ##Note I need to convert n and m to pointers only because the R example of .C wants only pointers so that's how I worte the mkNNIndx and mkNNIndxTree0
    nPtr = Cint[n];
    mPtr = Cint[m];
    nThreadsPtr = Cint[nThreads];   

    if brute == 1.0
        ccall((:mkNNIndx, "../../../julia-R-nn-ccall2/nn.so"), Cvoid, (Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble}, 
            Ptr{Cint}, Ptr{Cdouble}, Ptr{Cint}, Ptr{Cint}), nPtr, mPtr, coords, nnIndx, nnDist, nnIndxLU, nThreadsPtr)
    else
        ccall((:mkNNIndxCB, "../../../julia-R-nn-ccall2/nn.so"), Cvoid, (Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble}, 
            Ptr{Cint}, Ptr{Cdouble}, Ptr{Cint}, Ptr{Cint}), nPtr, mPtr, coords, nnIndx, nnDist, nnIndxLU, nThreadsPtr)
    end
    
    nnIndx = nnIndx .+ 1 
    nnIndxLU = nnIndxLU[1:(n + 1)]  .+ 1
    nnIndxLU[n + 1] = nIndx + 1
    return (nnIndx = nnIndx, nnDist = nnDist, nnIndxLU = nnIndxLU[1:(n + 1)])
end  

function getAD(coords, nnIndx, nnDist, nnIndxLU, ϕ, ν, A, D)
    # coords: n by 2 array,
    # nnIndx, nnDist, nnIndxLU output of the mkNNIndx function
    
    # return: A D
    # A = []
    # D = ones(size(coords)[1])
    
    # initialize memory
    NNdistM = []
    NNdistv = []
    n = length(nnIndxLU) - 1;
    if nnIndxLU[1] == nnIndxLU[2]
        D[1] = 1.0
        l = 2
    else 
        l = 1
    end
    for i in l:n
        if ν == 0.5
            NNdistM = cholesky(exp.(-ϕ .* pairwise(Euclidean(), 
                        coords[nnIndx[nnIndxLU[i]:(nnIndxLU[i + 1] - 1)], :], 
                        dims = 1)))
            NNdistv = NNdistM.L \ (exp.(-ϕ .* nnDist[nnIndxLU[i]:(nnIndxLU[i + 1] - 1)]))
            D[i] = 1.0 - NNdistv⋅NNdistv
            A[nnIndxLU[i]:(nnIndxLU[i + 1] - 1)] = NNdistM.U \ NNdistv
        end
    end
    # return (A = A, D = D)
end

function getAD_collapse(coords, nnIndx, nnDist, nnIndxLU, ϕ, ν, α, A, D)
    # coords: n by 2 array,
    # nnIndx, nnDist, nnIndxLU output of the mkNNIndx function
    # ϕ, ν, α: covariance parameter set
    # hold: whether AD is for holded locations or not
    
    # return: A D
    # A = []
    # D = ones(size(coords)[1])
    
    # initialize memory
    NNdistM = []
    NNdistv = []
    n = length(nnIndxLU) - 1;
    if nnIndxLU[1] == nnIndxLU[2]
        D[1] = 1.0 / α
        l = 2
    else 
        l = 1
    end
    for i in l:n
        if ν == 0.5
            NNdistM = cholesky(exp.(-ϕ .* pairwise(Euclidean(), 
                        coords[nnIndx[nnIndxLU[i]:(nnIndxLU[i + 1] - 1)], :], 
                        dims = 1)) + (1.0 / α - 1.0) * I)
            NNdistv = NNdistM.L \ (exp.(-ϕ .* nnDist[nnIndxLU[i]:(nnIndxLU[i + 1] - 1)]))
            D[i] = 1.0 / α - NNdistv⋅NNdistv
            A[nnIndxLU[i]:(nnIndxLU[i + 1] - 1)] = NNdistM.U \ NNdistv
        end
    end
    
    # return (A = A, D = D)
end


function kfoldperm(N,k)

    # k folder split of 1:N
    
    n,r = divrem(N,k)
    b = collect(1:n:N+1)
    for i in 1:length(b)
        b[i] += i > r ? r : i-1  
    end
    p = randperm(N)
    return [p[r] for r in [b[i]:b[i+1]-1 for i=1:k]]
end

colnorm(A) = [norm(A[:,i], 2) for i=1:size(A,2)]



