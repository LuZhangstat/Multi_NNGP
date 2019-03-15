# build nearest neighbor for observed location set S #
function BuildNN(coords, m, brute = 0.0)
    # coords: 2 by n array
    # m:      number of neighbors
    n = size(coords)[2]
    coords_fit = [coords[1, :] coords[2, :]]
    nIndx = trunc(Cint,(1 + m) / 2 * m + (n - m - 1) * m);
    nnIndx = zeros(Cint, nIndx);
    nnDist = zeros(Float64, nIndx);
    ##first of nnIndxLU column holds the nnIndx index for the i-th location and the second columns holds the number of neighbors the i-th location has (the second column is a bit of a waste but simplifies my life in the spNNGP)
    nnIndxLU = zeros(Cint, 2 * n);
    ##Note I need to convert n and m to pointers only because the R example of .C wants only pointers so that's how I worte the mkNNIndx and mkNNIndxTree0
    nPtr = Cint[n];
    mPtr = Cint[m];
    
    if brute == 1.0
        ccall((:mkNNIndx, "julia-R-nn-ccall/nn.so"), Cvoid, (Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble}, 
            Ptr{Cint}, Ptr{Cdouble}, Ptr{Cint}), nPtr, mPtr, coords_fit, nnIndx, nnDist, nnIndxLU)
    else
        ccall((:mkNNIndxTree0, "julia-R-nn-ccall/nn.so"), Cvoid, (Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble}, 
            Ptr{Cint}, Ptr{Cdouble}, Ptr{Cint}), nPtr, mPtr, coords_fit, nnIndx, nnDist, nnIndxLU)
    end
    
    nnIndx = nnIndx .+ 1 
    nnIndxLU = nnIndxLU[2:(n + 1)]  .+ 1
    nnIndxLU[n] = nIndx + 1
    return (nnIndx = nnIndx, nnDist = nnDist, nnIndxLU = nnIndxLU)
end  

function getAD(coords, nnIndx, nnDist, nnIndxLU, ϕ, ν, A, D)
    # coords: 2 by n array,
    # nnIndx, nnDist, nnIndxLU output of the mkNNIndx function
    
    # return: A D
    # A = []
    # D = ones(size(coords)[2])
    
    # initialize memory
    NNdistM = []
    NNdistv = []
    n = size(coords)[2]
    for i in 2:n
        if ν == 0.5
            NNdistM = cholesky(exp.(-ϕ .* pairwise(Euclidean(), 
                        coords[:, nnIndx[nnIndxLU[i - 1]:(nnIndxLU[i] - 1)]], 
                        dims = 2)))
            NNdistv = NNdistM.L \ (exp.(-ϕ .* nnDist[nnIndxLU[i - 1]:(nnIndxLU[i] - 1)]))
            D[i] = 1.0 - NNdistv⋅NNdistv
            A[nnIndxLU[i - 1]:(nnIndxLU[i] - 1)] = NNdistM.U \ NNdistv
        end
    end
    D[1] = 1.0;
    # return (A = A, D = D)
end