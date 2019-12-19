##rm *.o *.so && g++ -c -fpic nn.cpp -o nn.o && g++ -shared -o nn.so nn.o
rm *.o *.so
g++ -I/home/lu/.julia/packages/CxxWrap/sarOk/deps/usr/include -I/home/lu/julia/julia-1.1.0/include/julia -c -fopenmp -fpic nn.cpp -o nn.o
g++ -shared -o nn.so nn.o


using CSV

coords = convert(Array, CSV.read("coords.csv"));

n = size(coords)[1]
m = 5
nThreads = 2

nIndx = trunc(Cint,(1+m)/2*m+(n-m-1)*m);

nnIndx = zeros(Cint, nIndx, 1);

nnDist = zeros(Float64, nIndx, 1);

##first of nnIndxLU column holds the nnIndx index for the i-th location and the second columns holds the number of neighbors the i-th location has (the second column is a bit of a waste but simplifies my life in the spNNGP)
nnIndxLU = zeros(Cint, 2*n, 1);

##Note I need to convert n and m to pointers only because the R example of .C wants only pointers so that's how I worte the mkNNIndx and mkNNIndxTree0
nPtr = Cint[n];
mPtr = Cint[m];
nThreadsPtr = Cint[nThreads];

#ccall((:mkNNIndx, "./nn.so"), Cvoid, (Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble}, Ptr{Cint}, Ptr{Cdouble}, Ptr{Cint}), nPtr, mPtr, coords, nnIndx, nnDist, nnIndxLU)

using CxxWrap

# Load the module and generate the functions
module Cppnn2
  using CxxWrap
  @wrapmodule(joinpath("./","libnn2"))

#  function __init__()
#    @initcxx
#  end
end

# Call greet and show the result
@show Cppnn2.mkNNIndxCB(nPtr, mPtr, coords, nnIndx, nnDist, nnIndxLU, nThreadsPtr)

ccall((:mkNNIndxCB, "./libnn2.so"), Cvoid, (Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble}, Ptr{Cint}, Ptr{Cdouble}, Ptr{Cint}, Ptr{Cint}), nPtr, mPtr, coords, nnIndx, nnDist, nnIndxLU, nThreadsPtr)


