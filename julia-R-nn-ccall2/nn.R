##rm *.o *.so && g++ -c -fpic -fopenmp nn.cpp -o nn.o && g++ -shared -o nn.so nn.o
setwd("/home/lu/Documents/Github/Multi_NNGP/julia-R-nn-ccall2")
rm(list=ls())
dyn.load("nn.so")
#load("../RDA/data/raw_data.RData")

#coords <- cbind(raw_data$scaled_x[!is.na(raw_data$NDVI)], 
#                raw_data$scaled_y[!is.na(raw_data$NDVI)])

set.seed(1)

n <- 100
coords <- cbind(runif(n, 0, 1), runif(n, 0, 1))
write.csv(coords, "coords.csv", row.names=FALSE)

#n <- nrow(coords)

m <- 10

nThreads <- 4

nIndx <- (1+m)/2*m+(n-m-1)*m

nnIndx <- rep(0, nIndx)

nnDist <- rep(0, nIndx)

##first of nnIndxLU column holds the nnIndx index for the i-th location and the second columns holds the number of neighbors the i-th location has (the second column is a bit of a waste but simplifies my life in the spNNGP).
nnIndxLU <- matrix(0, n, 2)

#out <- .C("mkNNIndx", as.integer(n), as.integer(m), as.double(coords), as.integer(nnIndx), as.double(nnDist), as.integer(nnIndxLU), DUP=FALSE)

t <- proc.time()
outCB <- .C("mkNNIndxCB", as.integer(n), as.integer(m), as.double(coords), 
          as.integer(nnIndx), as.double(nnDist), as.integer(nnIndxLU), 
          as.integer(nThreads), DUP=FALSE)
proc.time() - t
##some plots to check
source("utils.R")

n.indx <- mk.n.indx.list(outCB[[4]], n, m)

for(i in 30:35){
    
    plot(coords[1:i, ], cex=2, xlab="Easting", ylab="Northing")
    abline(v=coords[i,1,drop=FALSE], lty=3, lwd=2)
    points(coords[i,,drop=FALSE], col="blue", pch=19, cex=2)
    points(coords[n.indx[[i]],,drop=FALSE], col="red", pch=19, cex=2)
    
    readline(prompt = "Pause. Press <Enter> to continue...")
}


save(outCB, file = "../RDA/data/outCBtest.RData")
