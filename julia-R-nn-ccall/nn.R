##rm *.o *.so && g++ -c -fpic nn.cpp -o nn.o && g++ -shared -o nn.so nn.o

rm(list=ls())
dyn.load("nn.so")

set.seed(1)

n <- 100
coords <- cbind(runif(n, 0, 1), runif(n, 0, 1))
coords <- coords[order(coords[,1]),]

write.csv(coords, "coords.csv", row.names=FALSE)

n <- nrow(coords)

m <- 5

nIndx <- (1+m)/2*m+(n-m-1)*m

nnIndx <- rep(0, nIndx)

nnDist <- rep(0, nIndx)

##first of nnIndxLU column holds the nnIndx index for the i-th location and the second columns holds the number of neighbors the i-th location has (the second column is a bit of a waste but simplifies my life in the spNNGP).
nnIndxLU <- matrix(0, n, 2)

#out <- .C("mkNNIndx", as.integer(n), as.integer(m), as.double(coords), as.integer(nnIndx), as.double(nnDist), as.integer(nnIndxLU), DUP=FALSE)

out <- .C("mkNNIndxTree0", as.integer(n), as.integer(m), as.double(coords), as.integer(nnIndx), as.double(nnDist), as.integer(nnIndxLU), DUP=FALSE)

##some plots to check
source("utils.R")

n.indx <- mk.n.indx.list(out[[4]], n, m)

for(i in 1:n){
    
    plot(coords, cex=2, xlab="Easting", ylab="Northing")
    abline(v=coords[i,1,drop=FALSE], lty=3, lwd=2)
    points(coords[i,,drop=FALSE], col="blue", pch=19, cex=2)
    points(coords[n.indx[[i]],,drop=FALSE], col="red", pch=19, cex=2)
    
    readline(prompt = "Pause. Press <Enter> to continue...")
}
