#include <iomanip>
#include <string>
#include <limits>
/*#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <R_ext/BLAS.h>
#include <R_ext/Utils.h>*/

///////////////////////////////////////////////////////////////////
//Brute force 
///////////////////////////////////////////////////////////////////
extern "C" {
  void mkNNIndx(int *nPtr, int *mPtr, double *coords, int *nnIndx, double *nnDist, int *nnIndxLU, int *nThreadsPtr);
}

///////////////////////////////////////////////////////////////////
//code book
///////////////////////////////////////////////////////////////////
double dmi(double *x, double *c, int inc);

double dei(double *x, double *c, int inc);

void fastNN(int m, int n, double *coords, int ui, double *u, int *sIndx, int *rSIndx, double *rSNNDist);

extern "C" {
  void mkNNIndxCB(int *nPtr, int *mPtr, double *coords, int *nnIndx, double *nnDist, int *nnIndxLU, int *nThreadsPtr);
}

//Description: given a location's index i and number of neighbors m this function provides the index to i and number of neighbors in nnIndx
void getNNIndx(int i, int m, int &iNNIndx, int &iNN);

double dist2(double a1, double a2, double b1, double b2);
