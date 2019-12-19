#include <iomanip>
#include <string>
#include <limits>
#include <iostream>
#include <cmath>
/*#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <R_ext/BLAS.h>
#include <R_ext/Utils.h>*/

#ifdef _OPENMP
#include <omp.h>
#endif

//#include "util.h"
#include "nn.h"


// replace rsort_with_index by fSort
void fSort(double *a, int *b, int n){
  
  int j, k, l;
  double v;
  
  for(j = 1; j <= n-1; j++){
    k = j;  
    while(k > 0 && a[k] < a[k-1]) {
      v = a[k]; l = b[k];
      a[k] = a[k-1]; b[k] = b[k-1];
      a[k-1] = v; b[k-1] = l;
      k--;
    }
  }
}

///////////////////////////////////////////////////////////////////
//Brute force 
///////////////////////////////////////////////////////////////////

void mkNNIndx(int *nPtr, int *mPtr, double *coords, int *nnIndx, double *nnDist, 
              int *nnIndxLU, int *nThreadsPtr){
  
  int i, j, iNNIndx, iNN;
  double d;
  
  int n = nPtr[0];
  int m = mPtr[0];
  int nThreads = nThreadsPtr[0];

  int BUCKETSIZE = 10;

#ifdef _OPENMP
  omp_set_num_threads(nThreads);
#else
  if(nThreads > 1){
    printf("n.omp.threads > %i, but source not compiled with OpenMP support.", nThreads);
    nThreads = 1;
  }
#endif

  int nIndx = static_cast<int>(static_cast<double>(1+m)/2*m+(n-m-1)*m);
  
  for(i = 0; i < nIndx; i++){
    nnDist[i] = std::numeric_limits<double>::infinity();
  }
  
#ifdef _OPENMP
#pragma omp parallel for private(j, iNNIndx, iNN, d)
#endif
  for(i = 0; i < n; i++){
    getNNIndx(i, m, iNNIndx, iNN);
    nnIndxLU[i] = iNNIndx;
    nnIndxLU[n+i] = iNN;   
    if(i != 0){  
      for(j = 0; j < i; j++){	
	d = dist2(coords[i], coords[n+i], coords[j], coords[n+j]);	
	if(d < nnDist[iNNIndx+iNN-1]){	  
	  nnDist[iNNIndx+iNN-1] = d;
	  nnIndx[iNNIndx+iNN-1] = j;
	  fSort(&nnDist[iNNIndx], &nnIndx[iNNIndx], iNN);
	  //rsort_with_index(&nnDist[iNNIndx], &nnIndx[iNNIndx], iNN);
	}	
      }
    }
  }
  return;
}


///////////////////////////////////////////////////////////////////
//code book
///////////////////////////////////////////////////////////////////

//Description: using the fast mean-distance-ordered nn search by Ra and Kim 1993
//Input:
//ui = is the index for which we need the m nearest neighbors
//m = number of nearest neighbors
//n = number of observations, i.e., length of u
//sIndx = the NNGP ordering index of length n that is pre-sorted by u
//u = x+y vector of coordinates assumed sorted on input
//rSIndx = vector or pointer to a vector to store the resulting nn sIndx (this is at most length m for ui >= m)
//rNNDist = vector or point to a vector to store the resulting nn Euclidean distance (this is at most length m for ui >= m)  

double dmi(double *x, double *c, int inc){
    return pow(x[0]+x[inc]-c[0]-c[inc], 2);
}

double dei(double *x, double *c, int inc){
  return pow(x[0]-c[0],2) + pow(x[inc]-c[inc],2);
}

void fastNN(int m, int n, double *coords, int ui, double *u, int *sIndx, 
            int *rSIndx, double *rSNNDist){
  
  int i,j,k;
  bool up, down;
  double dm, de;
  
  //rSNNDist will hold de (i.e., squared Euclidean distance) initially.
  for(i = 0; i < m; i++){
    rSNNDist[i] = std::numeric_limits<double>::infinity();
  }
  
  i = j = ui;
  
  up = down = true;
  
  while(up || down){
    
    if(i == 0){
      down = false;
    }

    if(j == (n-1)){
      up = false;
    }

    if(down){
      
      i--;
      
      dm = dmi(&coords[sIndx[ui]], &coords[sIndx[i]], n);
      
      if(dm > 2*rSNNDist[m-1]){
	down = false;
	
      }else{
	de = dei(&coords[sIndx[ui]], &coords[sIndx[i]], n);

	if(de < rSNNDist[m-1] && sIndx[i] < sIndx[ui]){
	  rSNNDist[m-1] = de;
	  rSIndx[m-1] = sIndx[i];
	  fSort(rSNNDist, rSIndx, m);
	  //rsort_with_index(rSNNDist, rSIndx, m);
	}
	
      }
    }//end down
    
    if(up){
      
      j++;
      
      dm = dmi(&coords[sIndx[ui]], &coords[sIndx[j]], n);
      
      if(dm > 2*rSNNDist[m-1]){
	up = false;
	
      }else{
	de = dei(&coords[sIndx[ui]], &coords[sIndx[j]], n);

	if(de < rSNNDist[m-1] && sIndx[j] < sIndx[ui]){
	  rSNNDist[m-1] = de;
	  rSIndx[m-1] = sIndx[j];
	  fSort(rSNNDist, rSIndx, m);
	  //rsort_with_index(rSNNDist, rSIndx, m);
	}
	
      }
      
    }//end up
    
  }
  
  for(i = 0; i < m; i++){
    rSNNDist[i] = sqrt(rSNNDist[i]);
  }

  return;
}

extern "C" {
  void mkNNIndxCB(int *nPtr, int *mPtr, double *coords, int *nnIndx, double *nnDist, 
                  int *nnIndxLU, int *nThreadsPtr){
    
    double d;
    
    int n = nPtr[0];
    int m = mPtr[0];
    int nThreads = nThreadsPtr[0];
    
    int BUCKETSIZE = 10; 

    
#ifdef _OPENMP
    omp_set_num_threads(nThreads);
#else
    if(nThreads > 1){
      printf("n.omp.threads > %i, but source not compiled with OpenMP support.", nThreads);
      nThreads = 1;
    }
#endif
    
    int i, iNNIndx, iNN;
    
    int *sIndx = new int[n];
    double *u = new double[n];
    
    for(i = 0; i < n; i++){
      sIndx[i] = i;
      u[i] = coords[i]+coords[n+i];
    }
    
    fSort(u, sIndx, n);
    //rsort_with_index(u, sIndx, n); 
    
    //make nnIndxLU and fill nnIndx and d
#ifdef _OPENMP
#pragma omp parallel for private(iNNIndx, iNN)
#endif  
    for(i = 0; i < n; i++){ //note this i indexes the u vector
      getNNIndx(sIndx[i], m, iNNIndx, iNN);
      nnIndxLU[sIndx[i]] = iNNIndx;
      nnIndxLU[n+sIndx[i]] = iNN;   
      fastNN(iNN, n, coords, i, u, sIndx, &nnIndx[iNNIndx], &nnDist[iNNIndx]);
    } 
    
    return;
  }
}

// functions from utils
void getNNIndx(int i, int m, int &iNNIndx, int &iNN){
  
  if(i == 0){
    iNNIndx = 0;//this should never be accessed
    iNN = 0;
    return;
  }else if(i < m){
    iNNIndx = static_cast<int>(static_cast<double>(i)/2*(i-1));
    iNN = i;
    return;
  }else{
    iNNIndx = static_cast<int>(static_cast<double>(m)/2*(m-1)+(i-m)*m);
    iNN = m;
    return;
  } 
}

double dist2(double a1, double a2, double b1, double b2){
  return(sqrt(pow(a1-b1,2)+pow(a2-b2,2)));
}


