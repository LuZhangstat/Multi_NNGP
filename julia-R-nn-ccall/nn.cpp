#include <limits>
#include <iostream>
#include <cmath>
#include "nn.h"
using namespace std;

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

//trees
Node *miniInsert(Node *Tree, double *coords, int index, int d,int n){

  int P = 2;
  
  if(Tree==NULL){
    return new Node(index);
  }
  
  if(coords[index]<=coords[Tree->index]&&d==0){
    Tree->left=miniInsert(Tree->left,coords,index,(d+1)%P,n);
  }
  
  if(coords[index]>coords[Tree->index]&&d==0){ 
    Tree->right=miniInsert(Tree->right,coords,index,(d+1)%P,n);
  }
  
  if(coords[index+n]<=coords[Tree->index+n]&&d==1){
    Tree->left=miniInsert(Tree->left,coords,index,(d+1)%P,n);
  }
  
  if(coords[index+n]>coords[Tree->index+n]&&d==1){ 
    Tree->right=miniInsert(Tree->right,coords,index,(d+1)%P,n);
  }
  
  return Tree;
}

void get_nn(Node *Tree, int index, int d, double *coords, int n, double *nnDist, int *nnIndx, int iNNIndx, int iNN, int check){

  int P = 2;
  
  if(Tree==NULL){
    return;
  }
  
  double disttemp= dist2(coords[index],coords[index+n],coords[Tree->index],coords[Tree->index+n]); 
  
  if(index!=Tree->index && disttemp<nnDist[iNNIndx+iNN-1]){
    nnDist[iNNIndx+iNN-1]=disttemp;
    nnIndx[iNNIndx+iNN-1]=Tree->index;
    fSort(&nnDist[iNNIndx], &nnIndx[iNNIndx], iNN);
    //rsort_with_index(&nnDist[iNNIndx], &nnIndx[iNNIndx], iNN);
  }
  
  Node *temp1=Tree->left;
  Node *temp2=Tree->right;
  
  if(d==0){
    
    if(coords[index]>coords[Tree->index]){
      std::swap(temp1,temp2);
    }
    
    get_nn(temp1,index,(d+1)%P,coords,n, nnDist, nnIndx, iNNIndx, iNN, check);
    
    if(fabs(coords[Tree->index]-coords[index])>nnDist[iNNIndx+iNN-1]){
      return;
    }
    
    get_nn(temp2,index,(d+1)%P,coords,n, nnDist, nnIndx, iNNIndx, iNN, check);
  }
  
  if(d==1){
    
    if(coords[index+n]>coords[Tree->index+n]){
      std::swap(temp1,temp2);
    }
    
    get_nn(temp1,index,(d+1)%P,coords,n, nnDist, nnIndx, iNNIndx, iNN,check);

    if(fabs(coords[Tree->index+n]-coords[index+n])>nnDist[iNNIndx+iNN-1]){
      return;
    }
    
    get_nn(temp2,index,(d+1)%P,coords,n, nnDist, nnIndx, iNNIndx, iNN,check);
  }

}


void mkNNIndxTree0(int *nPtr, int *mPtr, double *coords, int *nnIndx, double *nnDist, int *nnIndxLU){

  int n = nPtr[0];
  int m = mPtr[0];
  
  int i, iNNIndx, iNN;
  double d;
  int nIndx = static_cast<int>(static_cast<double>(1+m)/2*m+(n-m-1)*m);
  int BUCKETSIZE = 10;

  for(i = 0; i < nIndx; i++){
    nnDist[i] = std::numeric_limits<double>::infinity();
  }
  
  Node *Tree=NULL;
  int time_through=-1;
  
  for(i=0;i<n;i++){
    getNNIndx(i, m, iNNIndx, iNN);
    nnIndxLU[i] = iNNIndx;
    nnIndxLU[n+i] = iNN;
    if(time_through==-1){
      time_through=i;
    }
    
    if(i!=0){
      for(int j = time_through; j < i; j++){ 
	getNNIndx(i, m, iNNIndx, iNN);
	d = dist2(coords[i], coords[i+n], coords[j], coords[n+j]);
	if(d < nnDist[iNNIndx+iNN-1]){
	  nnDist[iNNIndx+iNN-1] = d;
	  nnIndx[iNNIndx+iNN-1] = j;
	  fSort(&nnDist[iNNIndx], &nnIndx[iNNIndx], iNN);
	  //rsort_with_index(&nnDist[iNNIndx], &nnIndx[iNNIndx], iNN);
	}
      }
      
      
      if(i%BUCKETSIZE==0){

	for(int j=time_through;j<time_through+BUCKETSIZE;j++){
	  
	  getNNIndx(j, m, iNNIndx, iNN);
	  get_nn(Tree,j,0, coords,n, nnDist,nnIndx,iNNIndx,iNN,i-BUCKETSIZE);
	}
	
	
	for(int j=time_through;j<time_through+BUCKETSIZE;j++){
	  Tree=miniInsert(Tree,coords,j,0, n);
	}
	
	time_through=-1;
      }
      if(i==n-1){
	
	for(int j=time_through;j<n;j++){
	  getNNIndx(j, m, iNNIndx, iNN);
	  get_nn(Tree,j,0, coords,n, nnDist,nnIndx,iNNIndx,iNN,i-BUCKETSIZE);
	}
	
      }
    } 
    if(i==0){
      Tree=miniInsert(Tree,coords,i,0,n);
      time_through=-1;
    }
  }

  
  delete Tree;
}


double dist2(double a1, double a2, double b1, double b2){
  return(sqrt(pow(a1-b1,2)+pow(a2-b2,2)));
}


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


void mkNNIndx(int *nPtr, int *mPtr, double *coords, int *nnIndx, double *nnDist, int *nnIndxLU){

  int n = nPtr[0];
  int m = mPtr[0];
  
  int i, j, iNNIndx, iNN;
  double d;
  
  int nIndx = static_cast<int>(static_cast<double>(1+m)/2*m+(n-m-1)*m);

  for(i = 0; i < nIndx; i++){
    nnDist[i] = std::numeric_limits<double>::infinity();
  }

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
  
}

