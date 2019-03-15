#include <limits>
#include <iostream>
#include <cmath>

//trees
struct Node{
	int index; // which point I am
	Node *left;
	Node *right; 
	Node (int i) { index = i; left = right = NULL; }
};

Node *miniInsert(Node *Tree, double *coords, int index, int d,int n);

void get_nn(Node *Tree, int index, int d, double *coords, int n, double *nnDist, int *nnIndx, int iNNIndx, int iNN, int check);

extern "C" {
  void mkNNIndxTree0(int *nPtr, int *mPtr, double *coords, int *nnIndx, double *nnDist, int *nnIndxLU);
}

double dist2(double a1, double a2, double b1, double b2);

//Description: given a location's index i and number of neighbors m this function provides the index to i and number of neighbors in nnIndx
void getNNIndx(int i, int m, int &iNNIndx, int &iNN);

extern "C" {
  void mkNNIndx(int *nPtr, int *mPtr, double *coords, int *nnIndx, double *nnDist, int *nnIndxLU);
}
