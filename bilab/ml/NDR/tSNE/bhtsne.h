/**
 * Author: Wei Cao
 * Date:   July 30, 2015
 *
 * Copyright (c) 2015, Wei Cao. All rights reserved.
 */
#ifndef BHTSNE_H
#define BHTSNE_H
#include <math.h>
#include <float.h>
#include <stdlib.h>
#include <stdio.h>
#include <cstring>
#include <time.h>


struct BHTSNE
{

  int N,D; // NxD matrix (N = #points, D = dimensionality)
  int no_dims = 2; // m' x n' after mapping
  double* X = nullptr;
  double* Y = nullptr;

  double perplexity = 30.0;
  double theta = 0.5;
  bool verbosity = false;
  bool scaling = true;
  int randseed = -1;
  
  BHTSNE(double* samples, int numOfsamples, int dims, int mapped_D, double perplex, 
         double th, int rseed, bool scale, bool verbose){
    verbosity = verbose;
    scaling = scale;
    if(rseed >= 0) {
      if(verbosity){
        printf("Using random seed: %d\n", rseed);
      }
      srand((unsigned int) rseed);
    } else {
      if(verbosity){
        printf("Using current time as random seed...\n");
      }
      srand(static_cast<unsigned int>(time(NULL)));
    }
    N = numOfsamples;
    D = dims;
    no_dims = mapped_D;
    perplexity = perplex;
    theta = th;
    randseed = rseed;
    
    Y = (double*) malloc(N * no_dims * sizeof(double));
    if(Y == NULL ) {
      printf("Memory allocation failed!\n"); 
      exit(1); 
    }
    for(int i = 0; i < N * no_dims; i++) {
      Y[i] = randn() * .0001;
    }
    X = samples;
  }
  
  ~BHTSNE(){
    if (Y != NULL)
    {
      free(Y);
    }
  }
 
  void symmetrizeMatrix(unsigned int** _row_P,
                        unsigned int** _col_P,
                        double** _val_P,
                        int N);
  
  void computeGradient(double* P, unsigned int* inp_row_P,
                                  unsigned int* inp_col_P, 
                                  double* inp_val_P,
                                  double *Y,
                                  int N,
                                  int D,
                                  double* dC,
                                  double theta);

  void computeExactGradient(double* P, double* Y, int N, int D, double* dC);
  
  double evaluateError(double* P, double* Y, int N, int D);

  double evaluateError(unsigned int* row_P,
                       unsigned int* col_P,
                       double* val_P,
                       double* Y,
                       int N, int D, double theta);

  void zeroMean(double* X, int N, int D);

  void computeGaussianPerplexity(double* X, int N, int D, double* P, double perplexity);
  
  void computeGaussianPerplexity(double* X, int N, int D,
                                 unsigned int** _row_P,
                                 unsigned int** _col_P,
                                 double** _val_P,
                                 double perplexity,
                                 int K);

  
  void computeSquaredEuclideanDistance(double *X, int N, int D, double* DD);

  double randn();
  void print_();
  void run();
};

#endif
