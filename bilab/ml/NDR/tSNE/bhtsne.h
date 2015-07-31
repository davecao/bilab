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
  double* cost = nullptr;
  double perplexity = 30.0;
  double theta = 0.5;
  bool verbosity = false;
  int randseed = -1;
  
  BHTSNE(double* samples, int numOfsamples, int mapped_D, double perplex, 
         double th, int rseed, bool verbose){
    verbosity = verbose;
    if(rseed >= 0) {
      printf("Using random seed: %d\n", rseed);
      srand((unsigned int) rseed);
    } else {
      printf("Using current time as random seed...\n");
      srand(time(NULL));
    }
    N = numOfsamples;
    no_dims = mapped_D;
    perplexity = perplex;
    theta = th;
    randseed = rseed;

    int* landmarks = (int*) malloc(N * sizeof(int));
    if(landmarks == NULL) { 
      printf("Memory allocation failed!\n"); 
      exit(1); 
    }
    for(int n = 0; n < N; n++) {
      landmarks[n] = n;
    }
    double* Y = (double*) malloc(N * no_dims * sizeof(double));
    double* costs = (double*) calloc(N, sizeof(double));
    if(Y == NULL || costs == NULL) { 
      printf("Memory allocation failed!\n"); 
      exit(1); 
    }
    for(int i = 0; i < N * no_dims; i++) {
      Y[i] = randn() * .0001;
    }
    X = samples;
    //run();

  }

  void symmetrizeMatrix(unsigned int** _row_P, unsigned int** _col_P, 
                                              double** _val_P, int N);
// void computeGradient(double* P, unsigned int* inp_row_P, 
//                                 unsigned int* inp_col_P, 
//                                 double* inp_val_P, 
//                                 double* Y, 
//                                 int N, int D, double* dC, double theta);

// void computeExactGradient(double* P, double* Y, int N, int D, double* dC);
  void computeGradient(double* P, unsigned int* inp_row_P, 
                                  unsigned int* inp_col_P, 
                                  double* inp_val_P, 
                                  //int N, int D, 
                                  double* dC);//, double theta);

  void computeExactGradient(double* P, 
                            //int N, int D, 
                            double* dC);
//  double evaluateError(double* P, double* Y, int N, int D);

//  double evaluateError(unsigned int* row_P, unsigned int* col_P, 
//                        double* val_P, double* Y, int N, int D, double theta);
  double evaluateError(double* P);

  double evaluateError(unsigned int* row_P, unsigned int* col_P, 
                        double* val_P);
  //void zeroMean(double* X, int N, int D);
  void zeroMean();
  void computeGaussianPerplexity(double* X, int N, int D, 
                                 double* P, double perplexity);
  
  void computeGaussianPerplexity(double* X, int N, int D, 
                                 unsigned int** _row_P, unsigned int** _col_P, 
                                 double** _val_P, double perplexity, int K);
  void computeSquaredEuclideanDistance(double* X, int N, int D, double* DD);

  double randn();

  void run();
};

#endif
