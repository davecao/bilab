
#ifndef BHTSNE_IMPL_H
#define BHTSNE_IMPL_H

#include "vptree.h"
#include "sptree.h"
#include "bhtsne.h"

#include <iostream>     // std::cout, std::fixed
#include <iomanip>      // std::setprecision


using namespace std;

static double sign(double x) {
  return (x == .0 ? .0 : (x < .0 ? -1.0 : 1.0));
}
/*
void moment(double** data, int double &ave, double &adev, double &sdev, double &var, double &skew, double &curt)
{
  int j;
  DP ep=0.0,s,p;
  
  int n=data.size();
  if (n <= 1) nrerror("n must be at least 2 in moment");
  s=0.0;
  for (j=0;j<n;j++) s += data[j];
  ave=s/n;
  adev=var=skew=curt=0.0;
  for (j=0;j<n;j++) {
    adev += fabs(s=data[j]-ave);
    ep += s;
    var += (p=s*s);
    skew += (p *= s);
    curt += (p *= s);
  }
  adev /= n;
  var=(var-ep*ep/n)/(n-1);
  sdev=sqrt(var);
  if (var != 0.0) {
    skew /= (n*var*sdev);
    curt=curt/(n*var*var)-3.0;
  } else nrerror("No skew/kurtosis when variance = 0 (in moment)");
}
*/

void mean_std_scaling(double* data, int nrows, int ncols){
  std::vector<double> ave_col(ncols);
  std::vector<double> std_col(ncols);
  
  for (int i=0; i<ncols; i++) {
    double s    = 0.0;
    double adev = 0.0;
    double var  = 0.0;
    double skew = 0.0;
    double curt = 0.0;
    double ep   = 0.0;
    for (int j=0; j<nrows; j++) {
      s += data[j*ncols + i];
    }
    double ave = s/nrows;
    double p = 0.0;
    for (int j=0; j<nrows; j++) {
      adev += fabs(s=data[j*ncols + i]-ave);
      ep += s;
      var += (p=s*s);
      skew += (p *= s);
      curt += (p *= s);
    }
    
//    ave_col.push_back(col_sum/nrows);
  }
}

void BHTSNE::print_(){
  
  cout.setf(ios::fixed, ios::floatfield); // set fixed floating format
  cout.precision(4); // for fixed format, two decimal places
  for (unsigned int i = 0; i != N; ++i)
  {
    for (unsigned int j = 0; j != no_dims; ++j)
    {
      if (j != 0) std::cout << " ";
      //std::cout << std::setw(5) << std::setfill(' ') << std::setprecision(5)
      cout << setw(8) << Y[i*no_dims + j];
    }
    std::cout << "\n";
  }
}


void BHTSNE::symmetrizeMatrix(unsigned int** _row_P, unsigned int** _col_P,
                              double** _val_P)
{
  // Get sparse matrix
  unsigned int* row_P = *_row_P;
  unsigned int* col_P = *_col_P;
  double* val_P = *_val_P;
  
  // Count number of elements and row counts of symmetric matrix
  int* row_counts = (int*) calloc(N, sizeof(int));
  if(row_counts == NULL) {
    printf("Memory allocation failed!\n");
    exit(1);
  }
  for(int n = 0; n < N; n++) {
    for(unsigned int i = row_P[n]; i < row_P[n + 1]; i++) {
      // Check whether element (col_P[i], n) is present
      bool present = false;
      for(unsigned int m = row_P[col_P[i]]; m < row_P[col_P[i] + 1]; m++) {
        if(col_P[m] == (unsigned int)n) {
          present = true;
        }
      }
      if(present) {
        row_counts[n]++;
      } else {
        row_counts[n]++;
        row_counts[col_P[i]]++;
      }
    }
  }
  int no_elem = 0;
  for(int n = 0; n < N; n++) {
    no_elem += row_counts[n];
  }
  // Allocate memory for symmetrized matrix
  unsigned int* sym_row_P = (unsigned int*) malloc((N + 1) * sizeof(unsigned int));
  unsigned int* sym_col_P = (unsigned int*) malloc(no_elem * sizeof(unsigned int));
  double* sym_val_P = (double*) malloc(no_elem * sizeof(double));
  if(sym_row_P == NULL || sym_col_P == NULL || sym_val_P == NULL) {
    printf("Memory allocation failed!\n");
    exit(1);
  }
  // Construct new row indices for symmetric matrix
  sym_row_P[0] = 0;
  for(int n = 0; n < N; n++) {
    sym_row_P[n + 1] = sym_row_P[n] + (unsigned int) row_counts[n];
  }
  // Fill the result matrix
  int* offset = (int*) calloc(N, sizeof(int));
  if(offset == NULL) {
    printf("Memory allocation failed!\n");
    exit(1);
  }
  for(int n = 0; n < N; n++) {
    for(unsigned int i = row_P[n]; i < row_P[n + 1]; i++) {
      // Check whether element (col_P[i], n) is present
      bool present = false;
      for(unsigned int m = row_P[col_P[i]]; m < row_P[col_P[i] + 1]; m++) {
        if(col_P[m] == (unsigned int)n) {
          present = true;
          if((unsigned int)n <= col_P[i]) {// make sure we do not add elements twice
            sym_col_P[sym_row_P[n]        + offset[n]]        = col_P[i];
            sym_col_P[sym_row_P[col_P[i]] + offset[col_P[i]]] = n;
            sym_val_P[sym_row_P[n]        + offset[n]]        = val_P[i] + val_P[m];
            sym_val_P[sym_row_P[col_P[i]] + offset[col_P[i]]] = val_P[i] + val_P[m];
          }
        }
      }
      // If (col_P[i], n) is not present, there is no addition involved
      if(!present) {
        sym_col_P[sym_row_P[n]        + offset[n]]        = col_P[i];
        sym_col_P[sym_row_P[col_P[i]] + offset[col_P[i]]] = n;
        sym_val_P[sym_row_P[n]        + offset[n]]        = val_P[i];
        sym_val_P[sym_row_P[col_P[i]] + offset[col_P[i]]] = val_P[i];
      }
      // Update offsets
      if(!present || (present && (unsigned int)n <= col_P[i])) {
        offset[n]++;
        if(col_P[i] != (unsigned int)n) {
          offset[col_P[i]]++;
        }
      }
    }
  }
  
  // Divide the result by two
  for(int i = 0; i < no_elem; i++) {
    sym_val_P[i] /= 2.0;
  }
  // Return symmetrized matrices
  free(*_row_P); *_row_P = sym_row_P;
  free(*_col_P); *_col_P = sym_col_P;
  free(*_val_P); *_val_P = sym_val_P;
  // Free up some memery
  free(offset); offset = NULL;
  free(row_counts); row_counts  = NULL;
}

// Compute gradient of the t-SNE cost function (using Barnes-Hut algorithm)
void BHTSNE::computeGradient(double* P, unsigned int* inp_row_P,
                             unsigned int* inp_col_P,
                             double* inp_val_P,
                             double* dC)
{
  // Construct space-partitioning tree on current map
  SPTree* tree = new SPTree(no_dims, Y, N);
  
  // Compute all terms required for t-SNE gradient
  double sum_Q = .0;
  double* pos_f = (double*) calloc(N * no_dims, sizeof(double));
  double* neg_f = (double*) calloc(N * no_dims, sizeof(double));
  if(pos_f == NULL || neg_f == NULL) {
    printf("Memory allocation failed!\n");
    exit(1);
  }
  tree->computeEdgeForces(inp_row_P, inp_col_P, inp_val_P, N, pos_f);
  for(int n = 0; n < N; n++) {
    tree->computeNonEdgeForces(n, theta, neg_f + n * no_dims, &sum_Q);
  }
  // Compute final t-SNE gradient
  for(int i = 0; i < N * no_dims; i++) {
    dC[i] = pos_f[i] - (neg_f[i] / sum_Q);
  }
  free(pos_f);
  free(neg_f);
  delete tree;
}
// Compute gradient of the t-SNE cost function (exact)
void BHTSNE::computeExactGradient(double* P,
                                  double* dC)
{
  // Make sure the current gradient contains zeros
  for(int i = 0; i < N * no_dims; i++) {
    dC[i] = 0.0;
  }
  // Compute the squared Euclidean distance matrix
  double* DD = (double*) malloc(N * N * sizeof(double));
  if(DD == NULL) {
    printf("Memory allocation failed!\n");
    exit(1);
  }
  computeSquaredEuclideanDistance(DD);
  
  // Compute Q-matrix and normalization sum
  double* Q = (double*) malloc(N * N * sizeof(double));
  if(Q == NULL) {
    printf("Memory allocation failed!\n");
    exit(1);
  }
  double sum_Q = .0;
  int nN = 0;
  for(int n = 0; n < N; n++) {
    for(int m = 0; m < N; m++) {
      if(n != m) {
        Q[nN + m] = 1 / (1 + DD[nN + m]);
        sum_Q += Q[nN + m];
      }
    }
    nN += N;
  }
  
  // Perform the computation of the gradient
  nN = 0;
  int nD = 0;
  for(int n = 0; n < N; n++) {
    int mD = 0;
    for(int m = 0; m < N; m++) {
      if(n != m) {
        double mult = (P[nN + m] - (Q[nN + m] / sum_Q)) * Q[nN + m];
        for(int d = 0; d < D; d++) {
          dC[nD + d] += (Y[nD + d] - Y[mD + d]) * mult;
        }
      }
      mD += D;
    }
    nN += N;
    nD += D;
  }
  
  // Free memory
  free(DD); DD = NULL;
  free(Q);  Q  = NULL;
}

//double BHTSNE::evaluateError(double* P, int N, int D)
double BHTSNE::evaluateError(double* P)
{
  // Compute the squared Euclidean distance matrix
  double* DD = (double*) malloc(N * N * sizeof(double));
  double* Q = (double*) malloc(N * N * sizeof(double));
  if(DD == NULL || Q == NULL) {
    printf("Memory allocation failed!\n");
    exit(1);
  }
  computeSquaredEuclideanDistance(DD);
  
  // Compute Q-matrix and normalization sum
  int nN = 0;
  double sum_Q = DBL_MIN;
  for(int n = 0; n < N; n++) {
    for(int m = 0; m < N; m++) {
      if(n != m) {
        Q[nN + m] = 1 / (1 + DD[nN + m]);
        sum_Q += Q[nN + m];
      } else {
        Q[nN + m] = DBL_MIN;
      }
    }
    nN += N;
  }
  for(int i = 0; i < N * N; i++) {
    Q[i] /= sum_Q;
  }
  // Sum t-SNE error
  double C = .0;
  for(int n = 0; n < N * N; n++) {
    C += P[n] * log((P[n] + FLT_MIN) / (Q[n] + FLT_MIN));
  }
  
  // Clean up memory
  free(DD);
  free(Q);
  return C;
}

double BHTSNE::evaluateError(unsigned int* row_P, unsigned int* col_P,
                             double* val_P)
{
  // Get estimate of normalization term
  SPTree* tree = new SPTree(no_dims, Y, N);
  double* buff = (double*) calloc(no_dims, sizeof(double));
  double sum_Q = .0;
  for(int n = 0; n < N; n++) {
    tree->computeNonEdgeForces(n, theta, buff, &sum_Q);
  }
  // Loop over all edges to compute t-SNE error
  int ind1, ind2;
  double C = .0, Q;
  for(int n = 0; n < N; n++) {
    ind1 = n * D;
    for(unsigned int i = row_P[n]; i < row_P[n + 1]; i++) {
      Q = .0;
      ind2 = col_P[i] * no_dims;
      for(int d = 0; d < no_dims; d++) buff[d]  = Y[ind1 + d];
      for(int d = 0; d < no_dims; d++) buff[d] -= Y[ind2 + d];
      for(int d = 0; d < no_dims; d++) Q += buff[d] * buff[d];
      Q = (1.0 / (1.0 + Q)) / sum_Q;
      C += val_P[i] * log((val_P[i] + FLT_MIN) / (Q + FLT_MIN));
    }
  }
  // Clean up memory
  free(buff);
  delete tree;
  return C;
}

void BHTSNE::zeroMean(double* X, int N, int D)
{
  // Compute data mean
  double* mean = (double*) calloc(D, sizeof(double));
  if(mean == NULL) {
    printf("Memory allocation failed!\n");
    exit(1);
  }
  int nD = 0;
  for(int n = 0; n < N; n++) {
    for(int d = 0; d < D; d++) {
      mean[d] += X[nD + d];
    }
    nD += D;
  }
  for(int d = 0; d < D; d++) {
    mean[d] /= (double) N;
  }
  // Subtract data mean
  nD = 0;
  for(int n = 0; n < N; n++) {
    for(int d = 0; d < D; d++) {
      X[nD + d] -= mean[d];
    }
    nD += D;
  }
  free(mean);
  mean = NULL;
}

// Compute input similarities with a fixed perplexity
void BHTSNE::computeGaussianPerplexity(double* P)
{
  // Compute the squared Euclidean distance matrix
  double* DD = (double*) malloc(N * N * sizeof(double));
  if(DD == NULL) {
    printf("Memory allocation failed!\n");
    exit(1);
  }
  computeSquaredEuclideanDistance(DD);
  
  // Compute the Gaussian kernel row by row
  int nN = 0;
  for(int n = 0; n < N; n++) {
    // Initialize some variables
    bool found = false;
    double beta = 1.0;
    double min_beta = -DBL_MAX;
    double max_beta =  DBL_MAX;
    double tol = 1e-5;
    double sum_P = 0.0;
    
    // Iterate until we found a good perplexity
    int iter = 0;
    while(!found && iter < 200) {
      // Compute Gaussian kernel row
      for(int m = 0; m < N; m++) {
        P[nN + m] = exp(-beta * DD[nN + m]);
      }
      P[nN + n] = DBL_MIN;
      
      // Compute entropy of current row
      sum_P = DBL_MIN;
      for(int m = 0; m < N; m++) {
        sum_P += P[nN + m];
      }
      double H = 0.0;
      for(int m = 0; m < N; m++) {
        H += beta * (DD[nN + m] * P[nN + m]);
      }
      H = (H / sum_P) + log(sum_P);
      
      // Evaluate whether the entropy is within the tolerance level
      double Hdiff = H - log(perplexity);
      if(Hdiff < tol && -Hdiff < tol) {
        found = true;
      } else {
        if(Hdiff > 0) {
          min_beta = beta;
          if(max_beta == DBL_MAX || max_beta == -DBL_MAX){
            beta *= 2.0;
          } else {
            beta = (beta + max_beta) / 2.0;
          }
        } else {
          max_beta = beta;
          if(min_beta == -DBL_MAX || min_beta == DBL_MAX){
            beta /= 2.0;
          } else{
            beta = (beta + min_beta) / 2.0;
          }
        }
      }
      // Update iteration counter
      iter++;
    }
    
    // Row normalize P
    for(int m = 0; m < N; m++) {
      P[nN + m] /= sum_P;
    }
    nN += N;
  }
  
  // Clean up memory
  free(DD);
  DD = NULL;
}

// Compute input similarities with a fixed perplexity using ball trees (this function allocates memory another function should free)
void BHTSNE::computeGaussianPerplexity(
                                       unsigned int** _row_P, unsigned int** _col_P,
                                       double** _val_P, int K)
{
  if(perplexity > K) {
    printf("Perplexity should be lower than K!\n");
  }
  // Allocate the memory we need
  *_row_P = (unsigned int*) malloc((N + 1) * sizeof(unsigned int));
  *_col_P = (unsigned int*) calloc(N * K, sizeof(unsigned int));
  *_val_P = (double*) calloc(N * K, sizeof(double));
  if(*_row_P == NULL || *_col_P == NULL || *_val_P == NULL) {
    printf("Memory allocation failed!\n");
    exit(1);
  }
  unsigned int* row_P = *_row_P;
  unsigned int* col_P = *_col_P;
  double* val_P = *_val_P;
  double* cur_P = (double*) malloc((N - 1) * sizeof(double));
  if(cur_P == NULL) {
    printf("Memory allocation failed!\n");
    exit(1);
  }
  row_P[0] = 0;
  for(int n = 0; n < N; n++) {
    row_P[n + 1] = row_P[n] + (unsigned int) K;
  }
  // Build ball tree on data set
  VpTree<DataPoint, euclidean_distance>* tree =
  new VpTree<DataPoint, euclidean_distance>();
  vector<DataPoint> obj_X(N, DataPoint(D, -1, X));
  for(int n = 0; n < N; n++) {
    obj_X[n] = DataPoint(D, n, X + n * D);
  }
  tree->create(obj_X);
  
  // Loop over all points to find nearest neighbors
  if (verbosity) {
    printf("Building tree...\n");
  }
  
  vector<DataPoint> indices;
  vector<double> distances;
  for(int n = 0; n < N; n++) {
    if(n % 10000 == 0 && verbosity) {
      printf(" - point %d of %d\n", n, N);
    }
    // Find nearest neighbors
    indices.clear();
    distances.clear();
    tree->search(obj_X[n], K + 1, &indices, &distances);
    
    // Initialize some variables for binary search
    bool found = false;
    double beta = 1.0;
    double min_beta = -DBL_MAX;
    double max_beta =  DBL_MAX;
    double tol = 1e-5;
    
    // Iterate until we found a good perplexity
    int iter = 0;
    double sum_P = 0.0;
    while(!found && iter < 200) {
      // Compute Gaussian kernel row
      for(int m = 0; m < K; m++) {
        cur_P[m] = exp(-beta * distances[m + 1]);
      }
      // Compute entropy of current row
      sum_P = DBL_MIN;
      for(int m = 0; m < K; m++) {
        sum_P += cur_P[m];
      }
      double H = .0;
      for(int m = 0; m < K; m++) {
        H += beta * (distances[m + 1] * cur_P[m]);
      }
      H = (H / sum_P) + log(sum_P);
      // Evaluate whether the entropy is within the tolerance level
      double Hdiff = H - log(perplexity);
      if(Hdiff < tol && -Hdiff < tol) {
        found = true;
      } else {
        if(Hdiff > 0) {
          min_beta = beta;
          if(max_beta == DBL_MAX || max_beta == -DBL_MAX){
            beta *= 2.0;
          } else {
            beta = (beta + max_beta) / 2.0;
          }
        } else {
          max_beta = beta;
          if(min_beta == -DBL_MAX || min_beta == DBL_MAX){
            beta /= 2.0;
          } else {
            beta = (beta + min_beta) / 2.0;
          }
        }
      }
      // Update iteration counter
      iter++;
    }
    
    // Row-normalize current row of P and store in matrix
    for(int m = 0; m < K; m++) {
      cur_P[m] /= sum_P;
    }
    for(int m = 0; m < K; m++) {
      col_P[row_P[n] + m] = (unsigned int) indices[m + 1].index();
      val_P[row_P[n] + m] = cur_P[m];
    }
  }
  if (verbosity) {
    printf("computeGaussianPerplexity using ball tree...finished.\n");
  }
  
  // Clean up memory
  obj_X.clear();
  free(cur_P);
  delete tree;
}
void BHTSNE::computeSquaredEuclideanDistance(double* DD)
{
  double* dataSums = (double*) calloc(N, sizeof(double));
  if(dataSums == NULL) {
    printf("Memory allocation failed!\n");
    exit(1);
  }
  int nD = 0;
  for(int n = 0; n < N; n++) {
    for(int d = 0; d < D; d++) {
      dataSums[n] += (X[nD + d] * X[nD + d]);
    }
    nD += D;
  }
  int nN = 0;
  for(int n = 0; n < N; n++) {
    for(int m = 0; m < N; m++) {
      DD[nN + m] = dataSums[n] + dataSums[m];
    }
    nN += N;
  }
  nN = 0;
  nD = 0;
  for(int n = 0; n < N; n++) {
    int mD = 0;
    DD[nN + n] = 0.0;
    for(int m = n + 1; m < N; m++) {
      DD[nN + m] = 0.0;
      for(int d = 0; d < D; d++) {
        DD[nN + m] += (X[nD + d] - X[mD + d]) * (X[nD + d] - X[mD + d]);
      }
      DD[m * N + n] = DD[nN + m];
      mD += D;
    }
    nN += N; nD += D;
  }
  free(dataSums);
  dataSums = NULL;
}

double BHTSNE::randn()
{
  double x, y, radius;
  do {
    x = 2 * (rand() / ((double) RAND_MAX + 1)) - 1;
    y = 2 * (rand() / ((double) RAND_MAX + 1)) - 1;
    radius = (x * x) + (y * y);
  } while((radius >= 1.0) || (radius == 0.0));
  radius = sqrt(-2 * log(radius) / radius);
  x *= radius;
  y *= radius;
  return x;
}
void BHTSNE::run()
{
  // Determine whether we are using an exact algorithm
  if(N - 1 < 3 * perplexity) {
    printf("Perplexity too large for the number of data points!\n");
    exit(1);
  }
  if (verbosity) {
    printf("Using no_dims = %d, perplexity = %f, and theta = %f\n", no_dims, perplexity, theta);
  }
  
  bool exact = (theta == .0) ? true : false;
  
  // Set learning parameters
  float total_time = .0;
  clock_t start, end;
  int max_iter = 1000, stop_lying_iter = 250, mom_switch_iter = 250;
  double momentum = .5, final_momentum = .8;
  double eta = 200.0;
  
  // Allocate some memory
  double* dY    = (double*) malloc(N * no_dims * sizeof(double));
  double* uY    = (double*) malloc(N * no_dims * sizeof(double));
  double* gains = (double*) malloc(N * no_dims * sizeof(double));
  if(dY == NULL || uY == NULL || gains == NULL) {
    printf("Memory allocation failed!\n");
    exit(1);
  }
  for(int i = 0; i < N * no_dims; i++){
    uY[i] =  .0;
  }
  for(int i = 0; i < N * no_dims; i++) {
    gains[i] = 1.0;
  }
  // Normalize input data (to prevent numerical problems)
  if (verbosity) {
    printf("Computing input similarities...\n");
  }
  
  start = clock();
  zeroMean(X, N, D);
  //  zeroMean();
  double max_X = .0;
  for(int i = 0; i < N * D; i++) {
    if(X[i] > max_X) max_X = X[i];
  }
  for(int i = 0; i < N * D; i++) {
    X[i] /= max_X;
  }
  // Compute input similarities for exact t-SNE
  double* P = nullptr;
  unsigned int* row_P = nullptr;
  unsigned int* col_P = nullptr;
  double* val_P = nullptr;
  if(exact) {
    // Compute similarities
    if(verbosity){
      printf("Exact?");
    }
    P = (double*) malloc(N * N * sizeof(double));
    if(P == NULL) {
      printf("Memory allocation failed!\n");
      exit(1);
    }
    computeGaussianPerplexity(P);
    
    // Symmetrize input similarities
    if (verbosity) {
      printf("Symmetrizing...\n");
    }
    int nN = 0;
    for(int n = 0; n < N; n++) {
      int mN = 0;
      for(int m = n + 1; m < N; m++) {
        P[nN + m] += P[mN + n];
        P[mN + n]  = P[nN + m];
        mN += N;
      }
      nN += N;
    }
    double sum_P = .0;
    for(int i = 0; i < N * N; i++) {
      sum_P += P[i];
    }
    for(int i = 0; i < N * N; i++) {
      P[i] /= sum_P;
    }
  } else { // Compute input similarities for approximate t-SNE
    
    // Compute asymmetric pairwise input similarities
    computeGaussianPerplexity(&row_P, &col_P, &val_P, (int) (3 * perplexity));
    
    // Symmetrize input similarities
    //symmetrizeMatrix(&row_P, &col_P, &val_P, N);
    symmetrizeMatrix(&row_P, &col_P, &val_P);
    double sum_P = .0;
    for(unsigned int i = 0; i < row_P[N]; i++) {
      sum_P += val_P[i];
    }
    for(unsigned int i = 0; i < row_P[N]; i++) {
      val_P[i] /= sum_P;
    }
  }
  end = clock();
  
  // Lie about the P-values
  if(exact) {
    for(int i = 0; i < N * N; i++) {
      P[i] *= 12.0;
    }
  } else {
    for(unsigned int i = 0; i < row_P[N]; i++) {
      val_P[i] *= 12.0;
    }
  }
  
  // Perform main training loop
  if(exact) {
    if (verbosity) {
      printf("Input similarities computed in %4.2f seconds!\nLearning embedding...\n",
             (float) (end - start) / CLOCKS_PER_SEC);
    }
  } else {
    if (verbosity) {
      printf("Input similarities computed in %4.2f seconds (sparsity = %f)!\nLearning embedding...\n",
             (float) (end - start) / CLOCKS_PER_SEC,
             (double) row_P[N] / ((double) N * (double) N));
    }
  }
  start = clock();
  
  for(int iter = 0; iter < max_iter; iter++) {
    
    // Compute (approximate) gradient
    if(exact) {
      //computeExactGradient(P, Y, N, no_dims, dY);
      computeExactGradient(P, dY);
    } else {
      //computeGradient(P, row_P, col_P, val_P, Y, N, no_dims, dY, theta);
      computeGradient(P, row_P, col_P, val_P, dY);
    }
    
    // Update gains
    for(int i = 0; i < N * no_dims; i++) {
      gains[i] = (sign(dY[i]) != sign(uY[i])) ? (gains[i] + .2) : (gains[i] * .8);
    }
    for(int i = 0; i < N * no_dims; i++) {
      if(gains[i] < .01) {
        gains[i] = .01;
      }
    }
    
    // Perform gradient update (with momentum and gains)
    for(int i = 0; i < N * no_dims; i++) {
      uY[i] = momentum * uY[i] - eta * gains[i] * dY[i];
    }
    
    for(int i = 0; i < N * no_dims; i++) {
      Y[i] = Y[i] + uY[i];
    }
    // Make solution zero-mean
    zeroMean(Y, N, no_dims);
    
    // Stop lying about the P-values after a while, and switch momentum
    if(iter == stop_lying_iter) {
      if(exact) {
        for(int i = 0; i < N * N; i++) {
          P[i] /= 12.0;
        }
      } else {
        for(unsigned int i = 0; i < row_P[N]; i++) {
          val_P[i] /= 12.0;
        }
      }
    }
    
    if(iter == mom_switch_iter) {
      momentum = final_momentum;
    }
    // Print out progress
    if(iter > 0 && (iter % 50 == 0 || iter == max_iter - 1)) {
      end = clock();
      double C = .0;
      
      if(exact) {
        //C = evaluateError(P, Y, N, no_dims);
        C = evaluateError(P);
      } else {// doing approximate computation here!
        //C = evaluateError(row_P, col_P, val_P, Y, N, no_dims, theta);
        C = evaluateError(row_P, col_P, val_P);
      }
      if (verbosity) {
        if(iter == 0){
          printf("Iteration %d: error is %f\n", iter + 1, C);
        } else {
          total_time += (float) (end - start) / CLOCKS_PER_SEC;
          printf("Iteration %d: error is %f (50 iterations in %4.2f seconds)\n",
                 iter, C, (float) (end - start) / CLOCKS_PER_SEC);
        }
      }
      start = clock();
    }
  }
  end = clock(); total_time += (float) (end - start) / CLOCKS_PER_SEC;
  
  // Clean up memory
  free(dY);
  free(uY);
  free(gains);
  if(exact) {
    free(P);
  } else {
    free(row_P); row_P = NULL;
    free(col_P); col_P = NULL;
    free(val_P); val_P = NULL;
  }
  if (verbosity) {
    printf("Fitting performed in %4.2f seconds.\n", total_time);
  }
  
}


#endif
