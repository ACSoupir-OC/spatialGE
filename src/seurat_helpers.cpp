// CODE TAKEN MODIFIED FROM THE SEURAT PACKAGE
//
// Copyright (c) 2021 Seurat authors
// Permission is hereby granted, free of charge, to any person obtaining a copy of this 
// software and associated documentation files (the "Software"), to deal in the Software 
// without restriction, including without limitation the rights to use, copy, modify, merge, 
// publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons 
// to whom the Software is furnished to do so, subject to the following conditions:
// The above copyright notice and this permission notice shall be included in all copies or 
// substantial portions of the Software.
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, 
// INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR 
// PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE 
// FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR 
// OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER 
// DEALINGS IN THE SOFTWARE.


#include <RcppEigen.h>

using namespace Rcpp;
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppProgress)]]

/* use this if you know the row means */
// [[Rcpp::export(rng = false)]]
NumericVector SparseRowVar2(Eigen::SparseMatrix<double> mat, NumericVector mu){
  mat = mat.transpose();
  NumericVector allVars = no_init(mat.cols());
  for (int k=0; k<mat.outerSize(); ++k){
    double colSum = 0;
    int nZero = mat.rows();
    for (Eigen::SparseMatrix<double>::InnerIterator it(mat,k); it; ++it) {
      nZero -= 1;
      colSum += pow(it.value() - mu[k], 2);
    }
    colSum += pow(mu[k], 2) * nZero;
    allVars[k] = colSum / (mat.rows() - 1);
  }
  return(allVars);
}

/* standardize matrix rows using given mean and standard deviation,
   clip values larger than vmax to vmax,
   then return variance for each row */
// [[Rcpp::export(rng = false)]]
NumericVector SparseRowVarStd(Eigen::SparseMatrix<double> mat, NumericVector mu, NumericVector sd, double vmax){
  mat = mat.transpose();
  NumericVector allVars(mat.cols());
  for (int k=0; k<mat.outerSize(); ++k){
    if (sd[k] == 0) continue;
    double colSum = 0;
    int nZero = mat.rows();
    for (Eigen::SparseMatrix<double>::InnerIterator it(mat,k); it; ++it)
    {
      nZero -= 1;
      colSum += pow(std::min(vmax, (it.value() - mu[k]) / sd[k]), 2);
    }
    colSum += pow((0 - mu[k]) / sd[k], 2) * nZero;
    allVars[k] = colSum / (mat.rows() - 1);
  }
  return(allVars);
}

// [[Rcpp::export(rng = false)]]
NumericVector STenrich_permutation_fast(NumericMatrix coords, int n_subsample, int n_samples){
  int n_rows = coords.nrow();
  NumericVector result(n_samples);
  
  for(int s = 0; s < n_samples; s++){
    // Randomly select n_subsample spots
    IntegerVector idx = sample(n_rows, n_subsample);
    
    // Calculate mean distance between all pairs in subsample
    double total_dist = 0;
    int n_pairs = 0;
    for(int i = 0; i < n_subsample; i++){
      for(int j = i+1; j < n_subsample; j++){
        double dx = coords(idx[i], 0) - coords(idx[j], 0);
        double dy = coords(idx[i], 1) - coords(idx[j], 1);
        total_dist += sqrt(dx*dx + dy*dy);
        n_pairs++;
      }
    }
    
    result[s] = n_pairs > 0 ? total_dist / n_pairs : 0;
  }
  
  return result;
}

// [[Rcpp::export(rng = false)]]
double STenrich_pvalue_calculation(NumericVector null_distribution, double observed_distance, int n_samples){
  int n_greater = 0;
  for(int i = 0; i < null_distribution.length(); i++){
    if(null_distribution[i] >= observed_distance){
      n_greater++;
    }
  }
  
  // Calculate p-value with continuity correction
  double pvalue = (n_greater + 1) / (n_samples + 1);
  return pvalue;
}

// Alias for STenrich_permutation_helper
// [[Rcpp::export(rng = false)]]
NumericVector STenrich_permutation_helper(NumericMatrix coords, int n_subsample, int n_samples){
  return STenrich_permutation_fast(coords, n_subsample, n_samples);
}

