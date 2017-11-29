#define ARMA_64BIT_WORD 1
#define _USE_MATH_DEFINES
#include <omp.h>
#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
#include <cmath>
#include <thread>
using namespace Rcpp;
using namespace arma;
using namespace std;
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]] 
// [[Rcpp::plugins(openmp)]]

void inline parallel(int ncores) {
  // int threads_available = std::thread::hardware_concurrency();
  // if (threads_available != 0) {
  //   Rcout << "Using " << threads_available << " cores" << endl;
  //   omp_set_num_threads(threads_available);
  //   // Rcout << "Using: " << omp_get_num_threads() << "cores" << endl;
  // }
  // else {
  //   Rcout << "Openmp coudnl't load, using 1 core.";
  // }
  Rcout << "Using " << ncores << " cores" << endl;
  omp_set_num_threads(ncores);
}

// [[Rcpp::export]]
arma::mat gaussian_filter(arma::mat x, arma::vec z, arma::rowvec bandwidth, int ncores = 1) {
  // setup
  parallel(ncores);
  uword n = x.n_rows, m = bandwidth.n_elem;
  mat out(n, m, fill::zeros);
  mat total_weight(n, m, fill::zeros);
  
  // constants
  const rowvec b2 = 2 * pow(bandwidth, 2);
  
  #pragma omp parallel for schedule(static) 
  for (uword i = 0; i < n; i++) {
    for (uword j = i; j < n; j++) {
      rowvec w = exp(- sum(square(x.row(i) - x.row(j))) / b2);
      total_weight.row(i) += w;
      total_weight.row(j) += w;
      out.row(i) += w * z(j);
      out.row(j) += w * z(i);
    }
  }
  out = out / total_weight;
  
  return out;
}


// // [[Rcpp::export]]
// void box_filter(arma::mat x, int level, int ncores = 1) {
//   // Load parallelism setup
//   parallel_setup(ncores);
// 
//   //
// 
// }

