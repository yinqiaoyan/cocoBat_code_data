#include <RcppArmadillo.h>
#include <RcppDist.h>
#include <stdio.h>
#include <cmath>
// [[Rcpp::depends(RcppArmadillo, RcppDist)]]
using namespace arma;
using namespace Rcpp;
using namespace std;


// [[Rcpp::export]]
arma::mat Lambda_func_rcpp(arma::mat dist_val, double eta2_val) {
  return(exp( -dist_val / (2 * eta2_val) ));
}


// [[Rcpp::export]]
arma::vec r_mvnorm_rcpp(arma::vec mu, arma::mat sigma) {
   int ncols = sigma.n_cols;
   arma::vec Y = arma::randn(ncols);
   return mu + arma::chol(sigma).t() * Y;
}

// [[Rcpp::export]]
arma::mat r_mvnorm_rcpp_n(int n, arma::vec mu, arma::mat sigma) {
   int ncols = sigma.n_cols;
   arma::mat Y = arma::randn(n, ncols);
   return arma::repmat(mu, 1, n).t() + Y * arma::chol(sigma);  // dim: n * ncols
}


// [[Rcpp::export]]
List GenerateYmat(int p, int G, int ni, arma::mat& dists,
                   arma::mat& mu_true_mat, arma::mat& gamma_true_mat, arma::mat& delta2_true_mat,
                   arma::mat& eta2_true_mat, List& true_labels,
                   double alpha_g_true, double sigma2_g) {
  List Y_mat(p);

  // #pragma omp parallel for shared(p, G) private(i, g)
  for (int i = 0; i < p; ++i) {
    arma::mat tmpMat(ni, G);
    for (int g = 0; g < G; ++g) {
      arma::mat tmpDelta;
      tmpDelta = delta2_true_mat(g, i) * Lambda_func_rcpp(dists, eta2_true_mat(g, i));
      arma::mat tmpSigma = tmpDelta * sigma2_g;
      arma::rowvec tmp_mu_g = mu_true_mat.row(g);
      arma::rowvec tmp_gamma_g = gamma_true_mat.row(g);
      arma::urowvec tmpIds = as<arma::urowvec>(wrap(true_labels[i])) - 1;
      arma::vec tmpMean = alpha_g_true + tmp_mu_g(tmpIds) + tmp_gamma_g(i);
      tmpMat.col(g) = r_mvnorm_rcpp(tmpMean, tmpSigma);
    }
    Y_mat[i] = tmpMat.t();
  }

  return Y_mat;
}


// [[Rcpp::export]]
List GenerateYmat_AllTrueVals(int p, int G, int ni, arma::mat& dists,
                   arma::vec& mu_true_vec, arma::vec& gamma_true_vec, arma::vec& delta2_true_vec,
                   arma::vec& eta2_true_vec, List& true_labels,
                   double alpha_g_true, double sigma2_g) {
  List Y_mat(p);

  // #pragma omp parallel for shared(p, G) private(i, g)
  for (int i = 0; i < p; ++i) {
    // arma::mat tmpMat(ni, G);
    arma::mat tmpDelta;
    tmpDelta = delta2_true_vec(i) * Lambda_func_rcpp(dists, eta2_true_vec(i));
    arma::mat tmpSigma = tmpDelta * sigma2_g;
    // arma::rowvec tmp_mu_g = mu_true_mat.row(g);
    // arma::rowvec tmp_gamma_g = gamma_true_mat.row(g);
    arma::uvec tmpIds = as<arma::uvec>(wrap(true_labels[i])) - 1;
    arma::vec tmpMean = alpha_g_true + mu_true_vec(tmpIds) + gamma_true_vec(i);
    // cout << tmpIds << endl;
    // cout << mu_true_vec << endl;
    // cout << mu_true_vec(tmpIds) << endl;
    Y_mat[i] = r_mvnorm_rcpp_n(G, tmpMean, tmpSigma);
    // Y_mat[i] = tmpMat.t();
  }

  return Y_mat;
}