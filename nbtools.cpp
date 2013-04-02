#include <RcppArmadillo.h>
using namespace Rcpp;

arma::mat prodonrow(arma::mat X, arma::colvec w)
{
  arma::mat wX = X;
  int n = X.n_rows;
  for(int i=0; i<n; i++)
    {
      wX.row(i) *= w(i);
    }
  return wX;
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List drawloadingsBlockC(NumericMatrix Yc, NumericMatrix Fc, NumericMatrix Omegac, NumericVector Fmeansc, NumericVector Intc, double overdispr, NumericMatrix priorprecr) {
  // each row of Yr(n,p) is a vector of observations for a single unit
  // each row of Fr(n,k) is a vector of factor scores for a single unit
  // each row of Omegar(n,p) is a vector of latent variables, corresponding to Yr'
  // Fmeans r is a vector of unit-level means (e.g. document lengths)
  // Intr is a vector of feature-level means (e.g. global word frequencies)
  // overdispr is the overdispersion parameter in the negative-binomial model
  // priorprecr is the prior precision matrix for each row of the loadings matrix
    int numobs = Yc.nrow(), numfeatures = Yc.ncol(), numfactors = Fc.ncol();
    double xi = overdispr;

    arma::mat Y(Yc.begin(), numobs, numfeatures, false);       // reuses memory and avoids extra copy
    arma::mat Omega(Omegac.begin(), numobs, numfeatures, false);
    arma::mat F(Fc.begin(), numobs, numfactors, false);
    arma::colvec Intercept(Intc.begin(), numfeatures, false);
    arma::colvec Fmeans(Fmeansc.begin(), numobs, false);

    arma::mat loadings(numfeatures, numfactors); loadings.fill(0.0);
    arma::colvec om;
    arma::colvec kap;
    arma::colvec offset;
    arma::colvec overdisp(numobs); overdisp.fill(xi);
    arma::mat OmF(numobs,numfactors);
    arma::mat tFOmF(numfactors, numfactors);
    arma::mat PostPrecision(numfactors, numfactors);
    for(int j=0; j<numfeatures; j++)
      {
	kap = (Y.col(j) - overdisp)/2.0;
	om = Omega.col(j);
	OmF = prodonrow(F,om);
	tFOmF = F.t() * OmF;
	PostPrecision = tFOmF;
	offset = prodonrow(Fmeans, om) + Intercept(j)*om;
	kap -= offset;
	loadings.row(j) = solve(PostPrecision, F.t() * kap).t();
      }
    return Rcpp::List::create(
			      Rcpp::Named("loadings") = loadings
    ) ;
}
