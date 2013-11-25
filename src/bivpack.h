#ifndef RCPPDPM_H
#define RCPPDPM_H

#include <sstream>
#include <string>
#include <RcppEigen.h>

using namespace std;
using namespace Eigen;
using namespace Rcpp;

// U(a, b)
double runif1d(double a=0., double b=1.);

// iid U(a, b)
VectorXd runifXd(int size, double a=0., double b=1.);

inline VectorXd runif(int size, double a=0., double b=1.) {
  return as<VectorXd>(Rcpp::runif(size, a, b));
}

// can't overload based on return type
//inline const VectorXd runif(int size, double a=0., double b=1.) {
//  return as< Map<VectorXd> >(Rcpp::runif(size, a, b));
//}

// N(mu, sd^2)
double rnorm1d(double mu=0., double sd=1.);

// iid N(mu, sd^2)
Vector2d rnorm2d(double mu=0., double sd=1.);

VectorXd rnormXd(int size, double mu=0., double sd=1.);

// MVN(mu, Sigma)
VectorXd rnormXd(VectorXd mu, MatrixXd Sigma);

Vector2d rnorm2d(Vector2d mu, Matrix2d Sigma);

// for truncation, e.g. I(-Inf, b) or I(a, Inf), see below
const double Inf=1./0., NaN=0./0.;

// N(0, 1)I(a, b)
double rtnorm(double a, double b);

// N(mu, sd^2)I(a, b)
double rtnorm(double a, double b, double mu, double sd=1.);

// Gam(shape, rate) where mean=shape/rate
double rgamma1d(double shape, double rate=1.);

// iid Gam(shape, rate) where mean=shape/rate
VectorXd rgammaXd(int size, double shape, double rate=1.);

// Chisq(df)
double rchisq1d(double df);

// iid Chisq(df)
VectorXd rchisqXd(int size, double df);

// W(Omega, nu) where mean=nu*Omega
MatrixXd rwishart(MatrixXd Omega, double nu);

MatrixXd rwishart(double nu, MatrixXd Omega);

// IW(Omega, nu)
MatrixXd rinvwishart(MatrixXd Omega, double nu);

MatrixXd rinvwishart(double nu, MatrixXd Omega);

// return a draw from Multinomial(1, prob)
// e.g. c(1:h) %*% rmultinom(1, 1, prob)
int rmultinom1(VectorXd prob);

// this is a little faster, but ...
// you need to remember to use transpose appropriately since it is not symmetric
Matrix2d inv_root_chol(Matrix2d &A);
  
MatrixXd inv_root_chol(MatrixXd &A);

// this is a little slower, but it is symmetric
Matrix2d inv_root_svd(Matrix2d &A);
  
MatrixXd inv_root_svd(MatrixXd &A);

double JKB(double rho);

double JKB_inverse(double zeta);

// returns zeta rather than rho; also prev=prev.zeta
double corrMH(double a, double b, int N, double prev, int m=10);

// discretely samples rho; preferable to corrMH
double corrMN(double a, double b, int N, int n=100);

// uses shifted/scaled Beta prior, i.e. 2*beta-1
//double corrBetaMH(double a, double b, int N, double prev, int m=5);

Vector2d rbvtruncnorm(Vector2d value, int quadrant, 
		      double rho, int burnin=0, int M=1, int thin=1);
/*
MatrixXd erase(MatrixXd A, int b);

VectorXi freq(VectorXi C);
*/

// mean of A applied over the columns: row vector of means returned
RowVectorXd mean(MatrixXd A);

/*
// mean of A applied over the columns: vector of means returned
VectorXd mean(MatrixXd A);
*/

MatrixXd SSCP(MatrixXd A);

double bnivLDF(double arg, double *myldfd);

RcppExport SEXP nniv(SEXP info, SEXP data, SEXP mcmc);

RcppExport SEXP bniv(SEXP info, SEXP data, SEXP mcmc);

RcppExport SEXP bnivDPM(SEXP info, SEXP data, SEXP mcmc);

RcppExport SEXP bbiv(SEXP info, SEXP data, SEXP mcmc);

//#define DEBUG_NEAL8 1
RcppExport SEXP bbivDPM(SEXP info, SEXP data, SEXP mcmc);

/*
double bvpF(VectorXd y, VectorXd phi, List prior);
VectorXd bvpG0(List prior);
VectorXd bvpP0(MatrixXd Y, List prior, VectorXd phi);
double bvp_alpha(int k, int N, double alpha, List alpha_prior);
*/

List neal8(const MatrixXd Y, VectorXi C, MatrixXd phi, VectorXi states, 
	   const int m, const double alpha, const List prior, 
	   double (*F)(VectorXd y, VectorXd phi, List prior),
	   VectorXd (*G0)(List prior),
	   VectorXd (*P0)(MatrixXd y, List prior, VectorXd phi));

#endif
