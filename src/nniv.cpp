#include "bivpack.h"

RcppExport SEXP nniv(SEXP arg1, SEXP arg2, SEXP arg3) {
  // 3 arguments
  // arg1 for parameters
  // arg2 for data
  // arg3 for Gibbs

  // data
  List list2(arg2); 

  const MatrixXd X=as< Map<MatrixXd> >(list2["X"]),
    Z=as< Map<MatrixXd> >(list2["Z"]);

  const VectorXd t=as< Map<VectorXd> >(list2["t"]),
    y=as< Map<VectorXd> >(list2["y"]);

  const int N=X.rows(), p=X.cols(), q=Z.cols(), r=p+q, s=p+r;

  // parameters
  List list1(arg1), 
    beta_info=list1["beta"], 
    Tprec_info=list1["Tprec"], 
    mu_info=list1["mu"], 
    theta_info=list1["theta"];

  // prior parameters
  List beta_prior=beta_info["prior"],
    Tprec_prior=Tprec_info["prior"],
    mu_prior=mu_info["prior"],
    theta_prior=theta_info["prior"];

  const double Tprec_prior_nu=as<double>(Tprec_prior["nu"]);
  const Matrix2d Tprec_prior_Psi=as< Map<MatrixXd> >(Tprec_prior["Psi"]);

  const double beta_prior_mean=as<double>(beta_prior["mean"]);
  const double beta_prior_prec=as<double>(beta_prior["prec"]);

  const Vector2d mu_prior_mean=as< Map<VectorXd> >(mu_prior["mean"]);
  const Matrix2d mu_prior_prec=as< Map<MatrixXd> >(mu_prior["prec"]);

  const VectorXd theta_prior_mean=as< Map<VectorXd> >(theta_prior["mean"]);
  const MatrixXd theta_prior_prec=as< Map<MatrixXd> >(theta_prior["prec"]);

  // initialize parameters
  double beta=as<double>(beta_info["init"]); 

  Matrix2d Tprec=as< Map<MatrixXd> >(Tprec_info["init"]);
  Vector2d mu   =as< Map<VectorXd> >(mu_info["init"]);
  VectorXd theta=as< Map<VectorXd> >(theta_info["init"]);

  // Gibbs
  List list3(arg3); //, save=list3["save"];

  const int burnin=as<int>(list3["burnin"]), M=as<int>(list3["M"]), 
    thin=as<int>(list3["thin"]), m=7+s;

  MatrixXd GS(M, m);

  // prior parameter intermediate values
  double beta_prior_prod=beta_prior_prec * beta_prior_mean;

  VectorXd theta_prior_prod=theta_prior_prec * theta_prior_mean;

  Vector2d mu_prior_prod=mu_prior_prec*mu_prior_mean;

  // parameter intermediate values
  Matrix2d Sigma, B_inverse;

  Sigma=Tprec.inverse();
  B_inverse.setIdentity();

  VectorXd gamma=theta.segment(0, p), 
    delta=theta.segment(p, q), 
    eta  =theta.segment(r, p);

  /*
    MatrixXd Theta(2, r);

    Theta.row(0)=theta.segment(0, r);
    Theta.bottomLeftCorner(1, q)=RowVectorXd::Zero(q);
    Theta.bottomRightCorner(1, p)=eta.transpose();
  */

  Vector2d eps, eps_sum;

  MatrixXd D(N, 2), theta_cond_var_root(s, s), W(2, s);

  W.setZero();

  MatrixXd theta_cond_prec(s, s);

  VectorXd theta_cond_prod(s), w(r);

  Matrix2d mu_cond_prec, mu_cond_var_root, E;

  Vector2d u, R, mu_cond_prod, mu_u;

  double beta_scale, beta_prec, beta_prod, beta_cond_var, beta_cond_mean;

  int h=0, i, l; 


  // Gibbs loop
  //for(int l=-burnin; l<=(M-1)*thin; ++l) {

  l=-burnin;

  do{
    // sample beta
    D.col(0).setConstant(-mu[0]);
    D.col(1).setConstant(-mu[1]);

    D.col(0) += (t - X*gamma - Z*delta);
    D.col(1) += (y - X*eta);

    beta_scale=1./(Sigma(0, 0)*t.dot(t));

    beta_prec=1./(beta_scale*Sigma.determinant());

    beta_prod=beta_prec*beta_scale
      *((Sigma(0, 0)*D.col(1)-Sigma(0, 1)*D.col(0)).array()*t.array()).sum();

    beta_cond_var=1./(beta_prec+beta_prior_prec);

    beta_cond_mean=beta_cond_var*(beta_prod+beta_prior_prod);

    beta=rnorm1d(beta_cond_mean, sqrt(beta_cond_var));

    B_inverse(1, 0)=-beta;

    // sample theta
    theta_cond_prec=theta_prior_prec;
    theta_cond_prod=theta_prior_prod;

    for(i=0; i<N; ++i) {
      /*
	W.topLeftCorner(1, p)=X.row(i);
	W.block(0, p, 1, q)=Z.row(i);
	W.bottomRightCorner(1, p)=X.row(i);
      */

      W.block(0, 0, 1, p)=X.row(i);
      W.block(0, p, 1, q)=Z.row(i);
      W.block(1, r, 1, p)=X.row(i);

      theta_cond_prec += (W.transpose() * Tprec * W);

      u[0]=t[i];
      u[1]=y[i];

      R=B_inverse*u-mu;

      theta_cond_prod += (W.transpose() * Tprec * R);
    }

    theta_cond_var_root=inv_root_chol(theta_cond_prec);
    // theta_cond_var_root=inv_root_svd(theta_cond_prec); // for validation only

    theta=theta_cond_var_root*(rnormXd(s)+theta_cond_var_root.transpose()*theta_cond_prod);

    gamma=theta.segment(0, p); 
    delta=theta.segment(p, q); 
    eta  =theta.segment(r, p);

    /*
      Theta.topRows(1)=theta.segment(0, r).transpose();
      Theta.bottomRightCorner(1, p)=eta.transpose();
    */

    // sample mu
    eps_sum.setZero();
    //W.setZero();

    E.setZero();

    for(i=0; i<N; ++i) {
      /*
	W.topLeftCorner(1, p)=X.row(i);
	W.block(0, p, 1, q)=Z.row(i);
	W.bottomRightCorner(1, p)=X.row(i);
      */

      W.block(0, 0, 1, p)=X.row(i);
      W.block(0, p, 1, q)=Z.row(i);
      W.block(1, r, 1, p)=X.row(i);

      /*
	w.segment(0, q)=Z.row(i);
	w.segment(q, p)=X.row(i);
      */

      u[0]=t[i];
      u[1]=y[i];

      //eps += B_inverse*u - Theta*w;
      eps = B_inverse*u - W*theta;
      eps_sum += eps;
      eps -= mu;
      E += eps*eps.transpose();
    }

    mu_cond_prod=Tprec*eps_sum+mu_prior_prod;

    mu_cond_prec=(N*Tprec+mu_prior_prec);

    mu_cond_var_root=inv_root_chol(mu_cond_prec);
    // mu_cond_var_root=inv_root_svd(mu_cond_prec); // for validation only

    mu=mu_cond_var_root*(rnormXd(2)+mu_cond_var_root.transpose()*mu_cond_prod);

    // sample Tprec
    Tprec = rwishart((E+Tprec_prior_Psi).inverse(), N+Tprec_prior_nu);
    Sigma = Tprec.inverse();

    if(l>=0 && l%thin == 0) {
      h = (l/thin);

      GS.block(h, 0, 1, s)=theta.transpose();

      GS(h, s)=beta;
      GS(h, s+1)=mu[0];
      GS(h, s+2)=mu[1];
      GS(h, s+3)=Tprec(0, 0);
      GS(h, s+4)=Tprec(0, 1);
      GS(h, s+5)=Tprec(0, 1);
      GS(h, s+6)=Tprec(1, 1);
    }

    l++;

  } while (l<=(M-1)*thin && beta==beta); 

  if(beta != beta) GS.conservativeResize(h+1, m);

  return wrap(GS);
}
