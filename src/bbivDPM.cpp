#include "bivpack.h"

// alpha random
double bbiv_alpha(int k, int N, double alpha, List alpha_prior) {
  double a=as<double>(alpha_prior["a"]), b=as<double>(alpha_prior["b"]), 
    nu=R::rbeta(alpha+1., N), u=runif1d(), OR=(a+k-1)/(N*(b-log(nu))),
    p=OR/(1.+OR);

  alpha = (u<p) ? rgamma1d(a+k, b-log(nu)) : rgamma1d(a+k-1., b-log(nu));

  return alpha;
}

// phi contains mu1, mu2, rho
// prior contains mu0, T0 and S0
double bbivF(VectorXd eps, VectorXd phi, List prior) {
  // BVN density in terms of rho 
  const double rho=phi[2], rho2=pow(rho, 2.);

  MatrixXd T;
  T.setConstant(2, 2, pow(1.-rho2, -1.));
  T(0, 1) *= -rho;
  T(1, 0) *= -rho;

  VectorXd mu=phi;
  mu.conservativeResize(2);
  eps.conservativeResize(2); // semi

  eps -= mu; 

  return pow(1.-rho2, -0.5)*exp(-0.5*(eps.transpose()*T*eps)[0]);
}

VectorXd bbivG0(List prior) {
  VectorXd phi=rnorm2d(as< Map<VectorXd> >(prior["mu0"]), 
		       as< Map<MatrixXd> >(prior["S0"]));

  phi.conservativeResize(3);

  //phi[2]=runif1d(0., 1.);//crazy check 3
  phi[2]=runif1d(-1., 1.);

  return phi;
}

VectorXd bbivP0(MatrixXd A, List prior, VectorXd phi) {
  const int n=A.rows(); //, rho_MH=as<int>(prior["rho_MH"]);

#ifdef DEBUG_NEAL8
  Rcout << "\npre:"  << phi;
#endif

  /* this Gibbs conditional is WRONG; see semi block below
  RowVectorXd mu=phi;
  mu.conservativeResize(2);

  MatrixXd D(n, 2);

  for(int i=0; i<n; ++i) D.row(i) = eps.row(i)-mu;

  MatrixXd DtD=D.transpose()*D;

  double rho=phi[2], zeta=corrMH(DtD(0, 0)+DtD(1, 1), DtD(0, 1), n, JKB(rho));

  rho=JKB_inverse(zeta);
  */

  MatrixXd eps=A.block(0, 0, n, 2); // semi

  const double rho=phi[2], zeta=JKB(rho); // semi

  Vector2d mean_eps=mean(eps).transpose(), mu0=as< Map<VectorXd> >(prior["mu0"]);

  Matrix2d T0=as< Map<MatrixXd> >(prior["T0"]), Tj, Sigma, Tprec;

  Tj << 1., -rho, -rho, 1.;

  Tj *= (double(n)/(1.-pow(rho, 2.)));

  Tprec=T0+Tj;

  Sigma=Tprec.inverse();

  Vector2d mean_mu=Sigma*(T0*mu0+Tj*mean_eps);

  phi=rnorm2d(mean_mu, Sigma);
  phi.conservativeResize(3);

  /* semi block begins */
  const double beta=as<double>(prior["beta"]);

  const VectorXd t=A.col(2), y=A.col(3), 
    gamma=as< Map<VectorXd> >(prior["gamma"]), 
    delta=as< Map<VectorXd> >(prior["delta"]), 
    eta  =as< Map<VectorXd> >(prior["eta"]); 

  const int p=gamma.size(), q=delta.size();

  MatrixXd D(n, 2), X=A.block(0, 4, n, p), Z=A.block(0, p+4, n, q);

  D.col(0).setConstant(phi[0]);
  D.col(1).setConstant(phi[1]);

  D.col(0) += (X*gamma + Z*delta);  //mutmargin
  D.col(1) += (beta*D.col(0)+X*eta);//muymargin

  D.col(0) = t - D.col(0);                //tdev=t-mutmargin        
  D.col(1) = y - D.col(1) - beta*D.col(0);//ycondtdev=y-muymargin-b*tdev

  Matrix2d SSCP=D.transpose()*D;

  phi[2]=corrMN(SSCP(0, 0)+SSCP(1, 1), SSCP(0, 1), n);

// #ifdef CORRBETAMH
//   phi[2]=corrBetaMH(SSCP(0, 0)+SSCP(1, 1), SSCP(0, 1), n, rho);
// #else
//   phi[2]=JKB_inverse(corrMH(SSCP(0, 0)+SSCP(1, 1), SSCP(0, 1), n, zeta, rho_MH));
// #endif
  /* semi block ends */

#ifdef DEBUG_NEAL8
  Rcout << "\npost:"  << phi   << "\nbeta:"  << beta 
	<< "\ngamma:" << gamma << "\ndelta:" << delta << "\neta:" << eta 
	<< "\nn:"     << n     << "\nA:\n"   << A     << '\n';
#endif
  
  return phi;
}

RcppExport SEXP bbivDPM(SEXP arg1, SEXP arg2, SEXP arg3) {
  // 3 arguments
  // arg1 for parameters
  // arg2 for data
  // arg3 for Gibbs

  // data
  List list2(arg2); 

  const MatrixXd X=as< Map<MatrixXd> >(list2["X"]),
    Z=as< Map<MatrixXd> >(list2["Z"]);

  const VectorXi v1=as< Map<VectorXi> >(list2["tbin"]),
    v2=as< Map<VectorXi> >(list2["ybin"]);

  const int N=X.rows(), p=X.cols(), q=Z.cols(), r=p+q, s=p+r;

#ifdef DEBUG_NEAL8
  List P, Phi, B;
  VectorXi S, one;
  one.setConstant(N, 1);
#endif

  // parameters
  List list1(arg1), 
    beta_info=list1["beta"], 
    rho_info=list1["rho"], 
    mu_info=list1["mu"], 
    theta_info=list1["theta"],
    dpm_info=list1["dpm"], // DPM
    alpha_info=dpm_info["alpha"], // alpha random
    alpha_prior; // alpha random

  const int m=as<int>(dpm_info["m"]); // DPM
  //  const double alpha=as<double>(dpm_info["alpha"]); // DPM
  const int alpha_fixed=as<int>(alpha_info["fixed"]); // alpha random
  double alpha=as<double>(alpha_info["init"]); // alpha random

  if(alpha_fixed==0) alpha_prior=alpha_info["prior"]; // alpha random

  VectorXi C=as< Map<VectorXi> >(dpm_info["C"]), 
    states=as< Map<VectorXi> >(dpm_info["states"]);
  
  /* checks done in neal8
  if(states.sum()!=N || C.size()!=N) { // limited reality check of C and states
    C.setConstant(N, 0); 
    states.setConstant(1, N); 
  }
  */

  // prior parameters
  List beta_prior=beta_info["prior"],
    // no prior parameters for rho
    dpm_prior=dpm_info["prior"],
    theta_prior=theta_info["prior"];

  const double beta_prior_mean=as<double>(beta_prior["mean"]);
  const double beta_prior_prec=as<double>(beta_prior["prec"]);

  const VectorXd theta_prior_mean=as< Map<VectorXd> >(theta_prior["mean"]);
  const MatrixXd theta_prior_prec=as< Map<MatrixXd> >(theta_prior["prec"]);

  // initialize parameters
  double beta    =as<double>(beta_info["init"]); 
  double rho_init=as<double>(rho_info["init"]);  // DPM
  //int    rho_MH  =as<int>(rho_info["MH"]);

  Vector2d mu_init=as< Map<VectorXd> >(mu_info["init"]); // DPM
  VectorXd theta  =as< Map<VectorXd> >(theta_info["init"]);
  VectorXd rho(N);  // DPM

  /*
  VectorXi C, states; // DPM
  C.setConstant(N, 0); // DPM
  states.setConstant(1, N); // DPM
  */

  MatrixXd mu(N, 2), phi(1, 3); // DPM
  phi(0, 0)=mu_init[0]; // DPM
  phi(0, 1)=mu_init[1]; // DPM
  phi(0, 2)=rho_init; // DPM

  // Gibbs
  List list3(arg3);

  const int burnin=as<int>(list3["burnin"]), M=as<int>(list3["M"]), 
    thin=as<int>(list3["thin"]);

  VectorXi quadrant(N);

  VectorXd t(N), y(N); // latents

  // prior parameter intermediate values
  double beta_prior_prod=beta_prior_prec * beta_prior_mean;

  VectorXd theta_prior_prod=theta_prior_prec * theta_prior_mean;

  Matrix2d Sigma, Tprec, B_inverse;

  B_inverse.setIdentity();

  VectorXd gamma=theta.segment(0, p), 
    delta=theta.segment(p, q), 
    eta  =theta.segment(r, p);

  MatrixXd eps(N, 2), D(N, 2), theta_cond_var_root(s, s), W(2, s), A(N, r+4); // DPM semi

  W.setZero();

  MatrixXd theta_cond_prec(s, s);

  VectorXd theta_cond_prod(s), w(r), mu_t(N), mu_y(N), sd_y(N);

  Vector2d u, R, mu_u; // DPM

  double beta_prec, beta_prod, beta_cond_var, beta_cond_mean, beta2; // DPM

  int h=0, i, l; 

  List GS(M); //DPM

  // assign quadrants
  for(i=0; i<N; ++i) {
    if(v1[i]==0 && v2[i]==0) quadrant[i] = 3;
    else if(v1[i]==0 && v2[i]==1) quadrant[i] = 2;
    else if(v1[i]==1 && v2[i]==0) quadrant[i] = 4;
    else if(v1[i]==1 && v2[i]==1) quadrant[i] = 1;
  }

  // Gibbs loop
  //for(int l=-burnin; l<=(M-1)*thin; ++l) {

  l=-burnin;

  do{
    // populate mu/rho //DPM
    for(i=0; i<N; ++i) {
      mu(i, 0)=phi(C[i], 0);
      mu(i, 1)=phi(C[i], 1);
      rho[i]=phi(C[i], 2);

      mu_t[i]=mu(i, 0);
      mu_y[i]=mu(i, 1);
    }

    // generate latents
    // mu_t = mu.col(0); //DPM
    mu_t += (Z*delta + X*gamma);

    // mu_y = mu.col(1); //DPM
    mu_y += (beta*mu_t + X*eta);

    beta2=pow(beta, 2.);

    for(i=0; i<N; ++i) {
      sd_y[i] = sqrt(beta2+2.*beta*rho[i]+1.); //DPM

      mu_u[0]=mu_t[i];
      mu_u[1]=mu_y[i]/sd_y[i]; //DPM
      //             z,    quadrant,    rho,                   burnin 
      u=rbvtruncnorm(mu_u, quadrant[i], (beta+rho[i])/sd_y[i], 10);
  
      t[i]=u[0];
      y[i]=sd_y[i]*u[1];
    }

    // sample beta
    D.col(0) = (t - mu.col(0) - X*gamma - Z*delta); //DPM
    D.col(1) = (y - mu.col(1) - X*eta); //DPM

    beta_prec=0.;
    beta_prod=0.;

    for(i=0; i<N; ++i) {
      double Sigma_det=1.-pow(rho[i], 2.); //DPM

      beta_prec += pow(t[i], 2.)/Sigma_det; //DPM
      beta_prod += -t[i]*(rho[i]*D(i, 0)-D(i, 1))/Sigma_det; //DPM
    }

    beta_cond_var=1./(beta_prec+beta_prior_prec);

    beta_cond_mean=beta_cond_var*(beta_prod+beta_prior_prod);

    beta=rnorm1d(beta_cond_mean, sqrt(beta_cond_var));

    B_inverse(1, 0)=-beta;

    // sample theta
    theta_cond_prec=theta_prior_prec;
    theta_cond_prod=theta_prior_prod;

    for(i=0; i<N; ++i) {
      double Sigma_det=1.-pow(rho[i], 2.); //DPM

      Tprec(0, 0)=1./Sigma_det;      Tprec(0, 1)=-rho[i]/Sigma_det; //DPM
      Tprec(1, 0)=-rho[i]/Sigma_det; Tprec(1, 1)=1./Sigma_det; //DPM

      mu_u[0]=mu(i, 0);  //DPM
      mu_u[1]=mu(i, 1);  //DPM

      W.block(0, 0, 1, p)=X.row(i);
      W.block(0, p, 1, q)=Z.row(i);
      W.block(1, r, 1, p)=X.row(i);

      theta_cond_prec += (W.transpose() * Tprec * W);

      u[0]=t[i];
      u[1]=y[i];

      R=B_inverse*u-mu_u; //DPM

      theta_cond_prod += (W.transpose() * Tprec * R);
    }

    theta_cond_var_root=inv_root_chol(theta_cond_prec);

    theta=theta_cond_var_root*(rnormXd(s)+theta_cond_var_root.transpose()*theta_cond_prod);

    gamma=theta.segment(0, p); 
    delta=theta.segment(p, q); 
    eta  =theta.segment(r, p);

    // sample mu and rho
    // this for block should be placed in P0
    // however, to keep changes minimal, we keep it here
    for(i=0; i<N; ++i) {

      W.block(0, 0, 1, p)=X.row(i);
      W.block(0, p, 1, q)=Z.row(i);
      W.block(1, r, 1, p)=X.row(i);

      u[0]=t[i];
      u[1]=y[i];

      eps.row(i) = (B_inverse*u - W*theta).transpose(); //DPM
    }

    /* semi block begins */
    A.block(0, 0, N, 2)=eps;
    A.col(2)=t;
    A.col(3)=y;
    A.block(0, 4, N, p)=X;
    A.block(0, p+4, N, q)=Z;

    List psi=List::create(Named("mu0")=as< Map<VectorXd> >(dpm_prior["mu0"]), 
			  Named("T0")=as< Map<MatrixXd> >(dpm_prior["T0"]),   
			  Named("S0")=as< Map<MatrixXd> >(dpm_prior["S0"]),  
			  Named("beta")=beta,   
			  Named("gamma")=gamma, 
			  Named("delta")=delta, 
			  Named("eta")=eta);   

    if(alpha_fixed==0) 
      alpha=bbiv_alpha(states.size(), N, alpha, alpha_prior); // alpha random

    List dpm_step=neal8(A, C, phi, states, m, alpha, psi, 
     			&bbivF, &bbivG0, &bbivP0);
    /* semi block end */

    /*
    C=as< Map<VectorXi> >(dpm_step[0]);
    phi=as< Map<MatrixXd> >(dpm_step[1]);
    states=as< Map<VectorXi> >(dpm_step[2]);
    */

    C=as< Map<VectorXi> >(dpm_step["C"]);
    phi=as< Map<MatrixXd> >(dpm_step["phi"]);
    states=as< Map<VectorXi> >(dpm_step["states"]);

#ifdef DEBUG_NEAL8
    S=as< Map<VectorXi> >(dpm_step["S"]);
    P=dpm_step["P"]; 
    Phi=dpm_step["Phi"];
    B=dpm_step["B"]; 
#endif
    
    if(l>=0 && l%thin == 0) {
      h = (l/thin);

#ifdef DEBUG_NEAL8
      GS[h]=List::create(Named("beta")=beta, Named("theta")=theta,
			 Named("C")=C+one, Named("phi")=phi, 
			 Named("states")=states, Named("alpha")=alpha, 
			 Named("S")=S+one, Named("P")=P, 
			 Named("Phi")=Phi, Named("B")=B);
#else
      if(alpha_fixed==0) // alpha random
	GS[h]=List::create(Named("beta")=beta, Named("theta")=theta, 
			   Named("C")=C, Named("phi")=phi, 
			   Named("states")=states, Named("alpha")=alpha);
      else GS[h]=List::create(Named("beta")=beta, Named("theta")=theta, 
                       Named("C")=C, Named("phi")=phi, Named("states")=states);
#endif

      // GS[h]=List::create(Named("beta")=beta, Named("theta")=theta, 
      //  			 Named("C")=C, Named("phi")=phi, Named("states")=states,
      // 			 Named("m")=m, Named("alpha")=alpha, Named("psi")=psi);
    }

    l++;

  } while (l<=(M-1)*thin); 

  return wrap(GS);
}
