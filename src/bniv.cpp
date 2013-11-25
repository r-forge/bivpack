#include "arms.h"
	
double bnivLDF(double arg, double *myldfd) {
  //alpha=myldfd[0], lambda=myldfd[1], c=myldfd[2]		
  return (2.*myldfd[0]-1.)*log(arg)-myldfd[1]*pow(arg, 2.)+myldfd[2]*arg;
}

RcppExport SEXP bniv(SEXP arg1, SEXP arg2, SEXP arg3) {
  //3 arguments
  //arg1 for parameters
  //arg2 for data
  //arg3 for parameters

  //data
  List list2(arg2);

  const MatrixXd X=as< Map<MatrixXd> >(list2["X"]),
    Z=as< Map<MatrixXd> >(list2["Z"]);

  const VectorXi tbin=as< Map<VectorXi> >(list2["tbin"]);

  const VectorXd y=as< Map<VectorXd> >(list2["y"]);

  const int N=X.rows(), p=X.cols(), q=Z.cols(), r=p+q, s=p+r;

  // parameters
  List list1(arg1),
    beta_info=list1["beta"],
    rho_info=list1["rho"],
    mu_info=list1["mu"],
    tau_info=list1["tau"],
    theta_info=list1["theta"];

  // prior parameters

  List beta_prior=beta_info["prior"],
    //rho_prior=rho_info["prior"],
    mu_prior=mu_info["prior"],
    tau_prior=tau_info["prior"],
    theta_prior=theta_info["prior"];

  const double beta_prior_mean=as<double>(beta_prior["mean"]);
  const double beta_prior_prec=as<double>(beta_prior["prec"]);

  const double tau_prior_alpha0=as<double>(tau_prior["alpha0"]);
  const double tau_prior_lambda0=as<double>(tau_prior["lambda0"]);

  const Vector2d mu_prior_mean=as< Map<VectorXd> >(mu_prior["mean"]);
  const Matrix2d mu_prior_prec=as< Map<MatrixXd> >(mu_prior["prec"]);

  const VectorXd theta_prior_mean=as< Map<VectorXd> >(theta_prior["mean"]);
  const MatrixXd theta_prior_prec=as< Map<MatrixXd> >(theta_prior["prec"]);

  // initialize parameters

  double beta=as<double>(beta_info["init"]);
  double rho=as<double>(rho_info["init"]);
  //double rho=0.;
  double tau=as<double>(tau_info["init"]);

  Vector2d mu=as< Map<VectorXd> >(mu_info["init"]);
  VectorXd theta=as< Map<VectorXd> >(theta_info["init"]);

  //Gibbs

  List list3(arg3); //, save=list3["save"];
  const int burnin=as<int>(list3["burnin"]), M=as<int>(list3["M"]), thin=as<int>(list3["thin"]),m=5+s;
  MatrixXd GS(M,m);
  VectorXd t(N); //latent variable

  //prior parameter intermediate values

  double beta_prior_prod=beta_prior_prec*beta_prior_mean;
  VectorXd theta_prior_prod=theta_prior_prec*theta_prior_mean;
  VectorXd mu_prior_prod=mu_prior_prec*mu_prior_mean;

  // parameter intermediate values

  double zeta=JKB(rho);

  Matrix2d Sigma, Tprec, B_inverse;
  Sigma.setIdentity();
  //double tau=1.;
  double s2sq=1.;
  double s2=sqrt(s2sq);
  double s12=rho*s2;
  Sigma(1,1)=s2sq;
  Sigma(0,1)=s12;
  Sigma(1,0)=s12;
  Tprec=Sigma.inverse();

  VectorXd gamma=theta.segment(0,p),
    delta=theta.segment(p,q),
    eta=theta.segment(r,p);

  //Gibbs loop preop

  //double *myldfd=(double *)malloc(3*sizeof(double));
  double myldfd[3];

  //myldfd.alpha = tau_prior_alpha0+N/2;
  myldfd[0]=tau_prior_alpha0+N/2.;

  MatrixXd dif(N,2);
  //MatrixXd B_inverse;
  B_inverse.setIdentity();
  MatrixXd W(2,s);
  W.setZero();
  MatrixXd theta_cond_var_root(s,s);
  MatrixXd theta_cond_prec(s,s);
  VectorXd theta_cond_prod(s), w(r);

  Matrix2d mu_cond_prec, mu_cond_var_root, SSCP;
  Vector2d u, R, mu_cond_prod, mu_u;

  VectorXd mutmargin(N), muymargin(N),endpoint(N);
  double numerator, denominator, sdtcondy;
  double b_scale,b_prec,b_prod,b_cond_var,b_cond_mean;
  VectorXd mutcondy;
  const double inf=1.0/0.0;
  //double lambda1, c;

  Vector2d eps;

  int h,j,k,l,b;

  l=-burnin;

  // arms prep

  //myldfd.alpha=as<double>(alpha);
  //myldfd.lambda=as<double>(lambda);
  //myldfd.c=as<double>(c);

  double xinits[5]={0.01, 0.1, 1., 10., 100.}; //{1.,2.,3.,4.};

  double lower=0.001, upper=1000.; // upper=inf DOES NOT WORK
  double convex=1.;
  double xprev, samp;
  double qcent[2], xcent[2];
  double roottau;
  int neval, err=0, maxenv=50, metro=0, nsamp=1;
  
BEGIN_RCPP
  do{
    //generate latents
    mutmargin.setConstant(mu[0]);
    mutmargin+=Z*delta+X*gamma;
    muymargin.setConstant(mu[1]);
    muymargin+=(beta*mutmargin+X*eta);
    numerator=beta+rho*sqrt(Sigma(1,1));
    denominator=pow(beta,2)+2*beta*rho*sqrt(Sigma(1,1))+Sigma(1,1);
    mutcondy=mutmargin+(numerator/denominator)*(y-muymargin);
    sdtcondy=sqrt(1-pow(numerator,2)/denominator);
    for(b=0;b<N;++b)
      {
	endpoint[b]=(0-mutcondy[b])/sdtcondy;
      }
	
    //endpoint=(0-mutcondy)/sdtcondy;

    for (j=0; j<N; ++j) {
      if(tbin[j]==0) t[j]=rtnorm(-inf, endpoint[j]);
      else t[j]=rtnorm(endpoint[j], inf);
    }
    t=sdtcondy*t+mutcondy;

    dif.col(0).setConstant(-mu[0]);
    dif.col(0)+=t-X*gamma-Z*delta;
    dif.col(1).setConstant(-mu[1]);
    dif.col(1)+=y-X*eta;

    // sample beta

    b_scale=1./(Sigma(0,0)*t.dot(t));
    b_prec=1./(b_scale*Sigma.determinant());
    b_prod=b_prec*b_scale*((Sigma(0,0)*dif.col(1)-Sigma(0,1)*dif.col(0)).array()*t.array()).sum();
    b_cond_var=1./(b_prec+beta_prior_prec);
    b_cond_mean=b_cond_var*(b_prod+beta_prior_prod);
    beta=rnorm1d(b_cond_mean, sqrt(b_cond_var));
    B_inverse(1,0)=-beta;

    // sample theta
	
    theta_cond_prec=theta_prior_prec;
    theta_cond_prod=theta_prior_prod;
	
    for(k=0; k<N; ++k) {
      W.block(0,0,1,p)=X.row(k);
      W.block(0,p,1,q)=Z.row(k);
      W.block(1,r,1,p)=X.row(k);
		
      theta_cond_prec += (W.transpose() * Tprec * W);
      u[0]=t[k];
      u[1]=y[k];

      R=B_inverse*u-mu;
      theta_cond_prod += (W.transpose()*Tprec * R);
    }

    //theta_cond_var_root=inv_root_chol(theta_cond_prec);
    theta_cond_var_root=inv_root_svd(theta_cond_prec);
    theta=theta_cond_var_root*(rnormXd(s)+theta_cond_var_root.transpose()*theta_cond_prod);

    gamma=theta.segment(0,p);
    delta=theta.segment(p,q);
    eta=theta.segment(r,p);

    // sample mu

    eps.setZero();
    for(k=0; k<N; ++k) {
      W.block(0,0,1,p)=X.row(k);
      W.block(0,p,1,q)=Z.row(k);
      W.block(1,r,1,p)=X.row(k);
      u[0]=t[k];
      u[1]=y[k];

      eps += B_inverse*u-W*theta;
    }

    mu_cond_prod=Tprec*eps+mu_prior_prod;
    mu_cond_prec=(N*Tprec+mu_prior_prec);
    mu_cond_var_root=inv_root_chol(mu_cond_prec);
    mu=mu_cond_var_root*(rnormXd(2)+mu_cond_var_root.transpose()*mu_cond_prod);

    // sample rho

    dif.col(0).setConstant(mu[0]);
    dif.col(0) += (X*gamma +Z*delta); //mutmargin
    dif.col(1).setConstant(mu[1]);
    dif.col(1) += (beta*dif.col(0)+X*eta); //muymargin
    dif.col(0)= t-dif.col(0); //tdev

    dif.col(1)= y-dif.col(1)-beta*dif.col(0); //ycondtdev

    SSCP=dif.transpose()*dif;

    zeta=corrMH(tau*SSCP(1,1)+SSCP(0,0), sqrt(tau)*SSCP(0,1), N, zeta);

    rho=JKB_inverse(zeta);
	
    // simulate tau

    myldfd[1]=dif.col(1).transpose()*dif.col(1);
    myldfd[1]/=(2.*(1.-pow(rho, 2.)));
    myldfd[1]+=tau_prior_lambda0;	

    myldfd[2]=dif.col(1).transpose()*dif.col(0);
    myldfd[2]*=rho;
    myldfd[2]/=(1.-pow(rho, 2.));

    double tmax=(myldfd[2]+sqrt(pow(myldfd[2], 2.)+8.*myldfd[1]*(2.*myldfd[0]-1.)))/(4.*myldfd[1]);

    double delta1=1./((2.*myldfd[0]-1.)/pow(tmax, 2.)-2.*myldfd[1]);

    double t1=max(0.001, tmax+delta1);
    double t2 = tmax-delta1;
	
    xprev=tau;	

    err=arms(&xinits[0], 5, &lower, &upper, &bnivLDF, &myldfd[0], 
	     &convex, maxenv, metro, &xprev, &samp, nsamp, NULL, NULL, 0, &neval);

    if(err!=0) samp=NA_REAL;

    roottau=samp;
    tau=pow(roottau, 2.);

    // reconstruct Sigma

    s2sq=1./tau;
    s2=sqrt(s2sq);
    s12=rho*s2;
    Sigma(1,1)=s2sq;
    Sigma(0,1)=s12;
    Sigma(1,0)=s12;
    Tprec=Sigma.inverse();

    if(l>=0 && l%thin==0) {
      h=l/thin;
      GS.block(h,0,1,s)=theta.transpose();
      GS(h,s)=beta;
      GS(h,s+1)=mu[0];
      GS(h,s+2)=mu[1];
      GS(h,s+3)=rho;
      GS(h,s+4)=s2sq;
			
    }
    ++l;	
  } while(l<=(M-1)*thin && err==0);

// free(myldfd);

  if (err != 0) {
    GS.conservativeResize(h+1, m);

    stringstream stream;
    string string;
  
    stream << err;
    stream >> string;

    throw range_error(string); 
  }

  return wrap(GS);
END_RCPP 	
}
