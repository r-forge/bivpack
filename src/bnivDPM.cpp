#include "arms.h"

double myldfd[3];

// already defined in bniv.cpp
//double bnivLDF(double arg, double *myldfd) {
//   //alpha=myldfd[0], lambda=myldfd[1], c=myldfd[2]		
//   return (2.*myldfd[0]-1.)*log(arg)-myldfd[1]*pow(arg, 2.)+myldfd[2]*arg;
// }

//alpha random
double bniv_alpha(int k, int N, double alpha, List alpha_prior) {
  double a=as<double>(alpha_prior["a"]), b=as<double>(alpha_prior["b"]), 
    nu=R::rbeta(alpha+1., N), u=runif1d(), OR=(a+k-1)/(N*(b-log(nu))),
    p=OR/(1.+OR);

  alpha = (u<p) ? rgamma1d(a+k, b-log(nu)) : rgamma1d(a+k-1., b-log(nu));

  return alpha;
}

//Functions required for DPM

double bnivF(VectorXd eps, VectorXd phi, List prior) {
	//BVN density in terms of rho and tau	
	const double rho=phi[2];
	const double tau=phi[3];
	double rho2=pow(rho,2.);
	
	MatrixXd T;
	T.setConstant(2,2, pow((1./tau)*(1.-rho2), -1.));
	T(0,1) *= -rho/sqrt(tau);
	T(1,0) *= -rho/sqrt(tau);
	T(0,0) *= 1./tau;

	VectorXd mu=phi;
	mu.conservativeResize(2);
	eps.conservativeResize(2); //semi
	eps -= mu;
	
	return pow((1./tau)*(1.-rho2), -0.5) * exp(-0.5*(eps.transpose()*T*eps)[0]);
}

VectorXd bnivG0(List a) {
	VectorXd phi=rnorm2d(as< Map<VectorXd> >(a["mu0"]),
			     as< Map<MatrixXd> >(a["S0"]));

	phi.conservativeResize(4);
	phi[2]=runif1d(-1.,1.);
	phi[3]=rgamma1d(as<double>(a["alpha0"]),as<double>(a["lambda0"]));
	return phi;			
}

VectorXd bnivP0(MatrixXd A, List prior, VectorXd phi) {
	const int n=A.rows();

	MatrixXd eps=A.block(0,0,n,2); //semi

	double rho=phi[2];

	const double zeta=JKB(rho);

	double tau=phi[3];

	Vector2d mean_eps=mean(eps).transpose(), mu0=as< Map<VectorXd> >(prior["mu0"]);

	Matrix2d T0=as< Map<MatrixXd> >(prior["T0"]), Tj, Sigma, Tprec;

	Tj << 1./tau, -rho/sqrt(tau), -rho/sqrt(tau), 1.;
	Tj *= (double(n) * tau)/(1.-pow(rho, 2.));
	
	Tprec=T0+Tj;

	Sigma=Tprec.inverse();

	Vector2d mean_mu=Sigma*(T0*mu0+Tj*mean_eps);

	phi=rnorm2d(mean_mu, Sigma);
	phi.conservativeResize(4);

	/* semi block begins */

 	const double beta=as<double>(prior["beta"]);

  	const VectorXd t=A.col(2), y=A.col(3), 
    		gamma=as< Map<VectorXd> >(prior["gamma"]), 
    		delta=as< Map<VectorXd> >(prior["delta"]), 
    		eta  =as< Map<VectorXd> >(prior["eta"]); 

  	const int p=gamma.size(), q=delta.size();

  	MatrixXd D(n, 2), X=A.block(0, 4, n, p), Z=A.block(0, p+4, n, q);

	D.col(0).setConstant(phi[0]);
	D.col(0) += (X*gamma +Z*delta); //mutmargin
	D.col(1).setConstant(phi[1]);
	D.col(1) += (beta*D.col(0)+X*eta); //muymargin
	D.col(0)= t-D.col(0); //tdev
	D.col(1)= y-D.col(1)-beta*D.col(0);

 
  	Matrix2d SSCP=D.transpose()*D;

	rho=JKB_inverse(corrMH(SSCP(0,0)+tau*SSCP(1,1), sqrt(tau)*SSCP(0,1), n, zeta));

	phi[2]=rho;

	double xinits[4]={0.1, 1., 10., 100.};
	double lower=0.01, upper=1000.; // upper=inf DOES NOT WORK
  	double convex=1.;
  	double xprev, samp;
  	double qcent[2], xcent[2];
  	double roottau;
  	int neval, err=0, maxenv=50, metro=0, nsamp=1;
	
	myldfd[0]=as<double>(prior["alpha0"]);
	myldfd[0] += (n/2.);

	myldfd[1]=D.col(1).transpose()*D.col(1);
    	myldfd[1]/=(2.*(1.-pow(rho, 2.)));
    	myldfd[1]+=as<double>(prior["lambda0"]);	

    	myldfd[2]=D.col(1).transpose()*D.col(0);
    	myldfd[2]*=rho;
    	myldfd[2]/=(1.-pow(rho, 2.));

    	double tmax=(myldfd[2]+sqrt(pow(myldfd[2], 2.)+8.*myldfd[1]*(2.*myldfd[0]-1.)))/(4.*myldfd[1]);
	double delta1=1./((2.*myldfd[0]-1.)/pow(tmax, 2.)-2.*myldfd[1]);
	double t1=max(0.001, tmax+delta1);
    	double t2 = tmax-delta1;
	xprev=tau;	

	//	if(n==1) metro=1;

    	err=arms(&xinits[0], 4, &lower, &upper, &bnivLDF, &myldfd[0], 
	     &convex, maxenv, metro, &xprev, &samp, nsamp, NULL, NULL, 0, &neval);

    	if(err!=0) samp=NA_REAL;
	
    	roottau=samp;
    	tau=pow(roottau, 2.);
	phi[3]=tau;

	return phi;
}

RcppExport SEXP bnivDPM(SEXP arg1, SEXP arg2, SEXP arg3) {
	//3 arguments
	//arg1 for parameters
	//arg2 for data
	//arg3 for parameters

	//data
	List list2(arg2);

	const MatrixXd X=as< Map<MatrixXd> >(list2["X"]),
		Z=as< Map<MatrixXd> >(list2["Z"]);

	const VectorXi tbin=as< Map<MatrixXi> >(list2["tbin"]);

	const VectorXd y=as< Map<MatrixXd> >(list2["y"]);

	const int N=X.rows(), p=X.cols(), q=Z.cols(), r=p+q, s=p+r;

	// parameters

	List list1(arg1),
		beta_info=list1["beta"],
		rho_info=list1["rho"],
		mu_info=list1["mu"],
		tau_info=list1["tau"],
		theta_info=list1["theta"],
		dpm_info=list1["dpm"], //DPM
		alpha_info=dpm_info["alpha"], //alpha random
		alpha_prior; //alpha random

	const int m=as<int>(dpm_info["m"]); //DPM
	//const double alpha=as<double>(dpm_info["alpha"]); //DPM
	const int alpha_fixed=as<int>(alpha_info["fixed"]); // alpha random
  	double alpha=as<double>(alpha_info["init"]); // alpha random

  	if(alpha_fixed==0) alpha_prior=alpha_info["prior"]; // alpha random

	VectorXi C=as< Map<VectorXi> >(dpm_info["C"]), 
	  states=as< Map<VectorXi> >(dpm_info["states"]);

	// prior parameters

	List beta_prior=beta_info["prior"],
		//no prior parameters for rho
	//mu_prior=mu_info["prior"],
	//tau_prior=tau_info["prior"],
	dpm_prior=dpm_info["prior"],	//DPM
	theta_prior=theta_info["prior"];

	const double beta_prior_mean=as<double>(beta_prior["mean"]);
	const double beta_prior_prec=as<double>(beta_prior["prec"]);

	const VectorXd theta_prior_mean=as< Map<VectorXd> >(theta_prior["mean"]);
	const MatrixXd theta_prior_prec=as< Map<MatrixXd> >(theta_prior["prec"]);

	// initialize parameters

	double beta=as<double>(beta_info["init"]);
	double rho_init=as<double>(rho_info["init"]);	//DPM
	double tau_init=as<double>(tau_info["init"]);	//DPM

	Vector2d mu_init=as< Map<VectorXd> >(mu_info["init"]);	//DPM

	VectorXd theta=as< Map<VectorXd> >(theta_info["init"]);

	VectorXd rho(N), tau(N);	//DPM

	/*
	VectorXi C, states; 	//DPM
	C.setConstant(N,0);	//DPM
	states.setConstant(1,N);	//DPM
	*/

	MatrixXd mu(N,2), phi(1,4);	//DPM
	phi(0,0)=mu_init[0];	//DPM
	phi(0,1)=mu_init[1];	//DPM
	phi(0,2)=rho_init;	//DPM
	phi(0,3)=tau_init;	//DPM

	//Gibbs

	List list3(arg3); 

	const int burnin=as<int>(list3["burnin"]), M=as<int>(list3["M"]), 
	  thin=as<int>(list3["thin"]);

	VectorXd t(N); //latent variable

	//prior parameter intermediate values

	double beta_prior_prod=beta_prior_prec*beta_prior_mean;
	VectorXd theta_prior_prod=theta_prior_prec*theta_prior_mean;


	// parameter intermediate values

	Matrix2d Sigma, Tprec, B_inverse;
	B_inverse.setIdentity();

	VectorXd gamma=theta.segment(0,p),
		delta=theta.segment(p,q),
		eta=theta.segment(r,p);

	MatrixXd eps(N,2), D(N,2), theta_cond_var_root(s,s), W(2,s); //DPM
	W.setZero();

	MatrixXd theta_cond_prec(s,s);
	VectorXd theta_cond_prod(s), w(r);

	Vector2d u, R, mu_u;

	VectorXd mutmargin(N), muymargin(N),endpoint(N), numerator(N), denominator(N);
	VectorXd sdtcondy(N);
	double b_scale,beta_prec,beta_prod,beta_cond_var,beta_cond_mean, beta2;
	VectorXd mutcondy(N);
	const double inf=1.0/0.0;

	int i, h,j,k,l,b;

	List GS(M);

	l=-burnin;

	do{
		//populate mu and rho and tau
		//DPM		
		for(i=0; i<N; ++i){
			mu(i,0)=phi(C[i], 0);
			mu(i,1)=phi(C[i],1);
			rho[i]=phi(C[i],2);
			tau[i]=phi(C[i], 3);
			mutmargin[i]=mu(i,0);
			muymargin[i]=mu(i,1);
		}
			

		//generate latents
		mutmargin+=Z*delta+X*gamma; //DPM

		muymargin+=(beta*mutmargin+X*eta); //DPM
	
		
		//DPM
		for(i=0; i<N; ++i) {
			numerator[i]=beta+rho[i]*sqrt(1/tau[i]);
			denominator[i]=pow(beta,2)+2*beta*rho[i]*sqrt(1/tau[i])+(1/tau[i]);
			mutcondy[i]=mutmargin[i]+(numerator[i]/denominator[i])*(y[i]-muymargin[i]);
			sdtcondy[i]=sqrt(1-pow(numerator[i],2)/denominator[i]);
			endpoint[i]=(0-mutcondy[i])/sdtcondy[i];
			if(tbin[i]==0) t[i]=rtnorm(-inf, endpoint[i]);
			else t[i]=rtnorm(endpoint[i], inf);
			t[i]=sdtcondy[i]*t[i]+mutcondy[i];
		}
	
	
		D.col(0)=t-mu.col(0)-X*gamma-Z*delta;	//DPM
		D.col(1)=y-mu.col(1)-X*eta;	//DPM
	
		// sample beta

		beta_prec=0;
		beta_prod=0;
		
		for(i=0; i<N; ++i) {
			double Sigma_det=(1/tau[i])*(1-pow(rho[i],2));	//DPM
			beta_prec +=  pow(t[i], 2.)/Sigma_det;	//DPM
			beta_prod += -t[i]*(sqrt(1/tau[i])*rho[i]*D(i,0)-D(i,1))/Sigma_det; //DPM	
		}

		beta_cond_var=1./(beta_prec+beta_prior_prec);
		beta_cond_mean=beta_cond_var*(beta_prod+beta_prior_prod);
		
		beta=rnorm1d(beta_cond_mean, sqrt(beta_cond_var));		

		B_inverse(1,0)=-beta;

		// sample theta
	
		theta_cond_prec=theta_prior_prec;
		theta_cond_prod=theta_prior_prod;
	
		for(k=0; k<N; ++k) {
			double Sigma_det=(1/tau[k])*(1-pow(rho[k],2));	//DPM
			Tprec(0,0)=(1/tau[k])/Sigma_det; 	//DPM
			Tprec(0,1)=(-rho[k]/sqrt(tau[k]))/Sigma_det;	//DPM
			Tprec(1,0)=(-rho[k]/sqrt(tau[k]))/Sigma_det;	//DPM
			Tprec(1,1)=1/Sigma_det;	//DPM

			mu_u[0]=mu(k,0);	//DPM
			mu_u[1]=mu(k,1);	//DPM
		
			W.block(0,0,1,p)=X.row(k);
			W.block(0,p,1,q)=Z.row(k);
			W.block(1,r,1,p)=X.row(k);
	
			theta_cond_prec += (W.transpose() * Tprec * W);
		
			u[0]=t[k];
			u[1]=y[k];

			R=B_inverse*u-mu_u;
		
			theta_cond_prod += (W.transpose()*Tprec * R);
		}

		theta_cond_var_root=inv_root_chol(theta_cond_prec);
		//theta_cond_var_root=inv_root_svd(theta_cond_prec);
		theta=theta_cond_var_root*(rnormXd(s)+theta_cond_var_root.transpose()*theta_cond_prod);

		gamma=theta.segment(0,p);
		delta=theta.segment(p,q);
		eta=theta.segment(r,p);

		// sample mu, rho, and tau

		for(k=0; k<N; ++k) {
			W.block(0,0,1,p)=X.row(k);
			W.block(0,p,1,q)=Z.row(k);
			W.block(1,r,1,p)=X.row(k);
	
			u[0]=t[k];
			u[1]=y[k];

			eps.row(k) = (B_inverse*u-W*theta).transpose();		//DPM
		}

		/* semi block begins */
    	 	MatrixXd A(N, r+4);

		A.block(0, 0, N, 2)=eps;
		A.col(2)=t;
    		A.col(3)=y;
    		A.block(0, 4, N, p)=X;
    		A.block(0, p+4, N, q)=Z;

    		List psi=List::create(Named("mu0")=as< Map<VectorXd> >(dpm_prior["mu0"]), 
			  Named("T0")=as< Map<MatrixXd> >(dpm_prior["T0"]),   
			  Named("S0")=as< Map<MatrixXd> >(dpm_prior["S0"]), 
			  Named("alpha0")=as<double>(dpm_prior["alpha0"]),
			  Named("lambda0")=as<double>(dpm_prior["lambda0"]), 
			  Named("beta")=beta,   
			  Named("gamma")=gamma, 
			  Named("delta")=delta, 
			  Named("eta")=eta);    

		if(alpha_fixed==0) 
     			 alpha=bniv_alpha(states.size(), N, alpha, alpha_prior);

    		List dpm_step=neal8(A, C, phi, states, m, alpha, psi, 
			&bnivF, &bnivG0, &bnivP0);
    		/* semi block end */

		C=as< Map<VectorXi> >(dpm_step["C"]);
		phi=as< Map<MatrixXd> >(dpm_step["phi"]);
		states=as< Map<VectorXi> >(dpm_step["states"]);
		
		if(l>=0 && l%thin==0) {
			h=(l/thin);

			GS[h]=List::create(Named("beta")=beta, 
					   Named("theta")=theta,
					   Named("C")=C, 
					   Named("phi")=phi, 
					   Named("states")=states, 
					   Named("alpha")=alpha);
		}
	
		++l;	
	} while(l<=(M-1)*thin);

	return wrap(GS);	
}
