#include "bivpack.h"

RNGScope scope;

// U(a, b)
double runif1d(double a, double b) {
  return as<double>(Rcpp::runif(1, a, b));
}

// iid U(a, b)
VectorXd runifXd(int size, double a, double b) {
  return as<VectorXd>(Rcpp::runif(size, a, b));
}

// N(mu, sd^2)
double rnorm1d(double mu, double sd) {
  return as<double>(rnorm(1, mu, sd));
}

// iid N(mu, sd^2)
Vector2d rnorm2d(double mu, double sd) {
  return as<VectorXd>(rnorm(2, mu, sd));
}

VectorXd rnormXd(int size, double mu, double sd) {
  return as<VectorXd>(rnorm(size, mu, sd));
}

// MVN(mu, Sigma)
VectorXd rnormXd(VectorXd mu, MatrixXd Sigma) {
  MatrixXd L=Sigma.llt().matrixL();

  return L*rnormXd(mu.size())+mu;
}

Vector2d rnorm2d(Vector2d mu, Matrix2d Sigma) {
  Matrix2d L=Sigma.llt().matrixL();

  return L*rnorm2d()+mu;
}

// N(0, 1)I(a, b)
double rtnorm(double a, double b) {
  if(a<4. && b>-4.) {
    NumericVector c=NumericVector::create(a), d=NumericVector::create(b),
      Pc=pnorm(c, 0., 1.), Pd=pnorm(d, 0., 1.), U=Rcpp::runif(1, 0., 1.),
      Q=qnorm(Pc+U*(Pd-Pc), 0., 1.);
    
    return as<double>(Q);
  }
  else {
    int sign=1;
    
    if(a<0.) {
      double c=a;
      a=-b;
      b=-c;
      sign=-1;
    }

    a=abs(a);

    double a2=pow(a, 2.), b2=pow(b, 2.), x;

    do { x=sqrt(a2-2.*log(1.-runif1d()*(1.-exp(-(b2-a2)/2.) ) )); } 
    while(runif1d() > (a/x));

    return sign*x;
  }
}

// N(mu, sd^2)I(a, b)
double rtnorm(double a, double b, double mu, double sd) {
  return mu+sd*rtnorm((a-mu)/sd, (b-mu)/sd);
}

// Gam(shape, rate) where mean=shape/rate
double rgamma1d(double shape, double rate) {
  return as<double>(rgamma(1, shape, 1./rate));
}

// iid Gam(shape, rate) where mean=shape/rate
VectorXd rgammaXd(int size, double shape, double rate) {
  return as<VectorXd>(rgamma(size, shape, 1./rate));
}

// Chisq(df)
double rchisq1d(double df) {
  return as<double>(rgamma(1, df/2., 2.));
}

// iid Chisq(df)
VectorXd rchisqXd(int size, double df) {
  return as<VectorXd>(rgamma(size, df/2., 2.));
}

// W(Omega, nu) where mean=nu*Omega
MatrixXd rwishart(MatrixXd Omega, double nu) {
  const int p=Omega.rows();

  MatrixXd T(p, p); T.setZero();

  // generates same stream as rwish() 
  for(int i=0; i<p; ++i) T(i, i) = sqrt(rchisq1d(nu-i));

  for(int j=0; j<(p-1); ++j) for(int i=j+1; i<p; ++i) T(i, j) = rnorm1d();

  MatrixXd A=T*T.transpose();
  MatrixXd L=Omega.llt().matrixL();

  return L*A*L.transpose();
}

MatrixXd rwishart(double nu, MatrixXd Omega) {
  return rwishart(Omega, nu);
}

// IW(Omega, nu)
MatrixXd rinvwishart(MatrixXd Omega, double nu) {
  return rwishart(Omega, nu).inverse();
}

MatrixXd rinvwishart(double nu, MatrixXd Omega) {
  return rinvwishart(Omega, nu);
}

Matrix2d inv_root_chol(Matrix2d &A) {
// this is a little faster, but ...
// you need to remember to use transpose appropriately since it is not symmetric
  Matrix2d L=A.llt().matrixL(), R=A.llt().solve(L);
  //  L=A.ldlt().matrixL(), 
  //  D=A.ldlt().vectorD().asDiagonal(), 
  //  R=A.ldlt().solve(L*(D.array().sqrt().matrix()));

  return R;
}
  
MatrixXd inv_root_chol(MatrixXd &A) {
// this is a little faster, but ...
// you need to remember to use transpose appropriately since it is not symmetric
  MatrixXd L=A.llt().matrixL(), R=A.llt().solve(L);
  return R;
}

Matrix2d inv_root_svd(Matrix2d &A) {
  // this is a little slower, but it is symmetric
  JacobiSVD<Matrix2d> SVD(A, ComputeThinU | ComputeThinV);
  Array2d sqrt_D=SVD.singularValues().array().sqrt();
  Vector2d inv_sqrt_D=(1./sqrt_D).matrix();
  return SVD.matrixV() * inv_sqrt_D.asDiagonal() * SVD.matrixU().transpose();
}
  
MatrixXd inv_root_svd(MatrixXd &A) {
  // this is a little slower, but it is symmetric
  JacobiSVD<MatrixXd> SVD(A, ComputeThinU | ComputeThinV);
  ArrayXd sqrt_D=SVD.singularValues().array().sqrt();
  VectorXd inv_sqrt_D=(1./sqrt_D).matrix();
  return SVD.matrixV() * inv_sqrt_D.asDiagonal() * SVD.matrixU().transpose();
}

double JKB(double rho) {return rho/sqrt(1.-pow(rho, 2.));}

double JKB_inverse(double zeta) {return zeta/sqrt(1.+pow(zeta, 2.));}

/*
double corrMH(double a, double b, int N, double x) {
  double x2=pow(x, 2.), y=as<double>(rcauchy(1, x, pow(a, -0.5))), y2=pow(y, 2.),
    p = pow((1.+y2)/(1.+x2), 0.5*N-1.5) 
    * exp(b*(y*sqrt(1.+y2)-x*sqrt(1.+x2))-0.5*a*(y2-x2));
         
  return runif1d()<p ? y : x;
}
*/

// returns zeta rather than rho; also prev=prev.zeta
double corrMH(double a, double b, int N, double prev, int m) {
  for(int i=0; i<m; ++i) {
    double prev2=pow(prev, 2.), 
      next=as<double>(rcauchy(1, prev, pow(a, -0.5))), next2=pow(next, 2.),
      p = pow((1.+next2)/(1.+prev2), 0.5*N-1.5) 
      * exp(b*(next*sqrt(1.+next2)-prev*sqrt(1.+prev2))-0.5*a*(next2-prev2));
         
    if(runif1d()<p) prev=next;  
  }

  return prev;
}

// discretely samples rho; preferable to corrMH
double corrMN(double a, double b, int N, int n) {
  const int m=2*n-1;

  const double step=1./n;

  VectorXd grid(m), prob(m);

  for(int i=0; i<m; ++i) {
    grid[i]=(1.-n+i)*step;

    double c=1.-pow(grid[i], 2.);

    prob[i]=1./(pow(c, N/2.)*exp((0.5*a-b*grid[i])/c));
  }
  
  return grid[rmultinom1(prob)];
}

/*
double corrBetaMH(double a, double b, int N, double prev_rho, int m) {
  for(int i=0; i<m; ++i) {
    double prev_rho2=pow(prev_rho, 2.), prev_zeta=(prev_rho+1.)/2., 
      p=NaN, next_rho;

    while(p!=p) {
      double next_zeta=R::rbeta(prev_zeta, 1.-prev_zeta);

      next_rho=2.*next_zeta-1.;

      double next_rho2=pow(next_rho, 2.);

      p=pow((1.+next_rho)/(1.-next_rho), prev_rho/2.)*
        pow((1.-prev_rho)/(1.+prev_rho), next_rho/2.)*
	pow((1.-prev_rho2)/(1.-next_rho2), (N+1.)/2.)*
	exp(((0.5*a-b*prev_rho)/(1.-prev_rho2))
           -((0.5*a-b*next_rho)/(1.-next_rho2)));
    }  

    if(runif1d()<p) prev_rho=next_rho;
  }

  return prev_rho;
}
*/

Vector2d rbvtruncnorm(Vector2d value, int quadrant, 
	       double rho, int burnin, int M, int thin) {
  Matrix2d L;

  L << 1., 0., rho, 1.;

  Vector2d x, mu(value);

  double sdy=sqrt(1.-pow(rho, 2.));

  mu[1]=mu[1]-mu[0]*rho;

    if (quadrant==1) {
      x[0]=1.; x[1]=1.;

      x[1]=x[1]-rho*x[0];

      for(int g=-burnin; g<=(M-1)*thin; ++g) {
	if(rho>0) {
          x[0]=mu[0]+rtnorm(max(0., -x[1]/rho)-mu[0], Inf);
          x[1]=mu[1]+sdy*rtnorm((-rho*x[0]-mu[1])/sdy, Inf);
        }
	else if(rho<0) {
          x[0]=mu[0]+rtnorm(-mu[0], -x[1]/rho-mu[0]);
          x[1]=mu[1]+sdy*rtnorm((-rho*x[0]-mu[1])/sdy, Inf);
        }
	else {
          x[0]=mu[0]+rtnorm(-mu[0], Inf);
          x[1]=mu[1]+sdy*rtnorm(-mu[1]/sdy, Inf);
        }
      }
    }  
    else if (quadrant==2) {
      x[0]=-1.; x[1]=1.;

      x[1]=x[1]-rho*x[0];

      for(int g=-burnin; g<=(M-1)*thin; ++g) {
	if(rho>0) {
          x[0]=mu[0]+rtnorm(-x[1]/rho-mu[0], -mu[0]);
          x[1]=mu[1]+sdy*rtnorm((-rho*x[0]-mu[1])/sdy, Inf);
        }
	else if(rho<0) {
          x[0]=mu[0]+rtnorm(-Inf, min(0., -x[1]/rho)-mu[0]);
          x[1]=mu[1]+sdy*rtnorm((-rho*x[0]-mu[1])/sdy, Inf);
        }
	else {
          x[0]=mu[0]+rtnorm(-Inf, -mu[0]);
	  x[1]=mu[1]+sdy*rtnorm(-mu[1]/sdy, Inf);
        }
      }    
    }
    else if (quadrant==3) {
      x[0]=-1.; x[1]=-1.;

      x[1]=x[1]-rho*x[0];

      for(int g=-burnin; g<=(M-1)*thin; ++g) {
	if(rho>0) {
          x[0]=mu[0]+rtnorm(-Inf, min(0., -x[1]/rho)-mu[0]);
          x[1]=mu[1]+sdy*rtnorm(-Inf, (-rho*x[0]-mu[1])/sdy);
        }
	else if(rho<0) {
          x[0]=mu[0]+rtnorm(-x[1]/rho-mu[0], -mu[0]);
          x[1]=mu[1]+sdy*rtnorm(-Inf, (-rho*x[0]-mu[1])/sdy);
        }
	else {
          x[0]=mu[0]+rtnorm(-Inf, -mu[0]);
	  x[1]=mu[1]+sdy*rtnorm(-Inf, -mu[1]/sdy);
        }
      }        
    }
    else {
      x[0]=1.; x[1]=-1.;

      x[1]=x[1]-rho*x[0];

      for(int g=-burnin; g<=(M-1)*thin; ++g) {
	if(rho>0) {
          x[0]=mu[0]+rtnorm(-mu[0], -x[1]/rho-mu[0]);
          x[1]=mu[1]+sdy*rtnorm(-Inf, (-rho*x[0]-mu[1])/sdy);
        }
	else if(rho<0) {
          x[0]=mu[0]+rtnorm(max(0., -x[1]/rho)-mu[0], Inf);
          x[1]=mu[1]+sdy*rtnorm(-Inf, (-rho*x[0]-mu[1])/sdy);
        }
	else {
          x[0]=mu[0]+rtnorm(-mu[0], Inf);
	  x[1]=mu[1]+sdy*rtnorm(-Inf, -mu[1]/sdy);
        }
      }            
    }

    return L*x;
}

/*
MatrixXd erase(MatrixXd A, int b) {
  // erase the b-th row of a matrix
  const int n=A.rows()-1;

  MatrixXd C=A;

  C.conservativeResize(n, A.cols());

  for(int i=b; i<n; ++i) C.row(i)=A.row(i+1);
  
  return C;
}

VectorXi freq(VectorXi C) {
  const int n=C.size(), k=C.maxCoeff()+1; // 0, ..., max

  VectorXi a;

  a.setZero(k);

  for(int i=0; i<k; ++i) {
    for(int j=0; j<n; ++j) if(C[j]==i) a[i]=a[i]+1;
  }

  return(a);
}
*/

RowVectorXd mean(MatrixXd A) {
  // mean of A applied over the columns: row vector of means returned
  const int n=A.rows();

  RowVectorXd x;

  x.setConstant(n, 1./n);

  return x*A;
}

/*
VectorXd mean(MatrixXd A) {
  // mean of A applied over the columns: vector of means returned
  const int n=A.rows();

  VectorXd x;

  x.setConstant(n, 1./n);

  return A.transpose()*x;
}
*/

MatrixXd SSCP(MatrixXd A) {
  // return a matrix of the centered SSCP
  RowVectorXd mean_A=mean(A);
  const int n=A.rows();

  for(int i=0; i<n; ++i) A.row(i) -= mean_A;

  return A.transpose()*A;
}

int rmultinom1(VectorXd prob) {
  // return a draw from Multinomial(1, prob)
  // e.g. c(1:h) %*% rmultinom(1, 1, prob)
  const int n=prob.size(), m=n-1;

  prob /= prob.sum();

  int r=m;

  double cum=prob[0], u=runif1d();

  for(int i=1; i<n && r==m; ++i) {
    if(u<cum) r=i-1;

    cum += prob[i];
  }

  return r;
}

