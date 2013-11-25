#include "bivpack.h"

List neal8(const MatrixXd Y, VectorXi C, MatrixXd phi, VectorXi states, 
	   const int m, const double alpha, const List prior,
	   double (*F)(VectorXd y, VectorXd phi, List prior),
	   VectorXd (*G0)(List prior),
	   VectorXd (*P0)(MatrixXd y, List prior, VectorXd phi)) {

  const int N=Y.rows(), M=Y.cols(), p=phi.cols();  
  int k=states.size(); // number of states

BEGIN_RCPP
  // superficially reality check C/states/phi/m
  if(states.sum()!=N || C.size()!=N || phi.rows()!=k || m<1) {
    stringstream stream;
    string string;
 
    stream << "States sum:" << states.sum() << '\n'
	   << "C size:" << C.size() << '\n'
	   << "States size:" << k << '\n'
           << "phi rows:" << phi.rows() << '\n'
           << "m:" << m << '\n';
    stream >> string;

    throw range_error(string); 
  }

#ifdef DEBUG_NEAL8
  VectorXi S(N); 
  List P(N), Phi(N);
#endif

  if(alpha>0.) {
    const double alpha_m=alpha/m;

    for (int i=0; i<N; ++i) {
      int c_i = C[i];

      int h = k+m, singleton=-1;

      phi.conservativeResize(h, p);

      if(states[c_i]==1) {
	singleton = k-1;
	h--;

	if(c_i<singleton) {
	  for (int j=0; j<N; ++j) if(c_i<C[j]) C[j]=C[j]-1;

	  // phi=erase(phi, c_i); // erase() replaced with 2 lines below
	  phi.row(k)=phi.row(c_i);

	  for (int c=c_i; c<singleton; ++c) {
	    states[c]=states[c+1];
	    phi.row(c)=phi.row(c+1); // now erase() not needed
	  }

	  c_i = singleton; // bug-fix:  this line added
	  C[i] = singleton;
	  states[singleton]=1;
	  phi.row(singleton)=phi.row(k); // now erase() not needed
	}
      }

#ifdef DEBUG_NEAL8
      S[i]=singleton;
#endif

      VectorXd prob(h);

      for(int c=0; c<h; ++c) {
	if((c < (k-1)) || (singleton==(-1) && c==(k-1))) {
	  int n_c = (c==c_i) ? states[c_i]-1 : states[c];

	  prob[c] = n_c* (*F) (Y.row(i), phi.row(c), prior);
	}
	else{
	  if(singleton!=c) phi.row(c) = (*G0) (prior);

	  prob[c] = alpha_m* (*F) (Y.row(i), phi.row(c), prior);
	}
      }

#ifdef DEBUG_NEAL8
      P[i]=List::create(Named("prob")=prob/prob.sum());
#endif

      c_i = rmultinom1(prob);

      if(singleton != -1) {
	if (singleton <= c_i) phi.row(singleton)=phi.row(c_i);
	else {
	  C[i]=c_i; 
	  k--;

	  states.conservativeResize(k);
	  states[c_i]=states[c_i]+1;
	}
      }
      else if((k-1)<c_i) {
	states[C[i]]=states[C[i]]-1;
	C[i] = k;
	phi.row(k) = phi.row(c_i);
	c_i=k;
	k++;

	states.conservativeResize(k);
	states[c_i]=1;
      }
      else if(C[i] != c_i) {
	states[C[i]]=states[C[i]]-1;
	states[c_i]=states[c_i]+1;
      
	C[i] = c_i;
      }

#ifdef DEBUG_NEAL8
      Phi[i]=List::create(Named("phi")=phi.block(0, 0, h, p));
#endif

      phi.conservativeResize(k, p);
    }
  }

#ifdef DEBUG_NEAL8
  List B(k);
#endif

  for(int c=0; c<k; c++) {
    int k_c=states[c];

    MatrixXd A(k_c, M);
  
    for(int i=0, j=0; i<N && j<k_c; ++i) if(C[i] == c) {
	A.row(j)=Y.row(i);
	j++;
      }
  
    phi.row(c) = (*P0) (A, prior, phi.row(c)); 

#ifdef DEBUG_NEAL8
    B[c]=List::create(Named("A")=A);
#endif
  }
  
  //return List::create(C, phi, states);

#ifdef DEBUG_NEAL8
  return (alpha>0.) 
    ? List::create(Named("C")=C, Named("phi")=phi, Named("states")=states,
      Named("S")=S, Named("P")=P, Named("Phi")=Phi, Named("B")=B)
    : List::create(Named("C")=C, Named("phi")=phi, Named("states")=states,
      Named("B")=B);
#else
  return List::create(Named("C")=C, Named("phi")=phi, Named("states")=states);
#endif
END_RCPP
}

