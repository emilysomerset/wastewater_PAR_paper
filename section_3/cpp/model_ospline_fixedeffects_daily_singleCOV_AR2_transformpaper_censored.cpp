#include <TMB.hpp>                                // Links in the TMB libraries
//#include <fenv.h>

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_VECTOR(y); // response variable
  // DATA_VECTOR(t); // need this to calculate delta
  DATA_VECTOR(denom); // response variable
  DATA_SPARSE_MATRIX(X); // Design matrix (for fixed effects)
  DATA_SPARSE_MATRIX(B); // Design matrix (for random effects)
  DATA_SPARSE_MATRIX(P); // Penalty matrix
  DATA_SPARSE_MATRIX(daily);
  DATA_SPARSE_MATRIX(obs); //observed data indexes
  DATA_IVECTOR(stationsizes);
  DATA_IVECTOR(n1);
  DATA_IVECTOR(y_ind_obs);
  DATA_IVECTOR(y_ind_cens);
  
  //Prior stuff
  DATA_SCALAR(betaprec); // beta ~iid N(0,1/betaprec)
  DATA_SCALAR(logPdet); // Determinant of (fixed) penalty matrix
  DATA_SCALAR(u1); // pc prior, u1 param
  DATA_SCALAR(alpha1); // pc prior, alpha1 param
  DATA_SCALAR(u2); // pc prior, u2 param
  DATA_SCALAR(alpha2);
  DATA_SCALAR(lambda_phi); //psi = 1/range
  DATA_SCALAR(lambda_tau);
  DATA_SCALAR(lambda_cov);
  
  // Parameter
  PARAMETER_VECTOR(W); // W = c(U,beta,Z), eta = B * U + X * beta + daily * Z
  
  int d = P.cols(); // Number of B-Spline coefficients
  int betadim = X.cols();
  int ndays = daily.cols();
  int nstation = stationsizes.size();
  int ncens = y_ind_cens.size();
  int nfull = obs.cols();
  
  vector<Type> U(d);
  vector<Type> beta(betadim);
  vector<Type> Z(ndays);
  matrix<Type> alpha(2,nfull); // alpha
  
  
  for (int i=0;i<d;i++) U(i) = W(i);
  for (int i=0;i<betadim;i++) beta(i) = W(i+d);
  for (int i=0;i<ndays;i++) Z(i) = W(i+d + betadim);
  for (int i=0;i<nfull;i++) alpha(0,i) = W(i + d + ndays + betadim);
  for (int i=0;i<nfull;i++) alpha(1,i) = W(i + d + ndays + betadim + nfull);
  
  PARAMETER(theta1); // Associated with u1, alpha1, and the O-splines
  PARAMETER(theta2); // Associated with u2, alpha2, and the daily random effect
  PARAMETER(cov_log);
  PARAMETER(theta3); // Associated with u_phi, alpha_phi, associated with Matern phi
  PARAMETER(theta4); //Associated with u_tau, alpha_tau, associated with Matern tau
  
  // Transformations
  vector<Type> alpha_0 = alpha.row(0);
  Type phi = exp((theta3-theta4/2)/sqrt(2));
  Type tau = exp((theta3+theta4/2)/sqrt(2));
  Type c = pow(12,0.5)*pow(phi,-1);
  Type sigmastat = tau*sqrt(2)*pow(c,1.5);
  Type delta = 1;
  matrix<Type> obsdaily = obs * daily;
  vector<Type> eta = X * beta + B * U + obsdaily * Z + obs * alpha_0;
  Type cov = exp(cov_log);
  vector<Type> scale = exp(eta) * pow(cov,2) * denom;
  Type shape = 1/pow(cov,2);
  
  
  // Log likelihood
  Type ll = 0;
  ll += sum(dgamma(y(y_ind_obs), shape, scale(y_ind_obs), TRUE));
  for (int i=0;i<ncens;i++) {
    ll += log(pgamma(y(y_ind_cens(i)), shape, scale(y_ind_cens(i))));
  }
  
  // Log prior on W
  Type lpW = 0;
  
  // Cross product
  vector<Type> PU = P*U;
  Type UPU = (U * PU).sum();
  lpW += -0.5 * exp(theta1) * UPU; // U part
  Type bb = (beta * beta).sum();
  lpW += -0.5 * betaprec * bb; // Beta part
  Type zz = (Z * Z).sum();
  lpW += -0.5 * exp(theta2) * zz; // Z part
  
  // Log determinant
  Type logdet1 = d * theta1 + logPdet;
  lpW += 0.5 * logdet1; // P part
  
  Type logdet2 = ndays * theta2;
  lpW += 0.5 * logdet2; // Z part
  
  // Log prior on alpha
  Type lpA = 0;
  
  //prior for the first column of alpha
  matrix<Type> T(2,2);
  matrix<Type> Sigma(2,2);
  matrix<Type> Sigma_inv(2,2);
  matrix<Type> Sigma_init(2,2);
  matrix<Type> Sigma_init_inv(2,2);
  
  T(0,0) = exp(-c*delta)*(1+c*delta);
  T(0,1) = exp(-c*delta)*delta;
  T(1,0) = -exp(-c*delta)*pow(c,2)*delta;
  T(1,1) = exp(-c*delta)*(1-c*delta);
  
  matrix<Type> alpha_mu = T * alpha; // T*alpha
  
  Sigma_init(0,0) = pow(sigmastat,2)/(4*pow(c,3));
  Sigma_init(0,1) = 0;
  Sigma_init(1,0) = 0;
  Sigma_init(1,1) = pow(sigmastat,2)/(4*c);
  
  Sigma_init_inv(0,0) = Sigma_init(1,1);
  Sigma_init_inv(1,1) = Sigma_init(0,0);
  Sigma_init_inv(1,0) = 0;
  Sigma_init_inv(0,1) = 0;
  
  Sigma(0,0) = pow(sigmastat,2)/(4*pow(c,3))*(1-exp(-2*c*delta)*(-2*c*delta*(-delta*c-1)+1));
  Sigma(0,1) = pow(sigmastat,2)*0.5*pow(delta,2)*exp(-2*delta*c);
  Sigma(1,0) = pow(sigmastat,2)*0.5*pow(delta,2)*exp(-2*delta*c);
  Sigma(1,1) = pow(sigmastat,2)/(4*c)*(1-exp(-2*delta*c)*(-2*delta*c*(-delta*c+1)+1));
  
  Sigma_inv(0,0) = Sigma(1,1);
  Sigma_inv(1,1) = Sigma(0,0);
  Sigma_inv(1,0) = -Sigma(1,0);
  Sigma_inv(0,1) = -Sigma(0,1);
  
  Type det_Sigma_init = Sigma_init(0,0)*Sigma_init(1,1);
  Type det_Sigma = Sigma(0,0)*Sigma(1,1)-Sigma(1,0)*Sigma(0,1);
  
  for (int j=0; j<nstation; j++){
    vector<Type> PU = Sigma_init_inv*alpha.col(n1(j));
    vector<Type> alpha_col0 = alpha.col(n1(j));
    Type UPU = (alpha_col0 * PU).sum();
    lpA += -0.5*log(det_Sigma_init)-0.5*pow(det_Sigma_init,-1)*UPU;
    Rcout << "lpAInitial: " << lpA << "\n";
    
    for (int i=(n1(j)+1); i < n1(j+1); i++){
      vector<Type> resid_alpha = alpha.col(i)-alpha_mu.col(i-1);
      vector<Type> PU = Sigma_inv*resid_alpha;
      Type UPU = (resid_alpha * PU).sum();
      lpA += -0.5*log(det_Sigma)-0.5*pow(det_Sigma,-1)*UPU;
    }}
  
  
  // log prior for theta3 and theta4
  Type lplSlp = log(lambda_phi) + log(lambda_tau)-lambda_tau*tau  - lambda_phi*phi - log(2) + log(phi) + log(tau);
  
  
  REPORT(lpW);
  
  // Log prior for theta
  Type lpT = 0;
  Type phi1 = -log(alpha1) / u1;
  lpT += log(0.5 * phi1) - phi1*exp(-0.5*theta1) - 0.5*theta1;
  Type phi2 = -log(alpha2) / u2;
  lpT += log(0.5 * phi2) - phi2*exp(-0.5*theta2) - 0.5*theta2;
  lpT += log(lambda_cov) - lambda_cov*exp(cov_log) + cov_log;
  REPORT(lpT);
  
  // Final result!
  Type logpost = -1 * (lpW + lpT + ll + lplSlp+ lpA);
  REPORT(logpost);
  
  return logpost;
}