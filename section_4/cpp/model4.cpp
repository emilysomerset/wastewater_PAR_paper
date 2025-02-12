#include <TMB.hpp>                                // Links in the TMB libraries
//#include <fenv.h>

template<class Type>
Type objective_function<Type>::operator() ()
{
  //cases addition
  DATA_VECTOR(case_counts); //weekly case counts
  DATA_VECTOR(ratio); //same size as case counts. 
  DATA_INTEGER(lag);
  DATA_SCALAR(obs_start_case);
  DATA_SCALAR(u1); // pc prior, u1 param
  DATA_SCALAR(alpha1); // pc prior, alpha1 param
  DATA_SCALAR(mean_z0);
  DATA_SCALAR(sd_z0);
  
  
  // Parameter
  PARAMETER_VECTOR(W); 
  

  int ncounts = case_counts.size();
  
  vector<Type> logit_p(ncounts-lag);
  vector<Type> Z(ncounts-lag);
  
  
  for (int i=0;i<(ncounts-lag);i++) logit_p(i) = W(i);
  for (int i=0;i<(ncounts-lag);i++) Z(i) = W(i+ncounts-lag);
  
  PARAMETER(theta_p);

  
  // Transformations for counts
  
  vector<Type> p(ncounts - lag);
  vector<Type> pZ(ncounts-lag);
  vector<Type> ppZ(ncounts-lag);

  for (int i=0; i<(ncounts-lag);i++){
    p(i) = exp(logit_p(i))/(1+exp(logit_p(i)));
    pZ(i) = p(i)*Z(i); //pZ(0) = pZ1 = p1*Z1
    ppZ(i) = pow(p(i)*(1-p(i))*Z(i),0.5);
  }
  
  
  // Log likelihood
  Type ll = 0;
  for (int i=0; i<(ncounts-lag);i++){
  ll += dnorm(case_counts[i+lag], pZ[i], ppZ[i],TRUE);} // lag =0; Y(0)= Y1, pZ(0) = pZ1 
  Rcout << "ll: " << ll << "\n";

  // Log prior on W
  Type lpW = 0;
  
  lpW += -0.5 * 1 * pow(logit_p(0),2); // logit(p1)~N(0,1) 
  for (int i=1; i<(ncounts-lag);i++){
    lpW += -0.5 * exp(theta_p) * pow((logit_p(i)- logit_p(i-1)),2); //logit(p2)~N(logit(p1),exp(-theta_p))
  }
  lpW += 0.5 * (ncounts-lag-1)  * theta_p;  // Log determinant
  
  lpW += -0.5 * pow(sd_z0,-2) * pow((Z(0)-mean_z0),2);
  for (int i=1; i<(ncounts-lag);i++){
    lpW += -0.5 * pow((Z(i) - ratio(i) * Z(i-1)),2) * pow(ratio(i) * Z(i-1),-1); //Z1 ~ N(lambda1, lambda1)
    lpW += -0.5*log(ratio(i) * Z(i-1));  // Log determinant
  }
  Rcout << "lpW: " << lpW << "\n";
  REPORT(lpW);
  
  // Log prior for theta
  Type lpT = 0;
  Type phi1 = -log(alpha1) / u1;
  lpT += log(0.5 * phi1) - phi1*exp(-0.5*theta_p) - 0.5*theta_p; // theta_p ~ Exp(-log(alpha1)/u1)
  REPORT(lpT);
  Rcout << "lpT: " << lpT << "\n";
  
  // Final result!
  Type logpost = -1 * (lpW + lpT + ll);
  Rcout << "logpost: " << logpost << "\n";
  REPORT(logpost);
  
  return logpost;
}