// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>

// [[Rcpp::export]]
arma::colvec ldot(arma::colvec tau, double alpha,
	       double lambda, arma::mat D,
	       double a, double b)
{


  int i;
  
  arma::mat E = exp( - D/lambda);
  for(i = 0; i < E.n_rows; ++i)
    {
      E(i,i) = E(i,i) + 1e-5;
    }

  arma::mat E_inv = inv(E);
  arma::mat F = 1.0/(lambda * lambda) * D % E;
  arma::mat M  = E_inv * (-F) * E_inv;
  arma::mat G  = -2/(lambda * lambda * lambda) * (D % E) + 1/(lambda * lambda) * (D % F);
  arma::mat L = M * F + E_inv * G;
  
  arma::mat N = M * (-F) * E_inv + E_inv * (-G) * E_inv + E_inv * (-F) * M;
  arma::mat res1 = -0.5 * trace(E_inv * F) - 0.5 * alpha * trans(tau) * M * tau - b + (a - 1) / (lambda * lambda);
  arma::mat res2 = -0.5 * trace(L) - 0.5 * alpha *  trans(tau) * N * tau - (a - 1) / (lambda * lambda);
  arma::colvec res(2);
  res(0) = res1(0,0);
  res(1) = res2(0,0);

  return(res);
}


// [[Rcpp::export]]
double j_double_prime_new(double tau, double tau_hat,
         double varsigma, double kappa, double xi_hat,
	       arma::colvec eps)
{
  double res;
  arma::colvec h = 1 + kappa * eps * (xi_hat + tau);
  if(any(h < 0))
    {
    res = -arma::datum::inf; 
    return(res);
    }

//  arma::colvec f_1 = (xi_hat + tau + 1)/(xi_hat + tau) * log(h); // Unused
  arma::colvec f_2 =  exp(-1/(xi_hat + tau) * log(h));
  
//  arma::colvec f_1_dot = -log(h)/pow(xi_hat + tau,2) + (xi_hat + tau + 1)/(xi_hat + tau) * pow(h,-1) % eps * kappa; // Unused
  arma::colvec f_2_dot = f_2 % ( log(h)/pow(xi_hat + tau,2) - pow(h,-1) % eps * kappa / (xi_hat + tau) );
  arma::colvec g_1_dot = -2 * pow(xi_hat + tau,-3) * log(h) + pow(xi_hat + tau,-2) * pow(h,-1) % eps * kappa;
  arma::colvec g_2_dot = -pow(h,-1) % eps *  pow(xi_hat + tau,-2) * kappa - (xi_hat + tau + 1)/(xi_hat + tau) * pow(h,-2) % pow(eps,2) * pow(kappa,2);
  arma::colvec g_3_dot_1 = f_2_dot % (log(h) * pow(xi_hat + tau,-2) );
  arma::colvec g_3_dot_2 = f_2 % ( -2 * log(h) * pow(xi_hat + tau,-3) + pow(h,-1) %  eps * kappa * pow(xi_hat + tau,-2) );
  arma::colvec g_3_dot = g_3_dot_1 + g_3_dot_2;
  arma::colvec g_4_dot_1 = f_2_dot % (pow(h,-1) % eps * kappa * pow(xi_hat + tau,-1));
  arma::colvec g_4_dot_2 = kappa * -f_2 % eps % ( pow(h,-1) * pow(xi_hat + tau,-2) + pow(h,-2) % eps * kappa * pow(xi_hat + tau,-1) );
  arma::colvec g_4_dot = g_4_dot_1 + g_4_dot_2;
  res = sum(g_1_dot) - sum(g_2_dot) - sum(g_3_dot) + sum(g_4_dot) - 1/varsigma;
return(res);
}

// [[Rcpp::export]]
double j_prime_new(double tau, double tau_hat,
         double varsigma, double kappa, double xi_hat,
         arma::colvec eps)
{
  double res;
  arma::colvec h = 1 + kappa * eps * (xi_hat + tau);
  if(any(h < 0))
    {
    res = -arma::datum::inf; 
    return(res);
    }

//  arma::colvec f_1 = (xi_hat + tau + 1)/(xi_hat + tau) * log(h); // Unused
  arma::colvec f_1_dot = -log(h)/pow(xi_hat + tau,2) + (xi_hat + tau + 1)/(xi_hat + tau) * pow(h,-1) % eps * kappa;
  arma::colvec f_2 =  exp(-log(h)/(xi_hat + tau));
  arma::colvec f_2_dot = f_2 % ( log(h)*pow(xi_hat + tau,-2) - pow(h,-1) % eps * kappa / (xi_hat + tau) );  
  res = -sum(f_1_dot) - sum(f_2_dot) - (tau - tau_hat)/varsigma;

return(res);
}


// [[Rcpp::export]]
arma::colvec gev_like_new(arma::colvec Y, double mu,
         double kappa, double xi)
{
  arma::colvec h;
  arma::colvec L;
  if (fabs(xi)< 1e-8){
    h = exp(-kappa*(Y-mu)); 
    L = log(kappa) + log(h) - h;
    } else {
    h = 1 + xi * kappa * (Y - mu); // No test to see if h is negative here. Include?
    L = log(kappa) - (xi+1)/xi * log(h) - pow(h,-1/xi);
    }
  return(L);
}

// [[Rcpp::export]]
double g_prime_new(double tau, double tau_hat,
         double varsigma, double xi, double kappa_hat,
         arma::colvec eps)
{
  double res;
  arma::colvec h;
  arma::colvec L;
  
  if (fabs(xi)< 1e-8){
    h = exp(-(kappa_hat + tau)*eps);
    L = 1/(kappa_hat + tau) + eps % (h-1);
  } else {
    h = 1 + xi * (kappa_hat + tau) * eps;
    if(any(h < 0))
      {
      res = -arma::datum::inf; 
      return(res);
      }
    L = 1/(kappa_hat + tau) - (xi + 1) * eps % pow(h,-1) + eps % pow(h,-1/xi - 1);
  }
  res = sum(L) - (tau - tau_hat)/varsigma;
  return(res);
}

// [[Rcpp::export]]
double g_double_prime_new(double tau, double tau_hat,
         double varsigma, double xi, double kappa_hat,
         arma::colvec eps)
{
  double res;
  arma::colvec h;
  arma::colvec L;
  
  if (fabs(xi)< 1e-8){
    h = exp(-(kappa_hat + tau)*eps);
    L = -pow(kappa_hat + tau,-2) - pow(eps,2) % h;
  } else {
  h = 1 + xi * (kappa_hat + tau) * eps;
  if(any(h < 0))
    {
    res = -arma::datum::inf; 
    return(res);
    }
  L = -pow(kappa_hat + tau,-2) + (xi + 1) * xi * pow(eps,2) % pow(h,-2) - pow(eps,2) %  pow(h,-1/xi - 2) * (xi + 1);
  }
  res = sum(L) - 1/varsigma;
  return(res);
}


// [[Rcpp::export]]
double f_prime_new(double tau, double tau_hat,
         double varsigma, double xi, double kappa,
         arma::colvec R)
{
  double res;
  arma::colvec h;
  arma::colvec L;
  
  if(fabs(xi) < 1e-8){
    h = exp(-kappa * (R - tau));
    L = kappa * (1 - h);
  }else{
    h = 1 + xi * kappa *(R-tau);
    if(any(h < 0)){
      res = -arma::datum::inf; 
      return(res);
    }
    L = (xi + 1) * kappa * pow(h,-1) - kappa * pow(h,-1/xi - 1);
  }
  res = sum(L) - (tau - tau_hat)/varsigma;
  return(res);
}


// [[Rcpp::export]]
double f_double_prime_new(double tau, double tau_hat,
         double varsigma, double xi, double kappa,
         arma::colvec R)
{
  double res;
  arma::colvec h;
  arma::colvec L;
  
  if(fabs(xi) < 1e-8){
    h = exp(-kappa * (R - tau));
    L = -pow(kappa,2) * h;
  }else{
    h = 1 + xi * kappa *(R-tau);
    if(any(h < 0)){
      res = -arma::datum::inf; 
      return(res);
    }
    L = xi * (xi + 1) * pow(kappa,2) * pow(h,-2) - (xi + 1) * pow(kappa,2) * pow(h,-1/xi - 2);
  }
  res = sum(L) - 1/varsigma;
  return(res);
}

// [[Rcpp::export]]
double g_eta_prime_new(double tau, double tau_hat,
                          double varsigma, double xi, double eta_hat,
                          arma::colvec eps)
{
  double res;
  arma::colvec h;
  arma::colvec L;
  
  if(fabs(xi) < 1e-8){
    h = exp(- exp(eta_hat + tau) * eps);
    L = 1  + log(h) - h % log(h);
  }else{
    h = 1 + xi * exp(eta_hat + tau) * eps;
    if(any(h < 0)){
      res = -arma::datum::inf; 
      return(res);
    }
    L = 1 - (xi + 1)/xi * (h-1) % pow(h,-1) + 1/xi * pow(h,-1/xi) - 1/xi * pow(h,-1/xi - 1);
  }
  res = sum(L) - (tau - tau_hat)/varsigma;
  return(res);
}

// [[Rcpp::export]]
double g_eta_double_prime_new(double tau, double tau_hat,
                       double varsigma, double xi, double eta_hat,
                       arma::colvec eps)
{
  double res;
  arma::colvec h;
  arma::colvec L;
  
  if(fabs(xi) < 1e-8){
    h = exp(- exp(eta_hat + tau) * eps);
    L = log(h) - h % pow(log(h),2) - h % log(h);
  }else{
    h = 1 + xi * exp(eta_hat + tau) * eps;
    if(any(h < 0)){
      res = -arma::datum::inf; 
      return(res);
    }
    L = -(xi + 1)/xi * (h-1) % pow(h,-2) - pow(h,-1/xi)/pow(xi,2) + (xi + 2)/pow(xi,2) * pow(h,-1/xi - 1) - (xi + 1)/pow(xi,2) * pow(h,-1/xi - 2);
  }
  res = sum(L) - 1/varsigma;
  return(res);
}


// [[Rcpp::export]]
arma::colvec gevUpdateM(arma::colvec Y, arma::mat X, arma::colvec M,
			double alpha, double lambda, arma::mat D,
			arma::colvec beta0, arma::mat Omega0)
{
  int i;
  
  //Alex is just messing around here.
  int p = X.n_cols;
  arma::mat M_curr = M;

  //---- Shared covariance ------
  arma::mat C = 1/alpha * exp(-D/lambda);
  for(i = 0; i < C.n_rows; ++i)
    {
      C(i,i) += 1e-5;
    }
  arma::mat C_inv = inv(C);
  //----------------------------

  //-- Stats for current model --
  int p_M_curr = M.n_rows;
  
  //----------------------------
}
