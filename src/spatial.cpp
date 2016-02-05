#include <Rcpp.h>

using namespace Rcpp;

// [[Rcpp::export]]
IntegerMatrix A(int p)
{
  IntegerMatrix A(p,p);
  return(A);
}
