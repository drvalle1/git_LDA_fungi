#include <Rcpp.h>
#include <iostream>
#include <ctime>
#include <fstream>
using namespace Rcpp;

/***************************************************************************************************************************/
/*********************************                      UTILS          *****************************************************/
/***************************************************************************************************************************/

//' This function converts vmat into theta
// [[Rcpp::export]]
NumericMatrix convertVtoTheta(NumericMatrix vmat,
                                NumericVector prod) {
  NumericMatrix res(vmat.nrow(),vmat.ncol());
  NumericVector prod1=clone(prod);
  
  for(int j=0; j<vmat.ncol();j++){
    res(_,j)=vmat(_,j)*prod1;    
    prod1=prod1*(1-vmat(_,j));
  }
  
  return (res);
}