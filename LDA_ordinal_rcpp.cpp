#include <Rcpp.h>
#include <iostream>
#include <ctime>
#include <fstream>
using namespace Rcpp;

/***************************************************************************************************************************/
/*********************************                      UTILS          *****************************************************/
/***************************************************************************************************************************/

// This function helps with multinomial draws
// [[Rcpp::export]]
int whichLessDVPresence(double value, NumericVector prob) {
  int res=prob.length()-1;
  double probcum = 0;
  
  for (int i = 0; i < prob.length(); i++) {
    probcum = probcum + prob(i);
    if (value < probcum) {
      res = i;
      break;
    }
  }
  return res;
}

//' This function helps the sample.theta function
// [[Rcpp::export]]
NumericVector CalcSqDiff(NumericVector z, NumericMatrix media) {
  NumericVector res(media.nrow());
  double tmp;
  
  for(int j=0; j<media.nrow(); j++){ //loop over potential thetas
    tmp=0;
    for(int k=0; k<z.length(); k++){ //loop over species
      tmp=tmp+pow(z[k]-media(j,k),2);
    }
    res[j]=tmp;
  }
  
  return (res);
}

//' This function helps the sample.theta function
// [[Rcpp::export]]
IntegerVector SampleIndTheta(NumericMatrix z, NumericMatrix media, NumericVector lprior,
                             NumericVector runif1) {
  NumericVector res(media.nrow());
  IntegerVector fim(z.nrow());
  double tmp;
  
  for(int i=0; i<z.nrow(); i++){ //loop over locations
    for(int j=0; j<media.nrow(); j++){ //loop over potential thetas
      tmp=0;
      for(int k=0; k<z.ncol(); k++){ //loop over species
        tmp=tmp+pow(z[k]-media(j,k),2);
      }
      res[j]=tmp+lprior[j];
    }
    res=res-max(res);
    res=exp(res);
    res=res/sum(res);
    fim[i]=whichLessDVPresence(runif1[i], res);
  }
  
  return (fim);
}

