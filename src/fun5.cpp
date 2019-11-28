#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
SEXP fun5(NumericMatrix yy_g) {
    Rcpp::NumericMatrix y(yy_g);
    int n=y.nrow();
    Rcpp::NumericVector sumval(n);
    double totval=0;
    for (int j = 0; j < n; j++){
        Rcpp::NumericMatrix::Row yy = y(j,_);
        int y_gi = yy[0];
        double lamda_gi = yy[1];
        double tmp2_gi = yy[2];
        double alpha = yy[3];
        double tmp33;
        double tmp33tot;
        double bsquare =pow(pow(tmp2_gi*(y_gi-alpha)-lamda_gi, 2)-4*tmp2_gi*(-tmp2_gi*alpha*y_gi+lamda_gi), -1/2);
        int max = (tmp2_gi*(y_gi-alpha)-lamda_gi+bsquare)/(2*tmp2_gi);
        tmp33tot= Rf_dnbinom(max, alpha, 1-tmp2_gi,0)*Rf_dpois(y_gi-max,lamda_gi,0);
        for (int t = max+1; t < y_gi+1; t++) {
            tmp33= Rf_dnbinom(t, alpha, 1-tmp2_gi,0)*Rf_dpois(y_gi-t,lamda_gi,0);
            if(tmp33/ tmp33tot>0.01){
                tmp33tot+=tmp33;
            }
            if(tmp33/ tmp33tot<0.01){
                break;
            }
        }
        for (int t = max-1; t > -1; t--) {
            tmp33= Rf_dnbinom(t, alpha, 1-tmp2_gi,0)*Rf_dpois(y_gi-t,lamda_gi,0);
            if(tmp33tot/ tmp33>0.01){
                tmp33tot+=tmp33;
            }
            if(tmp33tot/ tmp33<0.01){
                break;
            }
        }
        sumval[j]=log(tmp33tot);
        totval+=sumval[j];
    }
    return Rcpp::wrap(totval);
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//

