require(Rcpp)
require(RcppArmadillo)
require(inline)

src <- '
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;
// [[Rcpp::export]]
SEXP testRcpp(NumericVector x_) {
 vec x(x_);
 uvec index = find(x);
 return wrap(x);
}
'
sourceCpp(code=src)

testRcpp(map0)

time0 <- c(0,10,0,20)
ind0 <- time0>0
map0 <- vector("integer",)
map0[ind0] <- as.integer(1:sum(ind0))
map0[!ind0] <- NaN
## which0 <- which(ind0)
which0 <- 1:length(time0)
which0[!ind0] <- NaN


##    uvec find0() { return map0(find(ind0)); }
##    uvec find0(uvec index) { return (map0(index))(find(ind0(index))); }

src <- '
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;
    uvec removeNaN(vec x) {
      vec newx = x;
      uvec index = find(newx == newx); // NaN != NaN
      newx = newx(index);
      return conv_to<uvec>::from(newx);
    }
// [[Rcpp::export]]
uvec testRcpp(NumericVector x_) {
  vec x = as<vec>(x_);
  uvec index = find(x==x);
  return index;
}
'
sourceCpp(code=src)

testRcpp(map0)
testRcpp(which0)



testRcpp(1:10)

src <- '
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;
// [[Rcpp::export]]
double testRcpp(double x) {
 return R::lgammafn(x);
}
'
sourceCpp(code=src)
testRcpp(3.0) - log(gamma(3.0))

