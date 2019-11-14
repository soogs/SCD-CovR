#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
arma::mat updateW_cpp (arma::mat Z,
	arma::mat W,
	arma::mat X,
	arma::mat P,
	Rcpp::List blockindex,
	Rcpp::NumericVector cd,
	arma::mat colsumX2, 
	arma::mat colsumP2,
	arma::mat lambda1,
	double lambda2,
	int R) {

	// for colsumX2, colsumP2, lambda1,
	// you have to input them as.matrix

	arma::mat res;
	arma::mat sumPRES;
	Rcpp::NumericVector updorder;
	Rcpp::NumericVector blockobject;
	arma::mat wols;
	arma::mat rkPxXk;
	
  	// initialize replicate index
  	int r;
  	int jj;
  	int j;

	for (r = 0; r < R; r++){
	  res = (Z - X * W * trans(P));
	  sumPRES = res * P.col(r);

	  blockobject = blockindex(cd(r) - 1);

	  updorder = sample(blockobject, blockobject.size(), false);
	  
	  for (jj = 0; jj < updorder.size(); jj++){
	  	
	  	j = updorder(jj) - 1;
	  	// minus 1 here because indexing in cpp starts from 0

	  	//arma::mat woldmat;
	  	arma::mat wold;
	  	arma::mat w;

	  	wold = W(j,r);

	  	//wold = as_scalar(wold_mat);

	  	rkPxXk = (trans(sumPRES) * X.col(j)) + 
	  			 (colsumP2(r) * colsumX2(j) * wold);

	  	wols = sign(rkPxXk) * (abs(rkPxXk) - lambda1(r)*0.5);

	  	w = wols / (lambda2 + (colsumP2(r) * colsumX2(j)));

	  	double absrk = as_scalar(abs(rkPxXk));

	  	if (absrk < lambda1(r)*0.5){
	  		w.zeros();
	  		sumPRES = sumPRES + X.col(j) * colsumP2(r) * wold;
	  	} else {
	  		sumPRES = sumPRES + X.col(j) * colsumP2(r) * (wold - w);

	  	}

	  	W(j,r) = as_scalar(w);
	  			 
	  }
	}
	return W;

}