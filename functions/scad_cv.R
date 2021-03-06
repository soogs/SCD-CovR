
scad_cv <- function(X, R, ridge, lasso, nrFolds, fixW,
                    MAXITER = 10000, stop_value = 1e-10, seed){
  
  set.seed(seed)
  
  shuffle_index <- sample(nrow(X))
  
  #Randomly shuffle the data
  data_shuffle <- X[shuffle_index,]
  
  #Create equally sized folds
  folds <- cut(seq(1,nrow(data_shuffle)),breaks=nrFolds,labels=FALSE)
  
  
  cve_k <- data.frame(error_x = NA)
  
  for(k in 1:nrFolds){
    
    # defining test and train data  
    test_index <- which(folds==k, arr.ind=TRUE)
    data_test <- matrix(data_shuffle[test_index, ], nrow = length(test_index))
    data_train <- data_shuffle[-test_index, ]
    
    # model estimation
    fit <- sparseSCAcpp(X = X, Q = R, 
                        RIDGE = ridge, LASSO = lasso, 
                        fixW = fixW, 
                        maxItrOuterloop = 100000, nStarts = 1,
                        print = FALSE, tol = 10^-8)
    
    # eigenvector crossvalidation #
    x_test <- matrix(data_test, nrow = nrow(data_test))
    
    pred_x <- matrix(NA, nrow(x_test), ncol(x_test))
    
    # x part first #
    for(j in 1:ncol(x_test)){
      TMinusVariableJ <- x_test[,-j] %*% fit$W[-j,]
      pred_x[,j] <- TMinusVariableJ %*% fit$P[j, ] 
    }
    
    cve_k[k,1] <- mean((pred_x - data_test)^2)
    
  }
  
  cve <- colMeans(cve_k)
  # cve = mean cv error
  
  se <- apply(cve_k,2,sd) / sqrt(nrFolds)
  # se = standard deviation of the mean cv errors per k / sqrt(nrFold)
  
  return(list(cve = cve, se = se, cve_k = cve_k))
}
