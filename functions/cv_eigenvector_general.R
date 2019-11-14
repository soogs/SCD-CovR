# 25-09-2019 #

# adjusting Niek's eigenvector CV function
# so that it suits what i am doing 
# a general CV function that can be used for R, alpha, ridge, lasso

sscovr_cv_general <- function(X, Y, R, ridge, lasso, lambda_y, nrFolds,
                  cdstructure, alpha, gamma = 3.7,
                  inits = c("rational", "oracle", "multistart"), nrstart = 1,
                  MAXITER = 10000, stop_value = 1e-10, blockcols, seed){
  
  set.seed(seed)
  
  #Randomly shuffle the data
  data_shuffle <- cbind(Y, X)[sample(nrow(X)),]
  
  #Create equally sized folds
  folds <- cut(seq(1,nrow(data_shuffle)),breaks=nrFolds,labels=FALSE)
  
  
  cve_k <- data.frame(error_x = NA, 
                      error_y = NA,
                      error_all = NA)
  
  for(k in 1:nrFolds){
    
    # defining test and train data  
    test_index <- which(folds==k, arr.ind=TRUE)
    data_test <- matrix(data_shuffle[test_index, ], nrow = length(test_index))
    data_train <- data_shuffle[-test_index, ]
    
    # model estimation
    fit <- sscovr2(X = data_train[,-1], Y = data_train[,1],
                   blockcols = blockcols, 
                   R = R, cdstructure = cdstructure, 
                   alpha = alpha, 
                   lambda1 = lasso, 
                   lambda2 = ridge, 
                   inits = "rational", nrstart = 1, MAXITER = 10000, 
                   stop_value = stop_value, lambda_y = lambda_y, ridge_in_between = T, ridge_end = T)
    
    # eigenvector crossvalidation #
    x_test <- matrix(data_test[,-1], nrow = nrow(data_test))
    y_test <- matrix(data_test[,1], nrow = nrow(data_test))
    
    pred_x <- matrix(NA, nrow(x_test), ncol(x_test))
    
    # x part first #
    for(j in 1:ncol(x_test)){
      TMinusVariableJ <- x_test[,-j] %*% fit$cd_results[[1]]$W[-j,]
      pred_x[,j] <- TMinusVariableJ %*% fit$cd_results[[1]]$P[j, ] 
    }
    
    # y part #
    # since y = XWPy, i think this part would be the same as normal cross-validation #
    pred_y <- x_test %*% fit$cd_results[[1]]$W %*% fit$cd_results[[1]]$P[1,] 
    
    cve_k[k,1] <- mean((pred_x - data_test[,-1])^2)
    cve_k[k,2] <- mean((pred_y - data_test[,1])^2)
    cve_k[k,3] <- mean((cbind(pred_y, pred_x) - data_test)^2)
    
  }
  
  cve <- colMeans(cve_k)
  # cve = mean cv error
  
  se <- apply(cve_k,2,sd) / sqrt(nrFolds)
  # se = standard deviation of the mean cv errors per k / sqrt(nrFold)

  return(list(cve = cve, se = se, cve_k = cve_k))
}
