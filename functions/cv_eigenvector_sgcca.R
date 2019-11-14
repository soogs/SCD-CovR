# 25-09-2019 #

# cross validation function for RGCCA::sgcca

# sgcca1 <- RGCCA::sgcca(A,C, c1=c(0.29,1),ncomp=c(3,1), scheme="factorial",verbose=TRUE)

sgcca_cv_general <- function(X, Y, C, c1, ncomp, scheme, verbose, 
                             nrFolds, seed){
  
  set.seed(seed)
  
  #Randomly shuffle the data
  data_shuffle <- cbind(Y, X)[sample(nrow(X)),]
  
  #Create equally sized folds
  folds <- cut(seq(1,nrow(data_shuffle)),breaks=nrFolds,labels=FALSE)
  
  cve_k <- data.frame(error_y = NA)
  
  for(k in 1:nrFolds){
    
    # defining test and train data  
    test_index <- which(folds==k, arr.ind=TRUE)
    data_test <- matrix(data_shuffle[test_index, ], nrow = length(test_index))
    data_train <- data_shuffle[-test_index, ]
    
    A <- list(data_train[,-1], data_train[,1])
    
    # model estimation
    fit <- RGCCA::sgcca(A = A, C = C, c1=c1, ncomp = ncomp, scheme = scheme, verbose= verbose)
    
    W <- fit$astar[[1]]
    
    Tmat <- data_train[,-1] %*% W
    
    reg <- lm(data_train[,1] ~ Tmat)
    
    # eigenvector crossvalidation #
    x_test <- matrix(data_test[,-1], nrow = nrow(data_test))
    y_test <- matrix(data_test[,1], nrow = nrow(data_test))
    
    # y part only #
    pred_y <- x_test %*% W %*% reg$coefficients[-1]
    
    cve_k[k,1] <- mean((pred_y - data_test[,1])^2)
    
  }
  
  cve <- colMeans(cve_k)
  # cve = mean cv error
  
  se <- apply(cve_k,2,sd) / sqrt(nrFolds)
  # se = standard deviation of the mean cv errors per k / sqrt(nrFold)
  
  return(list(cve = cve, se = se, cve_k = cve_k))
}
