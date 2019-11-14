# pretty much copying the spcovr matlab function from the bmc paper #

# this function works (checked against manual calculation) #
# but a number of things have been changed for the comparison, 
# so please alter those in order to use this function properly

# the loss is calculated without dividing by 'sum of squares of Z'

# edit 18042019: i do not allow a column full of zeros in the W matrix,
# because this will disable the computation of the P matrix using the kernel trick

# edit 23042019: i put different lambda values for each column of the component

# Rcpp::sourceCpp( "C:\\Users\\TSB-MTO\\Desktop\\Soogeun\\updateW.cpp")

sscovr2 <- function(X, Y, blockcols, R, cdstructure = NULL, alpha, 
                   lambda1, lambda2, gamma = 3.7, inits = c("rational", "oracle", "multistart"), 
                   nrstart = 10, MAXITER = 10000, stop_value = 1e-10, lambda_y, ridge_in_between, ridge_end){
  
  # blockcols = vector specifying the number of variables that each data block has
  
  # 1. define sub-functions ####
  
  # updateW function #
  # coordinate descent with different penalties #
  # adapted katrijn's matlab code to suit my purposes: scad, mcp #
  
  # if (!is.matrix(cdstructure)){
  #   print ("arg 'cdstructure' has to be a matrix: please provide a matrix")
  #   stop()
  # }
  
  if (is.vector(cdstructure)){
    cdstructure <- matrix(cdstructure, nrow = 1)
  }
  
  if (length(lambda1) != R){
    stop("Vector of length R is required as an input for the lasso penalty")
  }
  
  # updatePx function #
  # directly from katrijn's matlab code #
  
  # !!! check if answer identical to original procrustes rotation trick !!! # 
  updatePx <- function(wX, W, X){
    
    if (sum(colSums(W != 0) == 0) > 0){
      print ("ERROR: W matrix has zero-columns. This disallows the P computation. Try a lower l1 penalty.")
      stop()
    }
    
    Tmat <- X %*% W
    K1 <- t(wX) %*% Tmat
    K <- t(Tmat) %*% wX %*% t(wX) %*% Tmat
    eigs <- eigen(K)
    V <- eigs$vectors
    Ssq <- eigs$values
    
    S <- diag(Ssq^(-0.5))
    
    if(ncol(W) == 1){
      S <- matrix(Ssq^(-0.5), ncol = 1, nrow = 1)
    }
    
    Px <- K1 %*% V %*% S %*% t(V)
    
    # svdP <- svd(t(wX) %*% wX %*% W)
    # 
    # Px <- svdP$u %*% t(svdP$v)
    # 
    # svdP <- svd(t(W) %*% t(X) %*% wX)
    # 
    # Px <- svdP$v %*% t(svdP$u)
    
    return (Px)
  }
  
  # loss calculation function #
  losscal <- function(Z, X, W, P, lambda1, lambda2, ssZ){
    lambda1mat <- matrix(0, nrow = dim(W)[1], ncol = dim(W)[2])
    
    for (r in 1:length(lambda1)){
      lambda1mat[,r] <- lambda1[r]
    }
    
    result <- (sum((Z - X %*% W %*% t(P))^2) + (sum(lambda1mat * abs(W))) + (lambda2 * sum(W^2)) + (lambda_y * sum(P[1,]^2))) /ssZ
    return (result)
  }
  
  # pcovr function from bmc bioinformatics #
  
  pcovr <- function(X, Y, R, alpha){
    
    # if Y is provided as a vector.. #
    if (is.vector(Y)){
      Y <- matrix(data = Y, ncol = 1)
    }
    
    I <- nrow (X)
    Jx <- ncol (X)
    Jy <- ncol (Y)
    J <- Jx + Jy
    eps <- 1e-12
    
    # [Iy,Jy]=size(Y);
    # if Iy~=I, disp(' size Y and X do not match ');return;end;
    
    # weighting X and Y according to the alpha parameter
    w1 <- (I*J*alpha) / (sum(Y^2))
    w2 <- (I*J*(1-alpha)) / (sum(X^2))
    
    wY <- sqrt(w1) * Y
    wX <- sqrt(w2) * X
    
    # Z = [wY wX]
    Xstar <- cbind(wY, wX)
    
    # SVD data #
    if (J > I){
      XX <- X %*% t(X)
      eigs <- eigen(XX)
      Ux <- eigs$vectors
      Ssq <- eigs$values
      Sxsq <-Ssq
      Sxsq[Sxsq < 1e-6] <- 1e-6
      Sx <- sqrt(Sxsq)
      
      invSx <- Sx^(-1)
      
      Vx <- t(X) %*% Ux %*% diag(invSx)
      
      Uxstar <- t(Ux) %*% Xstar
      
      svd1 <- svd(Uxstar)
      U <- svd1$u[,1:R]
      S <- diag(svd1$d[1:R])
      V <- svd1$v[,1:R]
      
      
      if (R == 1){
        U <- matrix(svd1$u[,1:R], ncol = 1)
        S <- matrix(svd1$d[1:R], nrow = 1, ncol = 1)
        V <- matrix(svd1$v[,1:R], ncol = 1)
      }
      
    } else if (I >= J){
      XX <- t(X) %*% X
      eigs <- eigen(XX)
      
      Vx <- eigs$vectors
      Ssq <- eigs$values
      Sxsq <- Ssq
      Sxsq[Sxsq < 1e-6] <- 1e-6
      
      Sx <- sqrt(Sxsq)
      
      invSx <- Sx^(-1)
      
      Ux <- X %*% Vx %*% diag(invSx)
      Uxstar <- t(Ux) %*% Xstar
      
      svd1 <- svd(Uxstar)
      U <- svd1$u[,1:R]
      S <- diag(svd1$d[1:R])
      V <- svd1$v[,1:R]
      
      if (R == 1){
        U <- matrix(svd1$u[,1:R], ncol = 1)
        S <- matrix(svd1$d[1:R], nrow = 1, ncol = 1)
        V <- matrix(svd1$v[,1:R], ncol = 1)
      }
    }
    
    P <- V
    
    W <- Vx %*% diag(invSx) %*% U %*% S
    
    # reweighting step #
    Py <- P[1:Jy,] / (sqrt(w1))
    Px <- P[-c(1:Jy),] / (sqrt(w2))
    
    # if Py turns out to be a vector,
    # we need to make it into a matrix of one row 
    if (is.vector(Py)){
      Py <- matrix(data = Py, nrow = 1)
    }
    
    # fit measure #
    RsqX <- 1-sum(sum((X - X %*% W %*% t(Px))^2))/(sum(sum(X^2)))
    Rsqy <- 1-sum(sum((Y - X %*% W %*% t(Py))^2))/(sum(sum(Y^2)))
    
    return_list <- list(W = W, Px = Px, Py = Py, RsqX = RsqX, Rsqy = Rsqy)
    
    return (return_list)
  }
  
  # 2. define a few objects ####
  
  # Y could be provided as a column vector. transfer into matrix
  if (is.vector(Y)){
    Y <- matrix(data = Y, ncol = 1)
  }
  
  I <- nrow(X)
  Jx <- ncol(X)
  Jy <- ncol(Y)
  J <- Jx + Jy
  
  # weighting X and Y according to the alpha parameter
  w1 <- (I*J*alpha) / (sum(Y^2))
  w2 <- (I*J*(1-alpha)) / (sum(X^2))
  
  wY <- sqrt(w1) * Y
  wX <- sqrt(w2) * X
  
  # Z = [wY wX]
  Z <- cbind(wY, wX)
  
  # blockindex definition #
  blockcols2 <- cumsum(blockcols)
  
  blockindex <- list()
  
  blockindex[[1]] <- 1:blockcols2[1]
  
  for (i in 2:length(blockcols)){
    blockindex[[i]] <- (blockcols2[i-1]+1):blockcols2[i]
  }
  
  blockindex[[length(blockcols)+1]] <- 1:ncol(X)
  # the last object in the list is all of the column indices of the data 
  # (this represents a common component)
  
  # 3. cdstructure specification ####
  # creating the pattern of common and distinctive weights
  
  if (is.null(cdstructure)){
    # c, d1, d2 (d1 = distinctive for block 1..)
    # if there are 4 blocks and 2 components:
    # c, d1, d2, d3, d4.
    # possibilities are: cc, cd1, cd2, cd3, cd4, d1c, ...
    
    # 0 = common -> later becomes (ncol(X)+1) for indexing
    # 1 = distinctive 1
    # 2 = distinctive 2
    # ...
    
    # browser()
    
    # all of the combinations of common and distinctive 
    cd_all <- gtools::permutations(n = length(blockcols) + 1, 
                                   r = R, 
                                   v = 0:length(blockcols), 
                                   repeats.allowed = TRUE)
    
    # zeros, representing the common process, are replaced by (ncol(X)+1)
    # in this way, it can be indexed (because object[[0]] is impossible)
    cd_all[cd_all == 0] <- length(blockcols) + 1
    
  } else { # if a specific cdstructure is provided:
    cd_all <- matrix(data = cdstructure, nrow = nrow(cdstructure))
    
    cd_all[cd_all == 0] <- length(blockcols) + 1
  }
  
  cd_results <- list()
  cd_loss <- c()
  
  # 4. initial values generated ####
  for (cdindex in 1:nrow(cd_all)){
    
    # current common-distinctive pattern #
    cd <- cd_all[cdindex,]
    
    # initial value #
    # for both W and P, we have to provide initial values # 
    # borrowing the initial value system from the spca_adj function # 
    
    if (inits == "rational"){
      # W and P from pcovr
      pcovr_results <- pcovr(X = X, Y = Y, R = R, alpha = alpha)
      
      W <- pcovr_results$W
      # Px <- pcovr_results$Px
      # Py <- pcovr_results$Py
      # P <- rbind(Py, Px)
      
      nrstart <- 1
    } 
    
    if (inits == "oracle"){
      # user-specified matrices for W and P
      nrstart <- 1
    } 
    
    # vector to save results later
    multi_results <- list()
    multi_loss <- c()
    
    
    for (nr in 1:nrstart){
      
      if (inits == "multistart"){
        # initial values W, Px (orthogonal), Py
        
        W <- matrix(stats::runif(n = Jx*R, min = -1, max = 1), nrow = Jx, ncol = R)
        # Px <- matrix(stats::runif(n = Jx*R, min = -5, max = 5), nrow = Jx, ncol = R)
        # Px <- qr.Q(qr(Px))
        # 
        # Py <- matrix(stats::runif(n = Jy*R, min = -5, max = 5), nrow = Jy, ncol = R)
        # P <- rbind(Py, Px)
      }
      
      # initial values for P from SVD #
      svdwX <- svd(wX)
      Px <- svdwX$v[,1:R]
      
      if (R == 1){
        Px <- matrix(svdwX$v[,1:R], ncol = 1)
      }
      
      Py <- matrix(stats::runif(n = Jy*R, min = -1, max = 1), nrow = Jy, ncol = R)
      # Py <- matrix(glmnet::glmnet(x = X %*% W, y = wY, alpha = 0,  lambda = 0.5)$beta, ncol = 1)
      P <- rbind(Py, Px)
      
      # distinctive components become sparse
      for (r in 1:R){
        W[-blockindex[[cd[r]]],r] <- 0
      }
      
      colsumX2 <- colSums(X^2)
      
      ssZ <- sum(Z^2)
      
      # initial loss
      # loss0 <- losscal(Z, X, W, P, lambda1, lambda2, ssZ)
      
      loss0 <- 10000
      
      # convergence starting 
      conv <- 0
      iter <- 1
      
      loss_hist <- c(loss0)
      
      # 5. estimation ####
      while (conv == 0){
        
        colsumP2 <- colSums(P^2)
        
        # W given P #
        W <- updateW_cpp(Z = Z, W = W, X = X, P = P, R = R, lambda1 = as.matrix(lambda1), lambda2 = lambda2, colsumX2 = as.matrix(colsumX2), colsumP2 = as.matrix(colsumP2), blockindex = blockindex, cd = cd)
        
        hi <- colSums((X %*% W)^2)
        
        loss1 <- losscal(Z, X, W, P, lambda1, lambda2, ssZ)
        
        # current loss is always smaller than the previous loss
        if ((loss0 - loss1) < stop_value){
          conv <- 1
        }
        
        if (loss0 < loss1){
          print ("ERROR: current loss > previous loss, after weights estimation")
          stop()
        }
        
        iter <- iter + 1
        
        loss_hist[iter] <- loss1
        
        loss0 <- loss1
        
        # P given W #
        Px <- updatePx(wX = wX, X = X, W = W)
      
        
        # Py <- matrix(glmnet::glmnet(x = X %*% W, y = wY, alpha = 0,  lambda = 0.1)$beta, ncol = 1)
        
        # Py <- MASS::lm.ridge(wY ~ X %*% W, lambda = 0.1)
        
        if (ridge_in_between){
          jere <- solve(t(X %*% W) %*% (X %*% W) + diag(lambda_y, R)) %*% t(X %*% W) %*% wY
          
          Py <- jere
          
          } else {
          Py <- MASS::ginv(X %*% W) %*% wY
          
          }
        
                
        P <- rbind(t(Py), Px)
        
        
        
        # browser()
        
        loss1 <- losscal(Z, X, W, P, lambda1, lambda2, ssZ)
        
        # current loss is always smaller than the previous loss
        if ((loss0 - loss1) < stop_value){
          conv <- 1
        }
        
        if (loss0 < loss1){
          print ("ERROR: current loss > previous loss, after loadings estimation")
          stop()
        }
        
        iter <- iter + 1
        
        loss_hist[iter] <- loss1
        
        loss0 <- loss1
      }
      
      # one last step
      # W given P #
      W <- updateW_cpp(Z = Z, W = W, X = X, P = P, R = R, lambda1 = as.matrix(lambda1), lambda2 = lambda2, colsumX2 = as.matrix(colsumX2), colsumP2 = as.matrix(colsumP2), blockindex = blockindex, cd = cd)
      
      
      loss1 <- losscal(Z, X, W, P, lambda1, lambda2, ssZ)
      
      # current loss is always smaller than the previous loss
      if ((loss0 - loss1) < stop_value){
        conv <- 1
      }
      
      if (loss0 < loss1){
        print ("ERROR: current loss > previous loss, after weights estimation")
        stop()
      }
      
      iter <- iter + 1
      
      loss_hist[iter] <- loss1
      
      loss0 <- loss1
      
      # P given W #
      Px <- updatePx(wX = wX, X = X, W = W)
      Py <- MASS::ginv(X %*% W) %*% wY
      
      # Py <- matrix(glmnet::glmnet(x = X %*% W, y = wY, alpha = 0,  lambda = 0.1)$beta, ncol = 1)
      
      # Py <- MASS::lm.ridge(wY ~ X %*% W, lambda = 0.1)
      
      if (ridge_end){
        Py <- solve(t(X %*% W) %*% (X %*% W) + diag(lambda_y, R)) %*% t(X %*% W) %*% wY
        
        Py_reg <- MASS::ginv(X %*% W) %*% wY / (sqrt(w1))
      } else {
        Py_reg <- Py / (sqrt(w1))
      }
      
      P <- rbind(t(Py), Px)
      
      
      
      # reweighting step #
      Py <- matrix(P[1:Jy,] / (sqrt(w1)), ncol = R)
      Px <- P[-c(1:Jy),] / (sqrt(w2))
      
      if (R == 1){
        Px <- matrix(P[-c(1:Jy),] / (sqrt(w2)), ncol = 1)
      }
      
      P <- rbind(Py, Px)
      
      result_list <- list(W = W, P = P, loss = loss1, loss_hist = loss_hist, iter = iter, Py_reg = Py_reg)
      multi_loss[nr] <- loss1
      multi_results[[nr]] <- result_list
      
    }
    
    lossmin <- which.min(multi_loss)
    multi_results_min <- multi_results[[lossmin]]
    
    # browser()
    
    cd_results[[cdindex]] <- multi_results_min
    cd_loss[cdindex] <- multi_loss[lossmin]
  }
  
  return_results <- list(cd_results = cd_results, cd_loss = cd_loss, cdstructures = cd_all)
  return (return_results)
}
