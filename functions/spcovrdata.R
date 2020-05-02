# spcovrdata #
# last modified: 2 May 2020 (only comments) #

# generates data under a PCovR model 

# number of components = 3 #

spcovrdata <- function(VAFx, VAFr = c(), VAFy, I, J, R, sparseamount, Py = c(), modelseed, noiseseed, VAFsum){
  
  #divide the columns of matrix by 2norm
  divideByNorm <- function(mat){
    A <- 1 / matrix(apply(mat, 2, function(x){sqrt(sum(x^2))}),
                    nrow(mat), ncol(mat), byrow = T)
    return(mat * A)
  }
  
  #make a fixed percentage of the non zero coefficients zero (at random) 
  sparsify <- function(comdis, sparsity){
    amounts <- round(apply(comdis, 2, function(x){sum(x != 0)}) * sparsity)
    TF <- apply(comdis, 2, function(x){x != 0})
    where <- apply(TF, 2, which)
    for(i in 1:length(where)){
      comdis[sample(where[[i]], amounts[i]), i]  <- 0
    }
    return(comdis)
  }
  
  
  #orthogonalize two columns of a matrix by using gramm-schmid
  #only use intersection of the non-zero coefficients of the two vectors
  #to enforce that the zero's do not change
  orthogonalize <- function(A, index1, index2){
    int <- intersect(which(A[, index1] != 0), which(A[, index2] != 0))
    
    u <- A[int, index1]
    v <- A[int, index2]
    
    newv <- v - as.numeric((t(u) %*% v) / (t(u) %*% u)) * u
    A[int, index2] <- newv
    
    return(A)
  }
  
  #create an orthonormal matrix 
  #where the zero's are in fixed position
  makeP <- function(A){
    counter <- 0
    #while matrix A is not orthonormal and max iterations is not reached
    #do orthogonalize all columns with respect to each other
    #Do this until the columns are orthogonal, 
    #sometimes multiple passes are needed
    while(TRUE != all.equal(t(A) %*% A, diag(ncol(A))) &&  counter < 1000 ){
      for(i in 2:ncol(A)){
        for(j in 1:(i-1)){
          A <- orthogonalize(A, j, i)
        }
      }
      A <- divideByNorm(A)
      counter  <- counter + 1
      print(counter)
    }
    if(counter < 1000){
      return(list(A=A, status=1))
    } else {
      return(list(A=A, status=0))
    }
  }
  
  
  # scaleData function #
  scaleData <- function(X, value = 0){
    
    X <- scale(X, scale = FALSE)
    attr(X, "scaled:center") <- NULL
    sdX <-  apply(X, 2, function(x) sqrt( sum( x^2 ) / (length(x) - value )   ))  #compute the sd for each column
    
    sdX[sdX == 0] <- 1
    # to account for a column that is completely 0,
    # i make the sd into 1.
    
    sdX <- matrix(sdX, nrow(X), ncol(X), byrow = T)                     #put all the sd's in a matrix
    
    sdX
    
    X <- X * (1 / sdX)      #divide each entry in X by its sd
    return(X)
  }
  
  # VAFcontrol function #  
  VAFcontrol <- function(data, VAFx){
    
    datanoisebigger <- FALSE
    
    iii <- 0
    
    while(datanoisebigger == FALSE){
      iii <- iii + 1
      
      n <- nrow(data)
      p <- ncol(data)
      
      # sum of squares of the X2 dataset
      ssqXtrue <- sum(data^2)
      
      # sample from normal distribution (Ex = Error of X)
      Ex <- matrix(MASS::mvrnorm(n = n, mu = rep(1,p), Sigma = diag(p)),nrow = n, ncol = p)
      
      # centering and scaling the Ex matrix
      Ex <- t(t(Ex) - colMeans(Ex))
      Ex <- scaleData(Ex)
      Ex <- scale(Ex, center = T, scale = F)
      
      # sum of squares of EX
      ssqEx <- sum(Ex^2)
      
      # Rescale noise to desired level
      fx <- sqrt(ssqXtrue*(1-VAFx)/(VAFx * ssqEx))
      
      # 1. the VAFx and SSQ Ex are multiplied
      # 2. (1-VAFx) is multiplied with SSQ Xtrue
      # 3. ratio of these two are calculated
      # 4. square-rooted
      
      data_noise <- data + fx*Ex
      
      datanoisebigger <- (sum(data_noise^2) > sum(data^2))
      
      print(iii)
    }
    
    return(data_noise)
  }
  
  
  # VAFcontrol function #  
  VAFcontroly <- function(data, VAFy, ortho){
    
    datanoisebigger <- FALSE
    
    iii <- 0
    
    while(datanoisebigger == FALSE){
      iii <- iii + 1
      
      n <- nrow(data)
      p <- ncol(data)
      
      # sum of squares of the X2 dataset
      ssqXtrue <- sum(data^2)
      
      # centering the Ex matrix
      ortho <- scaleData(ortho)
      ortho <- scale(ortho, center = T, scale = F)
      
      # sum of squares of EX
      ssqEx <- sum(ortho^2)
      
      # Rescale noise to desired level
      fx <- sqrt(ssqXtrue*(1-VAFy)/(VAFy * ssqEx))
      # fx=sqrt(ssqXtrue*(1-VAFx(vx))/(VAFx(vx)*ssqEx));
      # (matlab code - here vx is an index because VAFx in the matlab code is a vector)
      
      # 1. the VAFx and SSQ Ex are multiplied
      # 2. (1-VAFx) is multiplied with SSQ Xtrue
      # 3. ratio of these two are calculated
      # 4. square-rooted
      
      data_noise <- data + fx*ortho
      
      datanoisebigger <- sum(data_noise^2) > sum(data^2)
      
      print(iii)
    }
    
    return(data_noise)
  }
  
  
  if (sparseamount > 1 | sparseamount < 0){
    print("sparseamount is the proportion of zeros: please provide values from 0 to 1")
    break()
  }
  
  if (is.vector(Py)){
    Py <- matrix(Py, nrow = 1)
  }
  
  # pca model formulation ####
  # ** step 1. specify the sparse and orthogonal matrix V **
  # note V = W = P
  set.seed(modelseed)
  
  V <- matrix(runif(J*R),J)
  
  # 3-component model 
  # 1D, 1D, 1C
  fixW <- matrix(0, nrow = J, ncol = 3)
  fixW[1:(J/2), 1] <- 1
  fixW[(J/2 + 1):J, 2] <- 1
  fixW[,3] <- 1
  
  comdis <- sparsify(fixW, (J*3 * sparseamount - (J))/(J*2))
  
  V[comdis == 0]  <- 0
  
  V <- makeP(V)$A
  
  V <- qr.Q(qr(V))
  
  colSums(abs(V) < 1e-7)
  
  Px <- V
  # you can observe that the off-diagonals are all zeros now
  
  # V'V = I
  t(Px) %*% Px
  
  # ** step 2. specify the U and D **
  # UD = T = XW
  # randomly generate U from multivariate normal
  # (columns of U are not correlated - remember that PCA components are uncorrelated)
  U <- MASS::mvrnorm(n = I, mu = rep(0,R+1), Sigma = diag(R+1), empirical=FALSE)
  
  U <- scale(U, center = T, scale = F)
  
  # orthogonalizing the U
  # by doing this, we achieve U'U = I which is a property of left singular vectors
  # (also means that our components are uncorrelated)
  
  U <- qr.Q(qr(U))
  
  ortho <- U[,R+1] # used for VAFy later
  
  U <- U[,1:R]
  
  # generate the U-matrix for the test set
  Utest <- MASS::mvrnorm(n = I, mu = rep(0,R+1), Sigma = diag(R+1), empirical=FALSE)
  
  Utest <- scale(Utest, center = T, scale = F)
  
  Utest <- qr.Q(qr(Utest))
  
  orthotest <- Utest[,R+1] # used for VAFy later
  
  Utest <- Utest[,1:R]
  
  
  # now we specify the D matPrix
  # within SVD: this is the diagonal matrix with singular values
  # this defines the amount of variance each corresponding principal component has
  # this depends on the VAFr: relative variance accounted for by each component 
  
  # UP TO THIS POINT, THE MODELING HAS TO BE THE SAME ACROSS DIFFERENT REPLICATES ####
  D <- diag(c(VAFsum * VAFr))
  
  
  # so now we have them all:
  # V: V'V = I and sparse
  # U: U'U = I 
  
  # D 
  
  X <- U %*% D %*% t(Px)
  
  y <- U %*% D %*% t(Py)
  
  dattrue <- cbind(y, X)
  
  Xtest <- Utest %*% D %*% t(Px)
  
  ytest <- Utest %*% D %*% t(Py)
  
  set.seed(noiseseed)
  
  Xout <- VAFcontrol(data = X, VAFx = VAFx)
  yout <- VAFcontroly(data = y, VAFy = VAFy, ortho = ortho)
  
  
  datout <- cbind(yout, Xout)
  
  Xtestout <- VAFcontrol(data = Xtest, VAFx = VAFx)
  ytestout <- VAFcontroly(data = ytest, VAFy = VAFy, ortho = orthotest)
  
  dattestout <- cbind(ytestout, Xtestout)
  
  
  returnobj <- list(datout = datout, dattestout = dattestout, dattrue = dattrue, U = U, Px = Px, D = D, VAFx = VAFx, VAFy = VAFy, VAFr = VAFr, Py = Py)
  return(returnobj)
}

