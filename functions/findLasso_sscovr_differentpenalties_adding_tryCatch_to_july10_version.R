findLasso_sscovr <- function(X, Y, zeros, R, init, ridge = 1e-6, blockcols, cdstructure, alpha, maxiterOut, maxiterIn){
  # zeros should be a vector
  
  if (length(zeros) != R){
    stop("'zeros' should be a vector of length R indicating the numbers of zeros corresponding to each component")
  }
  
  estimatedzeros <- rep(0,R)
  lasso <- rep(0,R)
  
  converged <- FALSE
  
  iterOut <- 0
  
  while(sum(abs(zeros - estimatedzeros)) > 0 && iterOut <= maxiterOut ){
    iterOut <- iterOut + 1
    
    # this mixes the order of the component
    mixedorder <- sample(R)
    
    
    
    out_round_count <- 0
    for (j in mixedorder){
      out_round_count <- out_round_count + 1
      
      iterIn <- 0
      
      up <- init
      down <- 0
      
      estimatedzero <- 0
      
      # lasso_old <- init
      
      lassochanging <- TRUE
      
      # onarow <- 0
      
      lasso_old <- 0
      
      while(abs(zeros[j] - estimatedzero) > 1 && iterIn <= maxiterIn && lassochanging){
        
        iterIn <- iterIn + 1
      
        # intitial lasso = 0
        lasso[j] <- (down + up)/2
        
        lasso_new <- lasso[j]
        
        # if sscovr yields error, return fit <- NA #
        fit_try <- tryCatch(sscovr2(X = X, Y = Y, blockcols = blockcols, R = R, cdstructure = cdstructure, 
                               alpha = alpha, lambda1 = lasso, lambda2 = ridge, inits = 'rational', 
                               nrstart = 1, MAXITER = 10000, stop_value = 1e-7, lambda_y = 1, ridge_in_between = T, ridge_end = T)
                            , error = function(e) NA)
        
        # if sscovr does not run because W has a zero-column,
        # pull down the 'up'
        if (sum(is.na(fit_try)) > 0){
          
          up  <- up - up/2
          iterIn <- iterIn - 1 
          print("error with sscovr function, lasso decreasing")
          
        } else {
          
          fit <- fit_try
          
          estimatedzero <- sum(abs(fit$cd_results[[1]]$W[,j]) < 1e-06)
          
          if(zeros[j] > estimatedzero){
            # lasso[j] <- lasso[j] + (up - lasso[j])/(10 + goingdown/2)
            
            down  <- lasso[j]
            
            # if the estimated zeros are not enough,
            # pull up the 'down'
          } else if (zeros[j] < estimatedzero){
            # lasso[j] <- lasso[j] - (lasso[j] - down)/(10 + goingup/2)
            
            up  <- lasso[j]
            # if the estimated zeros are more than enough,
            # pull down the 'up'
          }
          # else (don't do anything)
          
        }
        
        print("current zeros:")
        
        if (sum(is.na(fit_try)) > 0){
          print ("sscovr error, will calculate with a smaller lasso")
        } else {
          print(colSums(fit$cd_results[[1]]$W == 0))
        }
        
        print ("new lasso:")
        
        print(round(lasso,10))
        
        
        # if the lasso penalty fluctuates by a very small amount,
        # quit the thing
        if ((abs(lasso_new - lasso_old)) < 0.001){
          lassochanging <- FALSE
        }
        
        if (!lassochanging){
          print("lasso does not change much - next penalty")
        }
        
        if (!(abs(zeros[j] - estimatedzero) > 1)){
          print("difference from the aim is 1 - next penalty")
        }
        
        if (!(iterIn <= maxiterIn)){
          print("max number of iterations reached - next penalty")
        }
        
        
        
        lasso_old <- lasso_new
        
      }
      
      if(sum(is.na(fit_try)) == 0){
        if(out_round_count == R){
          estimatedzeros <- apply((abs(fit$cd_results[[1]]$W) < 1e-6), 2, sum)
        } 
      }
      
    }
    
  }
  
  if (sum(is.na(fit_try)) == 0 ){
    if(sum(apply((abs(fit$cd_results[[1]]$W) < 1e-6), 2, sum) == zeros) == R){
      converged <- TRUE
    }
  }
  
  return(list(lasso = lasso, converged = converged))
}