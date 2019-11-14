

# setwd("/home/soogs/Desktop/PhD 1 - SC regression/Project 2/SCaDs Niek/SCaDS-master/How_to_use_SCaDS/")
# 
# # Load in needed packages
# library(Rcpp)
# library(MASS)
# 
# source("./CVfunction.R")
# sourceCpp("./sparseSCA.cpp")
# 
# 
# dat <- spcovrdata(VAFx = 0.7, VAFr = c(0.9, 0.05, 0.05), VAFy = 0.2, I = 1000, 
#                   J = 200, R = 3, sparseamount = 0.4, Py = c(1, 0.6, -0.02), 
#                   modelseed = 111, seed = 22, VAFsum = 30)
# 
# 
# fixW <- matrix(1, nrow = ncol(dat$datout[,-1]), ncol = 3)
# 
# 
# hi <- sparseSCAcpp(X = dat$datout[,-1], Q = 3, RIDGE = 0.1, LASSO = rep(0.1, 3), 
#              fixW = fixW, maxItrOuterloop = 100000, nStarts = 1,
#              print = TRUE, tol = 10^-8)
# 
# hey <- findLasso_SCaDs(dat = dat$datout[,-1], zeros = colSums(dat$Px ==0), R = 3, init = 0.1, ridge = 1e-6, maxiterOut = 20, maxiterIn = 100)
# 
# hi2 <- sparseSCAcpp(X = dat$datout[,-1], Q = 3, RIDGE = 1e-6, LASSO = hey$lasso, 
#                    fixW = fixW, maxItrOuterloop = 100000, nStarts = 1,
#                    print = TRUE, tol = 10^-8)
# 
# colSums(hi2$W == 0)

# Rcpp::sourceCpp("./sparseSCA.cpp")

findLasso_SCaDs <- function(dat, zeros, R, fixW, init, ridge = 1e-6, maxiterOut, maxiterIn){
  # zeros should be a vector
  
  estimatedzeros <- rep(0,R)
  lasso <- rep(init,R)
  
  converged <- FALSE
  
  iterOut <- 0
  
  repeated <- 0
  
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
      
      lassochanging <- TRUE
      
      while(abs(zeros[j] - estimatedzero) > 0 && iterIn <= maxiterIn && lassochanging){
        
        # if the lasso penalty fluctuates by a very small amount,
        # quit the thing
        if ((up - down) < 1e-04){
          lassochanging <- FALSE
        }
        
        iterIn <- iterIn + 1
        lasso[j] <- (down + up)/2
        
        # lasso[-j] <- 0
        
        fit <- sparseSCAcpp(X = dat, Q = R, RIDGE = ridge, LASSO = lasso, 
                                      fixW = fixW, maxItrOuterloop = 100000, nStarts = 1,
                                      print = FALSE, tol = 10^-8)
        estimatedzero <- sum(abs(fit$W[,j]) < 1e-06)
        
        if(zeros[j] > estimatedzero){
          down  <- lasso[j]
          # if the estimated zeros are not enough,
          # pull up the 'down'
        } else if (zeros[j] < estimatedzero){
          up  <- lasso[j]
          # if the estimated zeros are more than enough,
          # pull down the 'up'
        }
        # else (don't do anything)
        
        print(round(lasso,10))
        
        if (!(abs(zeros[j] - estimatedzero) > 0)){
          print("number of zeros correct - next penalty")
        }
        
        if (!(iterIn <= maxiterIn)){
          print("max number of iterations reached - next penalty")
        }
        
        if (!lassochanging){
          print("lasso does not change much - next penalty")
        }
        
      }
      
      if(out_round_count == R){
      estimatedzeros <- apply((abs(fit$W) < 1e-6), 2, sum)
      }
    }
     
    
  }

  
  if( iterOut < maxiterOut && iterIn < maxiterIn ){
    converged <- TRUE
  }
  return(list(lasso = lasso, converged = converged))
}



# if(iterOut == maxiterOut){
#   if (sum((estimatedzeros == ncol(dat)) == 0) > 0){
#     iterOut <- iterOut - 1
#     repeated <- repeated + 1
#     
#     if(repeated == 3){
#       stop("too many repeats due to zero-column. adjust the fold shuffle")
#     }
#   }
# }