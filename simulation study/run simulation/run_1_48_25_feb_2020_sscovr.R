# revision of simulation for project 2 #
# 25-feb-2020 #
# last modified: 2 May 2020 #

# batch 1 (low dimensional): 1 ~ 48

# i test every condition 3 times: 

# 1. functions loaded ####
Rcpp::sourceCpp("./sparseSCA.cpp")
Rcpp::sourceCpp("./updateW.cpp")
source("./conditions_making.R")
# source("./evaluation_criteria.R")
source("./sscovr_looking_in.R")
source("./spcovrdata.R")
source("./cv_eigenvector_general.R")
source("./findLasso_sscovr_differentpenalties_adding_tryCatch_to_july10_version.R")
source("./findLasso_SCaDs.R")
source("./scad_cv.R")


# 2. replication starts ####
set.seed(98) # previously i used seed 6397, 2212
model_seed <- sample(x = 1:100000, size = nrow(condition_df))
noise_seed <- sample(x = 1:100000, size = nrow(condition_df))

reps_to_do <- 1:48

needed_objects <- ls()

needed_objects <- ls()

needed_objects <- append(needed_objects, "rrr")

for (rrr in reps_to_do){ 
  
  time1 <- Sys.time()
  
  modelseeding <- model_seed[rrr]
  noiseseeding <- noise_seed[rrr]
  
  cond <- condition_df[rrr,]
  
  # if (sum(cond[-length(cond)] == prevcond[-length(cond)]) != 6){
  #   model_selection <- TRUE
  # }
  
  if (cond$dimension == "low"){
    I <- 100
    J <- 20
  }
  
  if (cond$dimension == "high"){
    I <- 100
    J <- 200
  }
  
  if (cond$weak == "common"){
    VAFr <- c(0.5, 0.4, 0.1)
  }
  
  if (cond$weak == "distinctive"){
    VAFr <- c(0.1, 0.4, 0.5)
  }
  
  if (cond$components == 3){
    cd <- matrix(c(1,2,0), nrow = 1)
  }
  
  if (cond$relevant == "common"){
    
    if (cond$weak == "distinctive"){
      Py <- c(1.2, -0.3, 0.48)
    } 
    
    if (cond$weak == "common"){
      Py <- c(0.24, -0.3, 2.4)
    }
    
  }
  
  
  if (cond$relevant == "distinctive"){
    
    if (cond$weak == "common"){
      Py <- c(0.48, -0.3, 1.2)
    }
    
    if (cond$weak == "distinctive"){
      Py <- c(2.4, -0.3, 0.24)
    }
  }
  
  
  cond$components <- as.numeric(cond$components)
  cond$vafx <- as.numeric(cond$vafx)
  cond$vafy <- as.numeric(cond$vafy)
  cond$reps <- as.numeric(cond$reps)
  
  # data generating #
  dat <- spcovrdata(VAFx = cond$vafx, VAFr = VAFr, 
                    VAFy = cond$vafy, I = (I*2), J = J, R = 3,
                    sparseamount = 0.7, modelseed = modelseeding, 
                    noiseseed = noiseseeding, VAFsum = 100, Py = Py)
  
  # separating out train and test set
  train_index <- sample(1:nrow(dat$datout), size = I)
  
  X_train <- dat$datout[train_index, -1]
  X_test <- dat$datout[-train_index, -1]
  
  y_train <- dat$datout[train_index, 1]
  y_test <- dat$datout[-train_index, 1]
  
  
  # sscovr #
  
  # if estimating 3 components,
  # only CV for alpha,
  # if estimating 2 components,
  # CV for CD as well 
  
  alpha_range <- c(0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
  
  ridge_range <-  c(100, 50, 30, 10, 5, 1, 0.5, 0.1)
  
  ranges <- expand.grid(ridge_range, alpha_range)
  
  colnames(ranges) <- c("ridge", "alpha")
  
  alpha_seed <- sample(nrow(ranges))
  
  alpha_cve <- matrix(NA, ncol = 3+1+1, nrow = nrow(ranges))
  alpha_se <- matrix(NA, ncol = 3+1+1, nrow = nrow(ranges))
  
  for (i in 1:nrow(ranges)){
    cv_i <- sscovr_cv_general(X = X_train, Y = y_train, 
                              ridge = ranges[i,]$ridge, lasso = rep(0, cond$components), nrFolds = 10,
                              R = cond$components, 
                              cdstructure = matrix(rep(0, cond$components), nrow = 1), 
                              alpha = ranges[i,]$alpha, 
                              inits = "rational", nrstart = 1, 
                              blockcols = c(J/2, J/2), 
                              seed = alpha_seed[i], lambda_y = 1, stop_value = 1e-7)
    
    alpha_cve[i,] <- c(cv_i$cve, ranges[i,]$alpha, ranges[i,]$ridge)
    alpha_se[i,] <- c(cv_i$se, ranges[i,]$alpha, ranges[i,]$ridge)
    
    print(rep(i, 10))
  }
  
  # minimum alpha + 1se value:
  alpha_1se <- min(alpha_cve[,2]) + alpha_se[which.min(alpha_cve[,2]),2]
  
  # smallest alpha value, from the cv errors smaller than 1se
  alpha_ridge_chosen <- matrix(alpha_cve[(alpha_cve[,2] < alpha_1se), ], ncol = 5)[1,4:5]
  
  alpha_ridge_chosen <- data.frame(alpha = alpha_ridge_chosen[1],
                                   ridge = alpha_ridge_chosen[2])
  
  # alpha_chosen <- alpha_cve[min(which(alpha_cve[,2] < alpha_1se)),4]
  
  sscovr_alpha_ridge <- alpha_ridge_chosen
  
  # cross validation for common-distinctive:
  # only for estimating 2-components
  # consider only: 1D2D, 1D C, 2D C, C C
  
  if(cond$components == 2){
    cd_range <- matrix(c(1,2,
                         1,0,
                         2,0,
                         0,0), ncol = 2, byrow = T)
    
    cd_seed <- sample(1:1000, nrow(cd_range))
    
    cd_cve <- matrix(NA, ncol = 3+2, nrow = nrow(cd_range))
    cd_se <- matrix(NA, ncol = 3+2, nrow = nrow(cd_range))
    
    for (i in 1:nrow(cd_range)){
      cv_i <- sscovr_cv_general(X = X_train, Y = y_train, 
                                ridge = sscovr_alpha_ridge$ridge, lasso = rep(0,cond$components), nrFolds = 10,
                                R = cond$components, cdstructure = cd_range[i,], 
                                alpha = sscovr_alpha_ridge$alpha, 
                                inits = "rational", nrstart = 1, 
                                blockcols = c(J/2, J/2), 
                                seed = cd_seed[i], lambda_y = 1, stop_value = 1e-7)
      
      cd_cve[i,] <- c(cv_i$cve, cd_range[i,])
      cd_se[i,] <- c(cv_i$se, cd_range[i,] )
      
      print(rep(i, 10))
    }
    
    # for common-distinctive, i just go with the smallest cv error
    cd_chosen <- matrix(cd_range[which.min(cd_cve[,2]),], nrow = 1)
    
  }
  
  if (cond$components == 3){
    cd_chosen <- matrix(c(1,2,0), nrow = 1)
  }
  
  sscovr_cd <- cd_chosen
  
  # lasso-finding #
  
  # for estimating 2 components, 
  # the number of specified zeros is the total number of defined zeros * 1/3
  
  if (cond$components == 2){
    zeros <- rep(floor(sum(abs(dat$Px) < 1e-9) /3), 2)
  }
  
  if (cond$components == 3){
    zeros <- apply(abs(dat$Px) < 1e-9, 2, sum)
  }
  
  sscovrlasso <- findLasso_sscovr(X = X_train, Y = y_train
                                  , zeros = zeros
                                  , R = cond$components, 
                                  init = 50, ridge = sscovr_alpha_ridge$ridge
                                  , blockcols = c(J/2, J/2)
                                  , cdstructure = sscovr_cd
                                  , alpha = sscovr_alpha_ridge$alpha, 
                                  maxiterOut = 1, maxiterIn = 30)
  
  sscovr1 <- sscovr2(X = X_train, Y = y_train,
                     blockcols = c(J/2, J/2), 
                     R = cond$components, cdstructure = sscovr_cd,
                     alpha = sscovr_alpha_ridge$alpha, 
                     lambda1 = sscovrlasso$lasso, 
                     lambda2 = sscovr_alpha_ridge$ridge, 
                     inits = "rational", nrstart = 1, MAXITER = 10000, 
                     stop_value = 1e-7, lambda_y = 1, ridge_in_between = T, ridge_end = T)
  
  sscovr_nzeros <- apply(abs(sscovr1$cd_results[[1]]$W) < 1e-9, 2, sum)
  
  sum((X_test %*% sscovr1$cd_results[[1]]$W %*% sscovr1$cd_results[[1]]$P[1,] - y_test)^2) / sum(y_test^2)
  sum((X_train %*% sscovr1$cd_results[[1]]$W %*% sscovr1$cd_results[[1]]$P[1,] - y_train)^2) / sum(y_train^2)
  
  # spcovr # 
  # same alpha can be used for spcovr, because
  # when alpha cv is done, cd structure is assumed to be all common
  
  # lasso-finding #
  
  # for estimating 2 components, 
  # the number of specified zeros is the total number of defined zeros * 1/3
  
  if (cond$components == 2){
    zeros <- rep(floor(sum(abs(dat$Px) < 1e-9) /3), 2)
  }
  
  if (cond$components == 3){
    zeros <- apply(abs(dat$Px) < 1e-9, 2, sum)
  }
  
  spcovrlasso <- findLasso_sscovr(X = X_train, Y = y_train
                                  , zeros = zeros
                                  , R = cond$components, 
                                  init = 50, ridge = sscovr_alpha_ridge$ridge
                                  , blockcols = c(J/2, J/2)
                                  , cdstructure = matrix(rep(0, cond$components), nrow = 1)
                                  , alpha = sscovr_alpha_ridge$alpha, 
                                  maxiterOut = 1, maxiterIn = 30)
  
  spcovr1 <- sscovr2(X = X_train, Y = y_train,
                     blockcols = c(J/2, J/2), 
                     R = cond$components, 
                     cdstructure = matrix(rep(0,cond$components), nrow = 1),
                     alpha = sscovr_alpha_ridge$alpha, 
                     lambda1 = spcovrlasso$lasso, 
                     lambda2 = sscovr_alpha_ridge$ridge, 
                     inits = "rational", nrstart = 1, MAXITER = 10000, 
                     stop_value = 1e-7, lambda_y = 1, ridge_in_between = T, ridge_end = T)
  
  spcovr_nzeros <- apply(abs(spcovr1$cd_results[[1]]$W) < 1e-9, 2, sum)
  
  sum((X_test %*% spcovr1$cd_results[[1]]$W %*% spcovr1$cd_results[[1]]$P[1,] - y_test)^2) / sum(y_test^2)
  sum((X_train %*% spcovr1$cd_results[[1]]$W %*% spcovr1$cd_results[[1]]$P[1,] - y_train)^2) / sum(y_train^2)
  
  # pcr (scad) #
  
  ridge_range <-  c(100, 50, 30, 10, 5, 1, 0.5, 0.1)
  
  ridge_seed <- sample(length(ridge_range))
  
  ridge_cve <- matrix(NA, ncol = 1+1, nrow = length(ridge_range))
  ridge_se <- matrix(NA, ncol = 1+1, nrow = length(ridge_range))
  
  fixW <- matrix(1, ncol = cond$components, nrow = J)
  
  for (i in 1:length(ridge_range)){
    
    cv_i <- scad_cv(X = X_train, R = cond$components, ridge = ridge_range[i], 
                    lasso = rep(0, cond$components), nrFolds = 10, 
                    fixW = fixW, MAXITER = 10000, stop_value = 1e-7, seed = ridge_seed[i])
    
    ridge_cve[i,] <- c(cv_i$cve, ridge_range[i])
    ridge_se[i,] <- c(cv_i$se, ridge_range[i])
    
    print(rep(i, 10))
  }
  
  # minimum alpha + 1se value:
  ridge_1se <- min(ridge_cve[,1]) + ridge_se[which.min(ridge_cve[,1]),1]
  
  # smallest alpha value, from the cv errors smaller than 1se
  ridge_chosen <- matrix(ridge_cve[ridge_cve[,1] < ridge_1se,], ncol = 2)[1,2]
  
  scad_ridge <- ridge_chosen
  
  # cross validation for common-distinctive #
  if(cond$components == 3){
    
    fixW <- matrix(0, nrow = J, ncol = cond$components)
    fixW[1:(J/2), 1] <- 1
    fixW[(J/2 + 1):J, 2] <- 1
    fixW[,3:cond$components] <- 1
    
    scad_cd <- matrix(c(1,2,0), nrow = 1)
    
  }
  
  if(cond$components == 2){
    cd_range <- matrix(c(1,2,
                         1,0,
                         2,0,
                         0,0), ncol = 2, byrow = T)
    
    cd_seed <- sample(1:1000, nrow(cd_range))
    
    cd_cve <- matrix(NA, ncol = 1+2, nrow = nrow(cd_range))
    cd_se <- matrix(NA, ncol = 1+2, nrow = nrow(cd_range))
    
    for (i in 1:nrow(cd_range)){
      fixW <- matrix(0, nrow = J, ncol = cond$components)
      
      if (sum(cd_range[i,] == c(1,2)) == 2){
        fixW[1:(J/2), 1] <- 1
        fixW[(J/2 + 1):J, 2] <- 1
      }
      
      if (sum(cd_range[i,] == c(1,0)) == 2){
        fixW[1:(J/2), 1] <- 1
      }
      
      if (sum(cd_range[i,] == c(2,0)) == 2){
        fixW[(J/2 + 1):J, 1] <- 2
      }
      
      cv_i <- scad_cv(X = X_train, R = cond$components, ridge = scad_ridge, 
                      lasso = rep(0, cond$components), nrFolds = 10, 
                      fixW = fixW, MAXITER = 10000, stop_value = 1e-7, seed = cd_seed[i])
      
      cd_cve[i,] <- c(cv_i$cve, cd_range[i,])
      cd_se[i,] <- c(cv_i$se, cd_range[i,] )
      
      print(rep(i, 10))
    }
    
    # for common-distinctive, i just go with the smallest cv error
    cd_chosen <- matrix(cd_range[which.min(cd_cve[,1]),], nrow = 1)
    
    scad_cd <- cd_chosen
    
    fixW <- matrix(0, nrow = J, ncol = cond$components)
    
    if (sum(scad_cd == c(1,2)) == 2){
      fixW[1:(J/2), 1] <- 1
      fixW[(J/2 + 1):J, 2] <- 1
    }
    
    if (sum(scad_cd == c(1,0)) == 2){
      fixW[1:(J/2), 1] <- 1
    }
    
    if (sum(scad_cd == c(2,0)) == 2){
      fixW[(J/2 + 1):J, 1] <- 2
    }
    
  }
  
  # lasso finding #
  scadlasso <- findLasso_SCaDs(dat = X_train, 
                               zeros = apply(abs(dat$Px[,1:cond$components]) < 1e-9, 2, sum),
                               R = cond$components, init = 50, ridge = scad_ridge, 
                               maxiterOut = 2, maxiterIn = 50, fixW = fixW)
  
  # estimation #
  scad1 <- sparseSCAcpp(X = X_train, Q = cond$components, 
                        RIDGE = scad_ridge, LASSO = scadlasso$lasso, 
                        fixW = fixW, maxItrOuterloop = 100000, nStarts = 1,
                        print = FALSE, tol = 10^-8)
  
  pcr_nzeros <- apply(abs(scad1$W) < 1e-9, 2, sum)
  
  scadtmat <- X_train %*% scad1$W
  scadreg <- lm(y_train ~ scadtmat)
  scadreg$coefficients[is.na(scadreg$coefficients)] <- 0
  
  sum((X_test %*% scad1$W %*% scadreg$coefficients[-1] - y_test)^2) / sum(y_test^2)
  sum((X_train %*% scad1$W %*% scadreg$coefficients[-1] - y_train)^2) / sum(y_train^2)
  
  # PLS #
  
  A = list(X_train,y_train)
  
  C = matrix(c(0, 1, 1, 0), 2, 2)
  
  ssq1<-sum(X_train^2)
  ssq2<-sum(y_train^2)
  
  sgcca1 <- RGCCA::sgcca(A,C, c1=c(0.42,1),ncomp=c(cond$components,1), scheme="factorial",verbose=TRUE)
  #note that RGCCA scales all variables to unit var
  
  s1<-sgcca1$AVE[[1]]
  sum(s1[[1]])#VAF by the two components in pred. data reported in Table 1
  sum(s1[[2]])#VAF by the two components in outcome reported in Table 1
  
  sgcca_W <- sgcca1$astar[[1]]
  sgccatmat <- X_train %*% sgcca_W
  
  sgccareg <- lm(y_train ~ sgccatmat)$coefficients#find optimal regression weights to use for out-of-sample prediction
  
  pls_nzeros <- apply(abs(sgcca1$astar[[1]]) < 1e-9, 2, sum)
  
  sum((X_test %*% sgcca_W %*% sgccareg[-1] - y_test)^2) / sum(y_test^2)
  sum((X_train %*% sgcca_W %*% sgccareg[-1] - y_train)^2) / sum(y_train^2)
  
  print(rep(rrr, 20))
  
  time2 <- Sys.time()
  
  time_taken <- time2 - time1
  
  save(dat, train_index, X_train, y_train, X_test, y_test,
       sscovr1, sscovr_alpha_ridge, sscovr_cd, sscovrlasso,
       spcovr1, spcovrlasso, 
       scad1, scadlasso, scad_cd, scad_ridge, scadreg, fixW, 
       sgcca1, sgccareg, cond, time_taken,
       file = paste("../simulation study/results/results2/results_25_feb_2020_", rrr, ".Rdata",sep = ""))
  
  rm(list=setdiff(ls(), needed_objects))
  
  flush.console()
  
}
