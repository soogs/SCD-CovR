# obtaining evaluation criteria #

# i have these 5 criteria:

# 1. Prediction errors: "pred"
# 2. Correct classification rate: "correcthits"
# 3. Component identification: maximum value from the "type"
# 4. Component scores congruence: mean of the maximum tucker values
# 5. False positive rate (distinctive component): 
# mean of the "props" (always take account of how many distinctives the method found)

source("C:\\Users\\park\\Desktop\\paper2 writing for now\\oct02_new_simulation\\entire_results\\soogeun\\results_extract\\evaluation_criteria.R")

pred_df <- data.frame(ssc = rep(111, 3200), spc = rep(111, 3200),
                      pcr = rep(111, 3200), pls = rep(111, 3200))

hits_df <- data.frame(ssc = rep(111, 3200), spc = rep(111, 3200),
                      pcr = rep(111, 3200), pls = rep(111, 3200))

truepo_df <- data.frame(ssc = rep(111, 3200), spc = rep(111, 3200),
                       pcr = rep(111, 3200), pls = rep(111, 3200))

tucker_df <- data.frame(ssc = rep(111, 3200), spc = rep(111, 3200),
                        pcr = rep(111, 3200), pls = rep(111, 3200))

w_tucker_df <-  data.frame(ssc = rep(111, 3200), spc = rep(111, 3200),
                           pcr = rep(111, 3200), pls = rep(111, 3200))

w_truepo_df <-  data.frame(ssc = rep(111, 3200), spc = rep(111, 3200),
                           pcr = rep(111, 3200), pls = rep(111, 3200))

w_d1_truepo_df <-  data.frame(ssc = rep(111, 3200), spc = rep(111, 3200),
                           pcr = rep(111, 3200), pls = rep(111, 3200))

w_d2_truepo_df <-  data.frame(ssc = rep(111, 3200), spc = rep(111, 3200),
                              pcr = rep(111, 3200), pls = rep(111, 3200))


type_list <- list(ssc = data.frame(t1 = rep(111, 3200), t2 = rep(111, 3200),
                                   t3 = rep(111, 3200)),
                  spc = data.frame(t1 = rep(111, 3200), t2 = rep(111, 3200),
                                   t3 = rep(111, 3200)),
                  pcr = data.frame(t1 = rep(111, 3200), t2 = rep(111, 3200),
                                   t3 = rep(111, 3200)),
                  pls = data.frame(t1 = rep(111, 3200), t2 = rep(111, 3200),
                                   t3 = rep(111, 3200)))

truepo_list <- list(ssc = data.frame(t1 = rep(111, 3200), t2 = rep(111, 3200),
                                   t3 = rep(111, 3200)),
                  spc = data.frame(t1 = rep(111, 3200), t2 = rep(111, 3200),
                                   t3 = rep(111, 3200)),
                  pcr = data.frame(t1 = rep(111, 3200), t2 = rep(111, 3200),
                                   t3 = rep(111, 3200)),
                  pls = data.frame(t1 = rep(111, 3200), t2 = rep(111, 3200),
                                   t3 = rep(111, 3200)))


alpha_ridge_df <- data.frame(ssc_alpha = rep(111, 3200), ssc_ridge = rep(111, 3200),
                             pcr_ridge = rep(111, 3200))

ssc_lasso_df <- data.frame(w1 = rep(111, 3200), w2 = rep(111, 3200), w3 = rep(111, 3200))

pcr_lasso_df <- data.frame(w1 = rep(111, 3200), w2 = rep(111, 3200), w3 = rep(111, 3200))

spc_lasso_df <- data.frame(w1 = rep(111, 3200), w2 = rep(111, 3200), w3 = rep(111, 3200))

ssc_cd_df <- data.frame(w1 = rep(111, 3200), w2 = rep(111, 3200), w3 = rep(111, 3200))

pcr_cd_df <- data.frame(w1 = rep(111, 3200), w2 = rep(111, 3200), w3 = rep(111, 3200))


cond_df <- data.frame(components = rep(111, 3200), dimension = rep(111, 3200), 
                      vafx = rep(111, 3200), vafy = rep(111, 3200), 
                      weak = rep(111, 3200), relevant = rep(111, 3200), 
                      reps = rep(111, 3200))

# loaded <- ls()

# rm(list = loaded)

for (rrr in 1:3200){

  file_name <- paste("C:\\Users\\park\\Desktop\\paper2 writing for now\\oct02_new_simulation\\entire_results\\soogeun\\results_ridge_top\\ridge_top_", rrr, ".Rdata", sep = "")
  
  load(file_name)
  
  train_index <- train_index
  
  X_true <- dat$dattrue[-train_index, -1]
  y_true <- dat$dattrue[-train_index, 1]
  
  T_true <- dat$U[train_index, ]
  
  pls_W <- sgcca1$astar[[1]]
  pls_py <- sgccareg[-1]
  pls_T <- X_train %*% pls_W
  
  ssc_W <- sscovr1$cd_results[[1]]$W
  ssc_py <- sscovr1$cd_results[[1]]$P[1,]
  ssc_T <- X_train %*% ssc_W
  
  spc_W <- spcovr1$cd_results[[1]]$W
  spc_py <- spcovr1$cd_results[[1]]$P[1,]
  spc_T <- X_train %*% spc_W
  
  pcr_W <- scad1$W
  pcr_py <- scadreg$coefficients[-1]
  pcr_T <- X_train %*% pcr_W
  
  
  if (cond$components == 2){
    pls <- perf2(W_estimate = pls_W, W_defined = dat$Px, Py = pls_py, T_defined = T_true, X_train = X_train, X_true = X_true, 
          X_test = X_test, y_train = y_train, y_true = y_true,y_test = y_test)
    
    
    ssc <- perf2(W_estimate = ssc_W, W_defined = dat$Px, Py = ssc_py, T_defined = T_true, X_train = X_train, X_true = X_true, 
          X_test = X_test, y_train = y_train, y_true = y_true,y_test = y_test)
    
    
    spc <- perf2(W_estimate = spc_W, W_defined = dat$Px, Py = spc_py, T_defined = T_true, X_train = X_train, X_true = X_true, 
          X_test = X_test, y_train = y_train, y_true = y_true,y_test = y_test)
    
    
    pcr <- perf2(W_estimate = pcr_W, W_defined = dat$Px, Py = pcr_py, T_defined = T_true, X_train = X_train, X_true = X_true, 
          X_test = X_test, y_train = y_train, y_true = y_true,y_test = y_test)
    
  }
  
  
  if (cond$components == 3){
    pls <- perf(W_estimate = pls_W, W_defined = dat$Px, Py = pls_py, T_defined = T_true, X_train = X_train, X_true = X_true, 
                 X_test = X_test, y_train = y_train, y_true = y_true,y_test = y_test)
    
    
    ssc <- perf(W_estimate = ssc_W, W_defined = dat$Px, Py = ssc_py, T_defined = T_true, X_train = X_train, X_true = X_true, 
                 X_test = X_test, y_train = y_train, y_true = y_true,y_test = y_test)
    
    
    spc <- perf(W_estimate = spc_W, W_defined = dat$Px, Py = spc_py, T_defined = T_true, X_train = X_train, X_true = X_true, 
                 X_test = X_test, y_train = y_train, y_true = y_true,y_test = y_test)
    
    
    pcr <- perf(W_estimate = pcr_W, W_defined = dat$Px, Py = pcr_py, T_defined = T_true, X_train = X_train, X_true = X_true, 
                 X_test = X_test, y_train = y_train, y_true = y_true,y_test = y_test)
    
  }
  
  # 1. prediction error #
  pred_df$ssc[rrr] <- ssc$Error$test
  pred_df$spc[rrr] <- spc$Error$test
  pred_df$pcr[rrr] <- pcr$Error$test
  pred_df$pls[rrr] <- pls$Error$test
  
  # 2. classification rate #
  hits_df$ssc[rrr] <- ssc$W$correcthits
  hits_df$spc[rrr] <- spc$W$correcthits
  hits_df$pcr[rrr] <- pcr$W$correcthits
  hits_df$pls[rrr] <- pls$W$correcthits
  
  # 3. type identification #
  ssc$T_result <- split_check(ssc$T_result)
  spc$T_result <- split_check(spc$T_result)
  pcr$T_result <- split_check(pcr$T_result)
  pls$T_result <- split_check(pls$T_result)
  
  ssc_type <- lapply(ssc$T_result, function(x){x[which.max(abs(x[-4, 2])),]})
  spc_type <- lapply(spc$T_result, function(x){x[which.max(abs(x[-4, 2])),]})
  pcr_type <- lapply(pcr$T_result, function(x){x[which.max(abs(x[-4, 2])),]})
  pls_type <- lapply(pls$T_result, function(x){x[which.max(abs(x[-4, 2])),]})
  
  for (i in 1:3){
    
    # if the tucker congruence is smaller than 0.85, 
    # classify the component as noise component
    # if(abs(ssc_type[[i]]$value) < 0.85){
    #   ssc_type[[i]]$type <- "noise"
    # }
    # 
    # if(abs(spc_type[[i]]$value) < 0.85){
    #   spc_type[[i]]$type <- "noise"
    # }
    # 
    # if(abs(pcr_type[[i]]$value) < 0.85){
    #   pcr_type[[i]]$type <- "noise"
    # }
    # 
    # if(abs(pls_type[[i]]$value) < 0.85){
    #   pls_type[[i]]$type <- "noise"
    # }
    
    if(ssc_type[[i]]$value == 0){
      ssc_type[[i]]$type <- "2comp_estimated"
    }
    
    
    if(spc_type[[i]]$value == 0){
      spc_type[[i]]$type <- "2comp_estimated"
    }
    
    
    if(pcr_type[[i]]$value == 0){
      pcr_type[[i]]$type <- "2comp_estimated"
    }
    
    
    if(pls_type[[i]]$value == 0){
      pls_type[[i]]$type <- "2comp_estimated"
    }
    
    
    type_list$ssc[[i]][rrr] <- ssc_type[[i]]$type
    type_list$spc[[i]][rrr] <- spc_type[[i]]$type
    type_list$pcr[[i]][rrr] <- pcr_type[[i]]$type
    type_list$pls[[i]][rrr] <- pls_type[[i]]$type
  }
  
  # 4. tucker congruence #
  for (i in 1:3){
    
    if(ssc_type[[i]]$type == "2comp_estimated"){
      ssc_type[[i]]$value <- NA
    }
    
    
    if(spc_type[[i]]$type == "2comp_estimated"){
      spc_type[[i]]$value <- NA
    }
    
    
    if(pcr_type[[i]]$type == "2comp_estimated"){
      pcr_type[[i]]$value <- NA
    }
    
    
    if(pls_type[[i]]$type == "2comp_estimated"){
      pls_type[[i]]$value <- NA
    }
    
  }
  
  tucker_df$ssc[rrr] <- mean(abs(unlist(lapply(ssc_type, function(x){x[,2]}))), na.rm = TRUE)
  tucker_df$spc[rrr] <- mean(abs(unlist(lapply(spc_type, function(x){x[,2]}))), na.rm = TRUE)
  tucker_df$pcr[rrr] <- mean(abs(unlist(lapply(pcr_type, function(x){x[,2]}))), na.rm = TRUE)
  tucker_df$pls[rrr] <- mean(abs(unlist(lapply(pls_type, function(x){x[,2]}))), na.rm = TRUE)
  
  # 5. true positive rate #
  ssc_truepo <- unlist(lapply(ssc$T_result, function(x){x$value[4]}))
  spc_truepo <- unlist(lapply(spc$T_result, function(x){x$value[4]}))
  pcr_truepo <- unlist(lapply(pcr$T_result, function(x){x$value[4]}))
  pls_truepo <- unlist(lapply(pls$T_result, function(x){x$value[4]}))
  
  # true positive rate per component # 
  truepo_list$ssc[rrr,] <- ssc_truepo
  truepo_list$spc[rrr,] <- spc_truepo
  truepo_list$pcr[rrr,] <- pcr_truepo
  truepo_list$pls[rrr,] <- pls_truepo
  
  # taking the mean #
  truepo_df$ssc[rrr] <- mean(ssc_truepo[ssc_truepo != 111])
  truepo_df$spc[rrr] <- mean(spc_truepo[spc_truepo != 111])
  truepo_df$pcr[rrr] <- mean(pcr_truepo[pcr_truepo != 111])
  truepo_df$pls[rrr] <- mean(pls_truepo[pls_truepo != 111])
  
  truepo_df[rrr, is.na(truepo_df[rrr,])] <- 111
  
  # tuning parameters # 
  alpha_ridge_df$ssc_alpha[rrr] <- sscovr_alpha_ridge$alpha
  alpha_ridge_df$ssc_ridge[rrr] <- sscovr_alpha_ridge$ridge
  alpha_ridge_df$pcr_ridge[rrr] <- scad_ridge
  
  if (cond$components == 2){
    sscovrlasso$lasso <- append(sscovrlasso$lasso, 111)
    spcovrlasso$lasso <- append(spcovrlasso$lasso, 111)
    scadlasso$lasso <- append(scadlasso$lasso, 111)
    
    sscovr_cd <- cbind(sscovr_cd, 111)
    scad_cd <- cbind(scad_cd, 111)
  }
  
  # 6. tucker on W
  w_tucker_df$ssc[rrr] <- ssc$W$w_tucker
  w_tucker_df$spc[rrr] <- spc$W$w_tucker
  w_tucker_df$pcr[rrr] <- pcr$W$w_tucker
  w_tucker_df$pls[rrr] <- pls$W$w_tucker
  
  # 7. true positives based on tucker on W
  w_truepo_df$ssc[rrr] <- ssc$W$d_overall
  w_truepo_df$spc[rrr] <- spc$W$d_overall
  w_truepo_df$pcr[rrr] <- pcr$W$d_overall
  w_truepo_df$pls[rrr] <- pls$W$d_overall
  
  ssc_lasso_df[rrr,] <- sscovrlasso$lasso
  
  spc_lasso_df[rrr,] <- spcovrlasso$lasso
  
  pcr_lasso_df[rrr,] <- scadlasso$lasso
  
  ssc_cd_df[rrr,] <- sscovr_cd
  
  pcr_cd_df[rrr,] <- scad_cd
  
  cond_df[rrr,] <- cond
  
  # 8. distinctive 1 and distinctive 2: individual true positive rate
  w_d1_truepo_df$ssc[rrr] <- ssc$W$d1
  w_d1_truepo_df$spc[rrr] <- spc$W$d1
  w_d1_truepo_df$pcr[rrr] <- pcr$W$d1
  w_d1_truepo_df$pls[rrr] <- pls$W$d1
  
  w_d2_truepo_df$ssc[rrr] <- ssc$W$d2
  w_d2_truepo_df$spc[rrr] <- spc$W$d2
  w_d2_truepo_df$pcr[rrr] <- pcr$W$d2
  w_d2_truepo_df$pls[rrr] <- pls$W$d2
  
  
  print(rrr)
  
  rm(list = loaded)
  
}



