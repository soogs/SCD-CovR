# functions for the evaluation criteria # 

# dec-16 addition:
# number of distinctive components
# (i simply observe if this weights column indicates a distinctive component)

zero <- function(estimate, defined){
  nzeros <- sum(abs(estimate) < 1e-7)
  
  tucker <- RegularizedSCA::TuckerCoef(defined, estimate)
  
  estimate <- estimate[,tucker$perm]
  
  ratio <- sum((abs(estimate) < 1e-7) + 
                 (abs(defined) < 1e-7) == 2) / (nzeros)
  return(ratio)
}

nonzero <- function(estimate, defined){
  nonzeros <-  sum(abs(estimate) > 1e-7)
  
  tucker <- RegularizedSCA::TuckerCoef(defined, estimate)
  
  estimate <- estimate[,tucker$perm]
  
  ratio <- sum((abs(estimate) > 1e-7) + 
                 (abs(defined) > 1e-7) == 2) / (nonzeros)
  return(ratio)
}

corrects <- function(estimate, defined){
  total <- prod(dim(estimate))
  
  tucker <- RegularizedSCA::TuckerCoef(defined, estimate)
  
  estimate <- estimate[,tucker$perm]
  
  ratio <- (sum((abs(estimate) > 1e-7) + (abs(defined) > 1e-7) == 2) + 
              sum((abs(estimate) < 1e-7) + (abs(defined) < 1e-7) == 2)) / (total)
  
  return(ratio)
}

y_error <- function(X, W, Py, y){
  result <- sum((X %*% W %*% Py - y)^2) / (sum(y^2))
  return (result)
}

X_error <- function(X, W, Px){
  result <- sum((X %*% W %*% Px - X)^2) / (sum(X^2))
  return(result)
}

# nonzero_D function calculates the proportion of =
# number of nonzero elements in the corresponding datablock / total number of nonzero elements
# compares this to the first two defined coefficients
nonzero_D <- function(estimate, defined){
  tucker <- RegularizedSCA::TuckerCoef(defined, estimate)
  
  J <- nrow(defined)
  
  # using the defined weights matrix,
  # take a sum of the first or the second half of the weights coefficients,
  # which of the column sums equal to 0?
  # if the sum of the first half equals 0, then the component is distinctive to block 2
  
  d1 <- which(apply(abs(defined[(J/2 + 1):J,]), 2, sum) < 1e-7) 
  d2 <- which(apply(abs(defined[1:(J/2),]), 2, sum) < 1e-7)
  
  estimate <- estimate[,tucker$perm]
  # arrange the estimated weights in the same order as the defined matrix
  
  nonzero_d1 <- abs(estimate[,d1]) > 1e-7
  nonzero_d2 <- abs(estimate[,d2]) > 1e-7
  
  ratio_d1 <- sum(nonzero_d1[1:(J/2)]) / sum(nonzero_d1)
  ratio_d2 <- sum(nonzero_d2[(J/2 + 1):J]) / sum(nonzero_d2)
  
  ratio <- mean(c(ratio_d1, ratio_d2))
  
  ratio_list <- list(ratio = ratio, 
                     ratio_d1 = ratio_d1, ratio_d2 = ratio_d2)
  
  return(ratio_list)
}

t_type <- function(T_estimate, T_defined){
  
  tucker <- RegularizedSCA::TuckerCoef(T_defined, T_estimate)
  
  T_estimate <- T_estimate[,tucker$perm]
  
  type_list <- list(t1 = data.frame(type = rep(0,4), value = rep(0,4)), 
                    t2 = data.frame(type = rep(0,4), value = rep(0,4)),
                    t3 = data.frame(type = rep(0,4), value = rep(0,4)))
  
  for (i in 1:3){
    intact_result <- intact_tucker(t_vec = T_estimate[,i], T_defined = T_defined)
    type_list[[i]]$type[1] <- intact_result$which_intact
    type_list[[i]]$value[1] <- intact_result$tucker
  }
  
  for (i in 1:3){
    merge_result <- merge_tucker(t_vec = T_estimate[,i], T_defined = T_defined)
    type_list[[i]]$type[2] <- merge_result$which_merge
    type_list[[i]]$value[2] <- merge_result$tucker
  }
  
  split_result <- split_tucker(T_estimate = T_estimate, T_defined = T_defined)
  
  for (i in 1:3){
    type_list[[i]]$type[3] <- as.character(split_result[[i]]$type)
    type_list[[i]]$value[3] <- split_result[[i]]$value
  }
  
  return(type_list)
}

# check tucker congruence of two vectors
congruence <- function(a, b){
  
  if(anyNA(cbind(a, b))){
    return (0)
  }
  
  if (sum(abs(a)) < 1e-10 | sum(abs(b)) < 1e-10){
    result <- 0
    
    return(result)
  }
  
  result <- (t(a) %*% b) / sqrt(sum(a^2) * sum(b^2))
  return(c(result))
}


merge_tucker <- function(t_vec, T_defined){
  
  merge12 <- congruence(t_vec, T_defined[,1] + T_defined[,2])
  merge13 <- congruence(t_vec, T_defined[,1] + T_defined[,3])
  merge23 <- congruence(t_vec, T_defined[,2] + T_defined[,3])
  merge123 <- congruence(t_vec, T_defined[,1] + T_defined[,2] + T_defined[,3])
  
  merge_all <- c(merge12, merge13, merge23, merge123)
  merge_names <- c("merge12", "merge13", "merge23", "merge123")
  
  max_merge <- merge_all[which.max(abs(merge_all))]
  which_merge <- merge_names[which.max(abs(merge_all))]
  
  return(list(which_merge = which_merge, tucker = max_merge))
}


split_tucker <- function(T_estimate, T_defined){
  
  tucker <- RegularizedSCA::TuckerCoef(T_defined, T_estimate)
  
  T_estimate <- T_estimate[,tucker$perm]
  
  T_estimate[,colSums(T_estimate == 0) == nrow(T_estimate)] <- NA
  
  split1_12 <- congruence(T_estimate[,1] + T_estimate[,2], T_defined[,1])
  split1_13 <- congruence(T_estimate[,1] + T_estimate[,3], T_defined[,1])
  split1_23 <- congruence(T_estimate[,2] + T_estimate[,3], T_defined[,1])
  split1_123 <- congruence(rowSums(T_estimate), T_defined[,1])
  
  split2_12 <- congruence(T_estimate[,1] + T_estimate[,2], T_defined[,2])
  split2_13 <- congruence(T_estimate[,1] + T_estimate[,3], T_defined[,2])
  split2_23 <- congruence(T_estimate[,2] + T_estimate[,3], T_defined[,2])
  split2_123 <- congruence(rowSums(T_estimate), T_defined[,2])
  
  split3_12 <- congruence(T_estimate[,1] + T_estimate[,2], T_defined[,3])
  split3_13 <- congruence(T_estimate[,1] + T_estimate[,3], T_defined[,3])
  split3_23 <- congruence(T_estimate[,2] + T_estimate[,3], T_defined[,3])
  split3_123 <- congruence(rowSums(T_estimate), T_defined[,3])
  
  split_all <- c(split1_12, split1_13, split1_23, split1_123,
                 split2_12, split2_13, split2_23, split2_123,
                 split3_12, split3_13, split3_23, split3_123)
  
  split_names <- c("split1_12", "split1_13", "split1_23", "split1_123",
                   "split2_12", "split2_13", "split2_23", "split2_123",
                   "split3_12", "split3_13", "split3_23", "split3_123")
  
  split_results <- data.frame(names = split_names, tucker = split_all)
  
  # component 1 #
  est1_type <- c("split1_12", "split1_13", "split1_123",
                "split2_12", "split2_13", "split2_123",
                "split3_12", "split3_13", "split3_123")
              
  est1_value <- c(split1_12, split1_13, split1_123,
                 split2_12, split2_13, split2_123,
                 split3_12, split3_13, split3_123)
  
  est1 <- data.frame(type = est1_type, value = est1_value)
  
  est1_return <- est1[which.max(abs(est1$value)),]
  
  # component 2 #
  est2_type <- c("split1_12", "split1_23", "split1_123",
                 "split2_12", "split2_23", "split2_123",
                 "split3_12", "split3_23", "split3_123")
  
  est2_value <- c(split1_12, split1_23, split1_123,
                  split2_12, split2_23, split2_123,
                  split3_12, split3_23, split3_123)
  
  est2 <- data.frame(type = est2_type, value = est2_value)
  
  est2_return <- est2[which.max(abs(est2$value)),]
  
  
  # component 2 #
  est3_type <- c("split1_13", "split1_23", "split1_123",
                 "split2_13", "split2_23", "split2_123",
                 "split3_13", "split3_23", "split3_123")
  
  est3_value <- c(split1_13, split1_23, split1_123,
                  split2_13, split2_23, split2_123,
                  split3_13, split3_23, split3_123)
  
  est3 <- data.frame(type = est3_type, value = est3_value)
  
  est3_return <- est3[which.max(abs(est3$value)),]
  
  return_list <- list(t1 = est1_return,
                      t2 = est2_return,
                      t3 = est3_return)
  
  return(return_list)
}


intact_tucker <- function(t_vec, T_defined){
  
  intact1 <- congruence(t_vec, T_defined[,1])
  intact2 <- congruence(t_vec, T_defined[,2])
  intact3 <- congruence(t_vec, T_defined[,3])
  
  intact_all <- c(intact1, intact2, intact3)
  intact_names <- c("intact1", "intact2", "intact3")
  
  max_intact <- intact_all[which.max(abs(intact_all))]
  which_intact <- intact_names[which.max(abs(intact_all))]
  
  return(list(which_intact = which_intact, tucker = max_intact))
}

distinctive_quality <- function(t_type_list, W_estimate, T_estimate, X_train){
  
  J <- nrow(W_estimate)
  I <- nrow(T_estimate)
  
  # RegularizedSCA::TuckerCoef function cannot handle vectors filled only with 0
  # so i fill them up with 1
  if(ncol(W_estimate) != 3){
    T_estimate[,colSums(abs(T_estimate) == 0) == I] <- 1
  }
  
  if(ncol(W_estimate != 3)){
    W_estimate <- cbind(W_estimate, 1)
  }
  
  blockcols <- list(c(1:(J/2)), c((J/2+1):J), c(1:J))
  
  tucker <- RegularizedSCA::TuckerCoef(T_estimate, X_train %*% W_estimate)
  
  W_estimate <- W_estimate[,tucker$perm]
  
  types <- c(
    t_type_list$t1$type[which.max(abs(t_type_list$t1$value))],
    t_type_list$t2$type[which.max(abs(t_type_list$t2$value))],
    t_type_list$t3$type[which.max(abs(t_type_list$t3$value))]
  )
  
  
  distinctive_check <- function(w_vec, type){
    ratio <- 111
    
    if (type == "intact1"){
      ratio <- sum(abs(w_vec[blockcols[[1]]]) > 1e-7) / sum(abs(w_vec) > 1e-7)
    }
    
    if (type == "intact2"){
      ratio <- sum(abs(w_vec[blockcols[[2]]]) > 1e-7) / sum(abs(w_vec) > 1e-7)
    }
    
    if (sum(w_vec == 1) == length(w_vec)){
      ratio <- 111
    }
    
    return (ratio)
  }
  
  
  distinctive_ratio <- c()
  
  for (i in 1:3){
    distinctive_ratio[i] <- distinctive_check(w_vec = W_estimate[,i], type = types[i])
  }
  
  return(distinctive_ratio)
  
}



perf <- function(W_estimate, W_defined, 
                 Py, T_defined,
                 X_train, X_true, X_test, 
                 y_train, y_true, y_test){
  
  tucker <- RegularizedSCA::TuckerCoef(W_defined, W_estimate)
  
  W_estimate <- W_estimate[,tucker$perm]
  Py <- Py[tucker$perm]
  
  # first, only the criteria that examine the weights matrix
  zerohits <- zero(W_estimate, W_defined)
  nonzerohits <- nonzero(W_estimate, W_defined)
  correcthits <- corrects(W_estimate, W_defined)
  d_info <- nonzero_D(W_estimate, W_defined)
  d_overall <- d_info$ratio
  ratio_d1 <- d_info$ratio_d1
  ratio_d2 <- d_info$ratio_d2
  w_tucker <- tucker$tucker_value
  
  # errors on y
  fit <- y_error(X = X_train, W = W_estimate, Py = Py, y = y_train)
  pred_test <- y_error(X = X_test, W = W_estimate, Py = Py, y = y_test)
  pred_true <- y_error(X = X_train, W = W_estimate, Py = Py, y = y_true)
  
  # components examination
  # first, i flip the component scores by tucker congruence with the defined T
  # it is possible that the order is different from the orders of W. beware.
  T_estimate <- X_train %*% W_estimate
  
  T_tucker <- RegularizedSCA::TuckerCoef(T_defined, T_estimate)
  
  T_estimate <- T_estimate[,T_tucker$perm]
  
  T_type <- t_type(T_estimate = T_estimate, T_defined = T_defined)
  # the intact and merged are only determined through the component scores,
  # not examining the weights
  
  T_type <- split_check(T_type)
  
  T_quality <- distinctive_quality(t_type_list = T_type, W_estimate = W_estimate, 
                                   T_estimate = T_estimate, X_train = X_train)
  # at this point, if you had somesthing like
  # 0.5 for "intact2", 
  # this means that the weights column USED to calculate the component
  # that is identified with intact2 has 50% of the nonzero elements 
  # in the second block, as opposed to all the nonzero elements.
  
  for (i in 1:3){
    T_type[[i]]$type[4] <- "prop"
    T_type[[i]]$value[4] <- T_quality[i]
  }
  
  Ws <- data.frame(zerohits = zerohits, nonzerohits = nonzerohits, 
                   correcthits = correcthits, d_overall = d_overall, 
                   d1 = ratio_d1, d2 = ratio_d2,
                   w_tucker = w_tucker)
  
  Errors <- data.frame(fit = fit, test = pred_test, true = pred_true)
  
  returnobjects <- list(W = Ws, Error = Errors, T_result = T_type)
  
  return(returnobjects)
}



nonzero_D2 <- function(estimate, defined, compare_with){
  
  J <- nrow(estimate)
  
  blockcols <- list(c(1:(J/2)), c((J/2+1):J), c(1:J))
  
  ratios <- c()
  for (i in 1:2){
    ratio <- sum(abs(estimate[blockcols[[compare_with[i]]],i]) > 1e-7) / sum(abs(estimate[,i]) > 1e-7) 
    
    ratios[i] <- ratio
    
  }
  
  ratio_c <- 0
  ratio_d1 <- 0
  ratio_d2 <- 0
  
  c_index <- which(compare_with == 3)
  d1_index <- which(compare_with == 1)
  d2_index <- which(compare_with == 2)
  
  ratio_c <- ratios[c_index]
  ratio_d1 <- ratios[d1_index]
  ratio_d2 <- ratios[d2_index]
  
  if (sum(compare_with == 1) == 0){
    ratio_d1 <- 0
  } 
  
  if (sum(compare_with == 2) == 0){
    ratio_d2 <- 0
  } 
  
  ratio <- mean(ratios)
  
  ratio_list <- list(ratio = ratio, 
                     ratio_d1 = ratio_d1, ratio_d2 = ratio_d2)
  
  return(ratio_list)
}


# perf2 function for 2 components
perf2 <- function(W_estimate, W_defined, 
                 Py, T_defined,
                 X_train, X_true, X_test, 
                 y_train, y_true, y_test){
  
  congs <- matrix(NA, ncol = 3, nrow = 2)
  
  for (i in 1:2){
    cong1 <- congruence(W_estimate[,i], W_defined[,1])
    cong2 <- congruence(W_estimate[,i], W_defined[,2])
    cong3 <- congruence(W_estimate[,i], W_defined[,3])
    
    congs[i,] <- c(cong1, cong2, cong3)
  }
  
  compare_with <- apply(abs(congs), 1, which.max)
  
  # if the two estimated components have equal congruence 
  # with the defined components,
  # let the first estimated component be paired with the next best defined
  if (compare_with[1] == compare_with[2]){
    compare_with[1] <- order(congs[1,], decreasing = T)[2]
  }
  
  W_defined <- W_defined[, compare_with]
  
  tucker <- RegularizedSCA::TuckerCoef(W_defined, W_estimate)
  
  W_estimate <- W_estimate[,tucker$perm]
  Py <- Py[tucker$perm]
  
  # first, only the criteria that examine the weights matrix
  zerohits <- zero(W_estimate, W_defined)
  nonzerohits <- nonzero(W_estimate, W_defined)
  correcthits <- corrects(W_estimate, W_defined)
  d_info <- nonzero_D2(W_estimate, W_defined, compare_with = compare_with)
  d_overall <- d_info$ratio
  ratio_d1 <- d_info$ratio_d1
  ratio_d2 <- d_info$ratio_d2
  w_tucker <- tucker$tucker_value
  
  # errors on y
  fit <- y_error(X = X_train, W = W_estimate, Py = Py, y = y_train)
  pred_test <- y_error(X = X_test, W = W_estimate, Py = Py, y = y_test)
  pred_true <- y_error(X = X_train, W = W_estimate, Py = Py, y = y_true)
  
  # components examination
  T_estimate <- X_train %*% W_estimate
  
  T_estimate <- cbind(T_estimate, 0)
    
  T_tucker <- RegularizedSCA::TuckerCoef(T_defined, T_estimate)
  
  T_estimate <- T_estimate[,T_tucker$perm]
  
  T_type <- t_type2(T_estimate = T_estimate, T_defined = T_defined)
  # the intact and merged are only determined through the component scores,
  # not examining the weights
  
  T_type <- split_check(T_type)
  
  T_quality <- distinctive_quality(t_type_list = T_type, W_estimate = W_estimate,
                                   T_estimate = T_estimate, X_train = X_train)
  # at this point, if you had something like
  # 0.5 for "intact2", 
  # this means that the weights column USED to calculate the component
  # that is identified with intact2 has 50% of the nonzero elements 
  # in the second block, as opposed to all the nonzero elements.
  
  for (i in 1:3){
    T_type[[i]]$type[4] <- "prop"
    T_type[[i]]$value[4] <- T_quality[i]
  }
  
  Ws <- data.frame(zerohits = zerohits, nonzerohits = nonzerohits, 
                   correcthits = correcthits, d_overall = d_overall, 
                   d1 = ratio_d1, d2 = ratio_d2,
                   w_tucker = w_tucker)
  
  Errors <- data.frame(fit = fit, test = pred_test, true = pred_true)
  
  returnobjects <- list(W = Ws, Error = Errors, T_result = T_type, compare_with = compare_with)
  
  return(returnobjects)
}





t_type2 <- function(T_estimate, T_defined){
  
  tucker <- RegularizedSCA::TuckerCoef(T_defined, T_estimate)
  
  T_estimate <- T_estimate[,tucker$perm]
  
  type_list <- list(t1 = data.frame(type = rep(0,4), value = rep(0,4)), 
                    t2 = data.frame(type = rep(0,4), value = rep(0,4)),
                    t3 = data.frame(type = rep(0,4), value = rep(0,4)))
  
  for (i in 1:3){
    intact_result <- intact_tucker(t_vec = T_estimate[,i], T_defined = T_defined)
    type_list[[i]]$type[1] <- intact_result$which_intact
    type_list[[i]]$value[1] <- intact_result$tucker
  }
  
  for (i in 1:3){
    merge_result <- merge_tucker(t_vec = T_estimate[,i], T_defined = T_defined)
    type_list[[i]]$type[2] <- merge_result$which_merge
    type_list[[i]]$value[2] <- merge_result$tucker
  }
  
  split_result <- split_tucker(T_estimate = T_estimate, T_defined = T_defined)
  
  for (i in 1:3){
    type_list[[i]]$type[3] <- as.character(split_result[[i]]$type)
    type_list[[i]]$value[3] <- split_result[[i]]$value
  }
  
  return(type_list)
}


# function which correctly filters out splitting components
# if there is only one splitting component
# or if there are two or three but they don't match
# the tucker congruence values for these components become 0
split_check <- function(T_result){
  
  potential <- lapply(T_result, function(x){x[which.max(abs(x[-4, 2])),]})
  
  typecheck <- unlist(lapply(potential, function(x){x[,1]}))
  
  split_index <- grep(pattern = "split", x = typecheck)
  # index of the element that has the word "split" included
  
  # if only one retrieved component is categorized as a split component:
  if (length(split_index) == 1){
    T_result[[split_index]]$value[3] <- 0
  }
  
  # if two components are categorized as split components
  if (length(split_index) == 2){
    
    if(sum(duplicated(typecheck[split_index])) == 0){
      # if the splitting categories do not match:
      # let the congruence of these splitting components become 0 
      for (i in split_index){
        T_result[[i]]$value[3] <- 0
      }
    }
  }
  
  # if three components are categorized as split components
  if (length(split_index) == 3){
    
    if(sum(duplicated(typecheck[split_index])) == 0){
      # if none of the splitting categories match:
      # let the congruence of these splitting components become 0 
      for (i in split_index){
        T_result[[i]]$value[3] <- 0
      }
    }
    
    # if only one pair has a match
    if(sum(duplicated(typecheck[split_index])) == 1){
      which_match <- typecheck == typecheck[duplicated(typecheck[split_index])]
      no_match <- which(!which_match)
      
      T_result[[no_match]]$value[3] <- 0
      
    }
  }
  
  return(T_result)
}

distinctive <- function(estimate){
  J <- nrow(estimate)
  
  # if this is zero, that column is distinctive to block 2
  d2 <- colSums(abs(estimate[1:(J/2),]) > 1e-7)
  
  # if this is zero, that column is distinctive to block 1
  d1 <- colSums(abs(estimate[(J/2 + 1):(J),]) > 1e-7)
  
  # if the entire column is zero, do not count this column in
  # (this is considered as neither distinct nor common)
  no_index <- which(colSums(abs(estimate) > 1e-7) == 0)
  
  common <- as.numeric(((d1 != 0) + (d2 != 0)) > 1)
  
  d2 <- as.numeric(d2 == 0)
  d1 <- as.numeric(d1 == 0)
  
  d2[no_index] <- 0
  d1[no_index] <- 0
  common[no_index] <- 0
  
  result <- data.frame(d1 = d1, d2 = d2, common = common)
  
  return(result)
}
  
