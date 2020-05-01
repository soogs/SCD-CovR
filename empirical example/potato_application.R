
setwd("C:\\Users\\park\\Desktop\\paper2 writing for now\\oct02_new_simulation\\")

setwd("../potato_application/")
potato <- read.table("./potato.txt", sep = ",")

setwd("C:\\Users\\park\\Desktop\\paper2 writing for now\\paper2_multicore_results_complete\\")
Rcpp::sourceCpp("./sparseSCA.cpp")
Rcpp::sourceCpp("./updateW.cpp")

setwd("C:\\Users\\park\\Desktop\\paper2 writing for now\\oct02_new_simulation\\")
source("./evaluation_criteria.R")
source("./sscovr_looking_in.R")
source("./spcovrdata.R")
source("./cv_eigenvector_general.R")
source("./findLasso_sscovr_differentpenalties_adding_tryCatch_to_july10_version.R")
source("./findLasso_SCaDs.R")
source("./scad_cv.R")
source("./cv_eigenvector_sgcca.R")

# setting back the working directory
setwd("C:\\Users\\park\\Desktop\\Rstudio Github - Copy\\Project2_SPCovR\\paper2_writing\\")

library(RegularizedSCA)

var_label <-c(
  "Year",
  "storage",
  "organic",
  "DensityGroup",
  "Density",
  "DMM",
  "PEU",
  "citricacid",
  "starch",
  "TotalN",
  "phytic",
  "Ca",
  "Mg",
  "Na",
  "K",
  "his1",
  "his2",
  "his3",
  "his4",
  "his5",
  "his6",
  "FractureWork20",
  "BreakWork20",
  "stressT20",
  "strainH20",
  "modulus20",
  "slope20",
  "FractureWor100",
  "BreakWork100",
  "stressT100",
  "strainH100",
  "modulus100",
  "slope100",
  "FractureWor250",
  "BreakWork250",
  "stressT250",
  "strainH250",
  "modulus250",
  "slope250",
  "FractureWor500",
  "BreakWork500",
  "stressT500",
  "strainH500",
  "modulus500",
  "slope500",
  "FractureWor750",
  "BreakWork750",
  "stressT750",
  "strainH750",
  "modulus750",
  "slope750",
  "FractureWor1000",
  "BreakWork1000",
  "stressT1000",
  "strainH1000",
  "modulus1000",
  "slope1000",
  "ref",
  "hard",
  "firm",
  "elas",
  "adhes",
  "grainy",
  "mealy",
  "moist",
  "chewi"
)

colnames(potato) <- var_label

varieties <- c('v1',
               'v2',
               'v3',
               'v4',
               'v5',
               'v6',
               'v10',
               'v11',
               'v30',
               'v31',
               'v20',
               'v21',
               'v40',
               'v41',
               'v400',
               'v401',
               'v300',
               'v301',
               'v601',
               'v603',
               'v1',
               'v2',
               'v3',
               'v4',
               'v5',
               'v6')

potato <- cbind(varieties, potato)

dat <- potato[,7:length(potato)]

# so it's also possible to just exclude the citric acid variable:
dat <- dat[,-2]

sum(complete.cases(dat)) # this is now 20

dat <- dat[complete.cases(dat),]
dim(dat)

# pre-processing
# chemical:
colnames(dat)[1:14]
chem <- dat[,1:14]

chem <- pre_process(DATA = chem, weight = T)
sum(chem^2)

# uniaxial:
colnames(dat)[15:50]
uni <- dat[,15:50]

uni <- pre_process(DATA = uni, weight = T)
sum(uni^2)

# sensory:
sens <- dat[,51:59]
sens <- pre_process(DATA = sens, weight = T)

# using the principal component of the sensory variable as the criterion variable
sens_pc <- prcomp(sens)

pc <- sens_pc$x[,1]

sum(pc^2)

X <- cbind(chem, uni)

# first calculate the alpha by maximum likelihood 
# we need the estimates for variance of Ex and ey
# variance of Ex is the percentage of unexplained variance,
# when you do pca on X and choose the number of components
# through scree test

x_pc <- prcomp(X)

#scree plot
plot(x_pc$sdev^2, xlab = "Principal Component",
     ylab = "Proportion of Variance Explained",
     type = "b")


# 3 components explain most of the variance. 

var_Ex <- (sum(x_pc$sdev^2) - sum(x_pc$sdev[1:3]^2)) /  sum(x_pc$sdev^2) 
# proportion of unexplained variance by the dominant components = 22.4%

# for estimation of ey, this is the proportion of unexplained variance when y is regressed on X.
reg <- MASS::lm.ridge(formula = pc ~ X, lambda = 0.001)

var_ey <- sum((X %*% reg$coef - pc)^2) / sum(pc^2) # proportion of unexplained variance (1 - R^2)
# = 56.9%

# combining them together:
a_ml <- 1 -sum(X^2) / (sum(X^2) + sum(pc^2) * var_Ex / var_ey)  # 1 - alpha = 0.866
# therefore a lot of weight on X

Xy <- cbind(pc, X)

# now given this alpha, 
# we do the COMBI strategy, but mixing in what i need
# i do the ML-SCR and ML-RCV
# but each time doing CV for cd structure and lasso penalty
# then in the end I compare the two models and determine the end-resulting components

# we fix the ridge at a reasonable value. 

# ML-SCREE #
# we do scree test. 
# we observe where the elbow happens, calculating the loss at each different R value

R_range <- c(1:19)

losses <- c()

y_fit <- c()

for (i in 1:length(R_range)){
  fit1 <- sscovr2(X = Xy[,-1], Y = Xy[,1], blockcols = c(14, 36), 
                  R = R_range[i], lambda1 =  rep(0, R_range[i]), 
                  lambda2 = 0.01, 
                  cdstructure = matrix(rep(0,R_range[i]), nrow = 1), 
                  alpha = a_ml, inits = "rational", nrstart = 1, lambda_y = 1, MAXITER = 10000, stop_value = 1e-7, ridge_in_between = T, ridge_end = T)
  
  losses[i] <- fit1$cd_loss
  
  y_fit[i] <- sum((Xy[,-1] %*% fit1$cd_results[[1]]$W %*% fit1$cd_results[[1]]$P[1,] - Xy[,1])^2)
  
}

#scree plot
plot(losses, xlab = "Principal Component",
     ylab = "loss values",
     xlim = c(NULL), 
     type = "b")

#scree plot
plot(y_fit, xlab = "Principal Component",
     ylab = "loss values",
     xlim = c(NULL), 
     type = "b")

ml_screeplot <- 
  qplot(c(1:19), y_fit) + 
  geom_line() + 
  xlab("Number of principal components") + 
  ylab("in-sample-error") +
  theme_bw() 

# so according to the baseline alpha, we need 3 components
# now based on this 3 component solution
# we cross validate for ridge and alpha

alpha_range <- c(0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)

ridge_range <-  c(5, 3, 1, 0.5, 0.1, 0.05, 0.01, 0.005, 0.001)

ranges <- expand.grid(ridge_range, alpha_range)

colnames(ranges) <- c("ridge", "alpha")

set.seed(2566)
alpha_seed <- sample(nrow(ranges))

alpha_cve <- matrix(NA, ncol = 3+1+1, nrow = nrow(ranges))
alpha_se <- matrix(NA, ncol = 3+1+1, nrow = nrow(ranges))

for (i in 1:nrow(ranges)){
  cv_i <- sscovr_cv_general(X = Xy[,-1], Y = Xy[,1], 
                            ridge = ranges[i,]$ridge, 
                            lasso = rep(0, 3), nrFolds = nrFolds,
                            R = 3, 
                            cdstructure = matrix(rep(0, 3), nrow = 1), 
                            alpha = ranges[i,]$alpha, 
                            inits = "rational", nrstart = 1, 
                            blockcols = c(14, 36), 
                            seed = alpha_seed[i], lambda_y = 1, 
                            stop_value = 1e-7)
  
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

# CV for CD structure #

# repeated double cross validation from filzmoser #
# it's pretty simple. the mean CV and SD CV are all calculated the same way
# except, before the CV process we separate the set already into
# test and calibration set
# CV is then carried out in the calibration set 
# so calibration set divides into training set and validation set
# after CV on the calibration set, we can do 1SE rule
# to find the right model. Find that right model and store it
# do the process many times:
# number of partitions for calibration / test set 
# and also number of repetitions
# in the end, the model with the most frequency wins

# all of the combinations of common and distinctive 

cd_range <- gtools::combinations(n = 3, r = 3, v = 0:2, repeats.allowed = TRUE)

nrow(cd_range)

dim(Xy)

nrFolds <- 10 

set.seed(800)
cd_seed <- sample(1:10000, nrow(cd_range))

cd_cve <- matrix(NA, ncol = 3+3, nrow = nrow(cd_range))
cd_se <- matrix(NA, ncol = 3+3, nrow = nrow(cd_range))

for (i in 1:nrow(cd_range)){
  
  cv_i <- sscovr_cv_general(X = Xy[,-1], Y = Xy[,1], 
                            ridge = alpha_ridge_chosen$ridge, 
                            lasso = rep(0,3), nrFolds = nrFolds,
                            R = 3, cdstructure = cd_range[i,], 
                            alpha = alpha_ridge_chosen$alpha, 
                            inits = "rational", 
                            nrstart = 1, blockcols = c(14, 36), 
                            seed = cd_seed[i], lambda_y = 1, 
                            MAXITER = 10000, stop_value = 1e-7)
  
  cd_cve[i,] <- c(cv_i$cve, cd_range[i,])
  cd_se[i,] <- c(cv_i$se, cd_range[i,] )
  
  print(rep(i, 10))
}
    
cd_cve

# error on y #
cd_cv_df <- data.frame(cd = c(1:nrow(cd_range)), error=cd_cve[,2], 
                       lower = cd_cve[,2] - cd_se[,2], 
                       upper = cd_cve[,2] + cd_se[,2])

CD_y <- 
  ggplot() + 
  geom_errorbar(data=cd_cv_df, mapping=aes(x=cd, ymin=upper, ymax=lower), width=0.1, size=0.5, color="black") + 
  geom_point(data=cd_cv_df, mapping=aes(x=cd, y=error), size=1.5, shape=1, fill="white") +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) + 
  ylab("CV error") +
  xlab("Common and distinctive component structures")

which.min(cd_cve[,2])

cd_scr <- cd_range[which.min(cd_cve[,2]),]

# CV for lasso #
lasso_range <- seq(0.001, 4, length.out = 20)

nrFolds <- 10

lasso_cve <- matrix(NA, nrow = length(lasso_range), ncol = 4)
lasso_se <- matrix(NA, nrow = length(lasso_range), ncol = 4)

set.seed(5656)
lasso_seed <- sample(1:10000, length(lasso_range))
    
for (i in 1:length(lasso_range)){
  cv_i <- tryCatch(sscovr_cv_general(X = Xy[,-1], Y = Xy[,1], 
                                     ridge = alpha_ridge_chosen$ridge, 
                                     lasso = c(rep(lasso_range[i], 3)), 
                                     nrFolds = nrFolds,
                                     R = 3, 
                                     cdstructure = matrix(cd_scr, nrow = 1), 
                                     alpha = alpha_ridge_chosen$alpha, inits = "rational", 
                                     nrstart = 1, blockcols = c(14, 36), 
                                     seed = lasso_seed[i], lambda_y = 1, 
                                     MAXITER = 10000, stop_value = 1e-7), 
                   error = function(e) NA)
  
  if (anyNA(cv_i)){
    lasso_cve[i,] <- c(rep(NA,3),  lasso_range[i])
    lasso_se[i,] <- c(rep(NA,3), lasso_range[i])
  } else {
    lasso_cve[i,] <- c(cv_i$cve, lasso_range[i])
    lasso_se[i,] <- c(cv_i$se, lasso_range[i])
  }

  print(i)  
}
  
lasso_cve
lasso_se

# error from y only #
lasso_cv_df <- data.frame(lasso = lasso_cve[,4], error=lasso_cve[,2], 
                          lower = lasso_cve[,2] - lasso_se[,2], 
                          upper = lasso_cve[,2] + lasso_se[,2])

lasso_cv_df <- lasso_cv_df[!apply(lasso_cv_df,1,anyNA),]

lasso_y <- 
  ggplot() + 
  geom_errorbar(data=lasso_cv_df, mapping=aes(x=lasso, ymin=upper, ymax=lower), width=0.1, size=0.5, color="black") + 
  geom_point(data=lasso_cv_df, mapping=aes(x=lasso, y=error), size = 1.5, shape=1, fill="white") +
  theme_bw() +
  scale_x_continuous(name ="Lasso penalty", breaks = c(0:10)) +
  ylab("CV error") 

lasso_cv_df2 <- lasso_cv_df[lasso_cv_df[,2] < lasso_cv_df[which.min(lasso_cv_df[,2]),4],]

lasso_cv_df2[nrow(lasso_cv_df2),]

lasso_scr <- lasso_cv_df2[nrow(lasso_cv_df2),1]

model1 <- sscovr2(X = Xy[,-1], Y = Xy[,1], 
                  blockcols = c(14, 36), R = 3, 
                  lambda1 =  rep(lasso_scr, 3), 
                  lambda2 = alpha_ridge_chosen$ridge, 
                  cdstructure = matrix(cd_scr, nrow = 1), 
                  alpha = alpha_ridge_chosen$alpha, inits = "rational", 
                  nrstart = 1, lambda_y = 1, 
                  MAXITER = 10000, stop_value = 1e-7, 
                  ridge_in_between = T, ridge_end = T)

model1$cd_results[[1]]$W
model1$cd_results[[1]]$Py

fitted <- Xy[,-1] %*% model1$cd_results[[1]]$W %*% model1$cd_results[[1]]$P[1,] # fitted y value

(sscovr_rsq_y <- 1 - sum((Xy[,1] - fitted)^2) / sum(Xy[,1]^2))

# final cross validation error #
set.seed(77)
final_seed <- sample(1:1000, 1)

cv_i <- sscovr_cv_general(X = Xy[,-1], Y = Xy[,1], 
                          ridge = alpha_ridge_chosen$ridge, 
                          lasso = rep(lasso_scr, 3), nrFolds = 10,
                          R = 3, cdstructure = matrix(cd_scr, nrow = 1), 
                          alpha = alpha_ridge_chosen$alpha, 
                          inits = "rational", 
                          nrstart = 1, blockcols = c(14, 36), 
                          seed = final_seed, lambda_y = 1, 
                          MAXITER = 10000, stop_value = 1e-7)
  
(ssc_cv_final <- cv_i)

# PLS (sgcca) #
# number of zeros in the W:

A <- list(Xy[,-1], Xy[,1])

C <- matrix(c(0, 1, 1, 0), 2, 2)

# cross validation for l1 penalty #
nrFolds = 10

pls_l1_range <- seq(1/sqrt(50), 1, length.out = 20)

set.seed(25)
pls_l1_seed <- sample(1:1000, length(pls_l1_range))


pls_l1_cve <- matrix(NA, ncol = 2, nrow = length(pls_l1_range))
pls_l1_se <- matrix(NA, ncol = 2, nrow = length(pls_l1_range))

for (i in 1:length(pls_l1_range)){
  cv_i <- sgcca_cv_general(X = Xy[,-1], Y = Xy[,1], C = C, c1 = c(pls_l1_range[i], 1), 
                           ncomp = c(3,1), 
                           scheme = "factorial", verbose = TRUE, 
                           nrFolds = 10, seed = pls_l1_seed[i])
  
  pls_l1_cve[i,] <- c(cv_i$cve, pls_l1_range[i])
  pls_l1_se[i,] <- c(cv_i$se, pls_l1_range[i])
  
  print(rep(i, 10))
}


pls_l1_cve
pls_l1_se

# error from y only #
pls_l1_cv_df <- data.frame(lasso = pls_l1_cve[,2], error=pls_l1_cve[,1], 
                          lower = pls_l1_cve[,1] - pls_l1_se[,1], 
                          upper = pls_l1_cve[,1] + pls_l1_se[,1])

pls_l1 <- 
  ggplot() + 
  geom_errorbar(data=pls_l1_cv_df, mapping=aes(x=lasso, ymin=upper, ymax=lower), width=0.03, size=0.5, color="black") + 
  geom_point(data=pls_l1_cv_df, mapping=aes(x=lasso, y=error), size = 1.5, shape=1, fill="white") +
  theme_bw() +
  scale_x_continuous(name ="Lasso penalty", breaks = c(0:10)) +
  ylab("CV error") 

pls_l1_cv_df2 <- pls_l1_cv_df[pls_l1_cv_df[,2] < pls_l1_cv_df[which.min(pls_l1_cv_df[,2]),4],]

pls_l1_cv_df2[1,]

pls_l1_value <- pls_l1_cv_df2[1,]$lasso

# PLS: final model #

# setting the c1 parameter so that the number of sparse weights match the sscovr
sgcca1 <- RGCCA::sgcca(A,C, c1=c(pls_l1_value,1),ncomp=c(3,1), scheme="factorial",verbose=TRUE)
#note that RGCCA scales all variables to unit var

s1<-sgcca1$AVE[[1]]
sum(s1[[1]])#VAF by the two components in pred. data reported in Table 1
sum(s1[[2]])#VAF by the two components in outcome reported in Table 1

sgcca_w <- sgcca1$astar[[1]]

apply(sgcca_w == 0,2,sum)
sum(sgcca_w == 0)

sgccatmat <- Xy[,-1] %*% sgcca_w

sgccareg<-lm(Xy[,1]~sgccatmat)#find optimal regression weights to use for out-of-sample prediction

pls_fitted <- Xy[,-1] %*% sgcca_w %*% sgccareg$coefficients[-1]

(pls_rsq <- 1 - sum((Xy[,1] - pls_fitted)^2) / sum(Xy[,1]^2))
# R squared of pls

# cross validation for PLS #
set.seed(12)
sgcca_seed <- sample(1:1000, 1)

cv_i <- sgcca_cv_general(X = Xy[,-1], Y = Xy[,1], C = C, c1 = c(pls_l1_value, 1), 
                         ncomp = c(3,1), 
                         scheme = "factorial", verbose = TRUE, 
                         nrFolds = 10, seed = sgcca_seed)

(pls_cv_final <- cv_i)

# model selection for PCR

# cross validation for ridge
ridge_range <- c(seq(from = 5, to = 0.001, length.out = 20))

nrFolds <- 10

pcr_ridge_cve <- matrix(NA, nrow = length(ridge_range), ncol = 2)
pcr_ridge_se <- matrix(NA, nrow = length(ridge_range), ncol = 2)

set.seed(119)
ridge_seed <- sample(1:10000, length(ridge_range))

for (i in 1:length(ridge_range)){
  
  cv_i <- scad_cv(X = Xy[,-1], R = 3, ridge = ridge_range[i], 
                  lasso = rep(0, 3), nrFolds = 10, 
                  fixW = fixW, MAXITER = 10000, stop_value = 1e-7, 
                  seed = ridge_seed[i])
  
  pcr_ridge_cve[i,] <- c(cv_i$cve, ridge_range[i])
  pcr_ridge_se[i,] <- c(cv_i$se, ridge_range[i])
  
  print(i)  
}

pcr_ridge_cve
pcr_ridge_se

# error from y only #
pcr_ridge_cv_df <- data.frame(ridge = pcr_ridge_cve[,2], error=pcr_ridge_cve[,1], 
                              lower = pcr_ridge_cve[,1] - pcr_ridge_se[,1], 
                              upper = pcr_ridge_cve[,1] + pcr_ridge_se[,1])

pcr_ridge <- 
  ggplot() + 
  geom_errorbar(data=pcr_ridge_cv_df, mapping=aes(x=ridge, ymin=upper, ymax=lower), width=0.005, size=0.2, color="black") + 
  geom_point(data=pcr_ridge_cv_df, mapping=aes(x=ridge, y=error), size = 1.5, shape=1, fill="white") +
  theme_bw() +
  scale_x_continuous(name ="ridge penalty", breaks = round(ridge_range, 2)) +
  ylab("CV error") 

pcr_ridge_cv_df[pcr_ridge_cv_df[,2] < pcr_ridge_cv_df[which.min(pcr_ridge_cv_df[,2]),4],]

pcr_ridge_scr <- ridge_range[20]

# cross validation for common-distinctive #

cd_range <- gtools::combinations(n = 3, r = 3, v = 0:2, repeats.allowed = TRUE)

set.seed(999)
cd_seed <- sample(1:1000, nrow(cd_range))
  
pcr_cd_cve <- matrix(NA, ncol = 1+3, nrow = nrow(cd_range))
pcr_cd_se <- matrix(NA, ncol = 1+3, nrow = nrow(cd_range))
  

for (i in 1:nrow(cd_range)){
  fixW <- matrix(0, nrow = 50, ncol = 3)
  
  for (r in 1:3){
    if(cd_range[i, r] == 1){
      fixW[1:14, r] <- 1
    }
    
    if (cd_range[i, r] == 2){
      fixW[15:50, r] <- 1
    }
  }
    
   cv_i <- scad_cv(X = Xy[,-1], R = 3, ridge = pcr_ridge_scr, 
                   lasso = rep(0, 3), nrFolds = 10, 
                   fixW = fixW, MAXITER = 10000, stop_value = 1e-7, 
                   seed = cd_seed[i])
    
    pcr_cd_cve[i,] <- c(cv_i$cve, cd_range[i,])
    pcr_cd_se[i,] <- c(cv_i$se, cd_range[i,] )
    
    print(rep(i, 10))
  }

pcr_cd_cve
pcr_cd_se

# error on x #
pcr_cd_cv_df <- data.frame(cd = c(1:nrow(cd_range)), error=pcr_cd_cve[,1], 
                       lower = pcr_cd_cve[,1] - pcr_cd_se[,1], 
                       upper = pcr_cd_cve[,1] + pcr_cd_se[,1])
pcr_CD <- 
  ggplot() + 
  geom_errorbar(data=pcr_cd_cv_df, mapping=aes(x=cd, ymin=upper, ymax=lower), width=0.1, size=0.5, color="black") + 
  geom_point(data=pcr_cd_cv_df, mapping=aes(x=cd, y=error), size=1.5, shape=1, fill="white") +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) + 
  ylab("CV error") +
  xlab("Common and distinctive component structures")

which.min(pcr_cd_cve[,1])

pcr_cd_scr <- cd_range[9,]

fixW <- matrix(0, nrow = 50, ncol = 3)

for (r in 1:3){
  if(cd_range[9, r] == 1){
    fixW[1:14, r] <- 1
  }
  
  if (cd_range[9, r] == 2){
    fixW[15:50, r] <- 1
  }
}


# CV for lasso #
lasso_range <- c(seq(from = 0.000001, to = 0.2, length.out = 20))

nrFolds <- 10

pcr_lasso_cve <- matrix(NA, nrow = length(lasso_range), ncol = 2)
pcr_lasso_se <- matrix(NA, nrow = length(lasso_range), ncol = 2)

set.seed(119)
lasso_seed <- sample(1:10000, length(lasso_range))

for (i in 1:length(lasso_range)){
  
  cv_i <- scad_cv(X = Xy[,-1], R = 3, ridge = pcr_ridge_scr, 
                  lasso = rep(lasso_range[i], 3), nrFolds = 10, 
                  fixW = fixW, MAXITER = 10000, stop_value = 1e-7, 
                  seed = lasso_seed[i])
  
  pcr_lasso_cve[i,] <- c(cv_i$cve, lasso_range[i])
  pcr_lasso_se[i,] <- c(cv_i$se, lasso_range[i])
  
  print(i)  
}

pcr_lasso_cve
pcr_lasso_se

# error from y only #
pcr_lasso_cv_df <- data.frame(lasso = pcr_lasso_cve[,2], error=pcr_lasso_cve[,1], 
                          lower = pcr_lasso_cve[,1] - pcr_lasso_se[,1], 
                          upper = pcr_lasso_cve[,1] + pcr_lasso_se[,1])

pcr_lasso <- 
  ggplot() + 
  geom_errorbar(data=pcr_lasso_cv_df, mapping=aes(x=lasso, ymin=upper, ymax=lower), width=0.005, size=0.2, color="black") + 
  geom_point(data=pcr_lasso_cv_df, mapping=aes(x=lasso, y=error), size = 1.5, shape=1, fill="white") +
  theme_bw() +
  scale_x_continuous(name ="Lasso penalty", breaks = round(seq(min(lasso_range), 0.2, length.out = 5),2)) +
  ylab("CV error") 

pcr_lasso_cv_df[pcr_lasso_cv_df[,2] < pcr_lasso_cv_df[which.min(pcr_lasso_cv_df[,2]),4],]

pcr_lasso_scr <- lasso_range[2]

pcr_model <- sparseSCAcpp(X = Xy[,-1], Q = 3, 
             RIDGE = 0.01, LASSO = rep(pcr_lasso_scr,3), 
             fixW = fixW, maxItrOuterloop = 100000, nStarts = 1,
             print = FALSE, tol = 10^-8)

pcr_model$W

scadtmat <- Xy[,-1] %*% pcr_model$W
scadreg <- lm(Xy[,1] ~ scadtmat)
scadreg$coefficients[is.na(scadreg$coefficients)] <- 0

pcr_fitted <- Xy[,-1] %*% pcr_model$W %*% scadreg$coefficients[-1]

(pcr_rsq <- 1 - sum((Xy[,1] - pcr_fitted)^2) / sum(Xy[,1]^2))

# pcr final cross validation #

set.seed(414)

pcr_seed <- sample(1:1000, 1)

cv_i <- scad_cv_y(X = Xy[,-1], y = Xy[,1], R = 3, ridge = pcr_ridge_scr, 
                lasso = rep(pcr_lasso_scr, 3), nrFolds = 10, 
                fixW = fixW, MAXITER = 10000, stop_value = 1e-7, 
                seed = pcr_seed, lambda_y = 1)

(pcr_cv_final <- cv_i)


# SPCOVR: CV for lasso #
lasso_range <- seq(0.01, 8, length.out = 20)

nrFolds <- 10

spc_lasso_cve <- matrix(NA, nrow = length(lasso_range), ncol = 4)
spc_lasso_se <- matrix(NA, nrow = length(lasso_range), ncol = 4)

set.seed(100)
spc_lasso_seed <- sample(1:10000, length(lasso_range))

for (i in 1:length(lasso_range)){
  cv_i <- tryCatch(sscovr_cv_general(X = Xy[,-1], Y = Xy[,1], 
                                     ridge = alpha_ridge_chosen$ridge, 
                                     lasso = c(rep(lasso_range[i], 3)), 
                                     nrFolds = nrFolds,
                                     R = 3, 
                                     cdstructure = matrix(rep(0,3), nrow = 1), 
                                     alpha = alpha_ridge_chosen$alpha, inits = "rational", 
                                     nrstart = 1, blockcols = c(14, 36), 
                                     seed = spc_lasso_seed[i], lambda_y = 1, 
                                     MAXITER = 10000, stop_value = 1e-7), 
                   error = function(e) NA)
  
  if (anyNA(cv_i)){
    spc_lasso_cve[i,] <- c(rep(NA,3),  lasso_range[i])
    spc_lasso_se[i,] <- c(rep(NA,3), lasso_range[i])
  } else {
    spc_lasso_cve[i,] <- c(cv_i$cve, lasso_range[i])
    spc_lasso_se[i,] <- c(cv_i$se, lasso_range[i])
  }
  
  print(i)  
}

spc_lasso_cve
spc_lasso_se

# error from y only #
spc_lasso_cv_df <- data.frame(lasso = spc_lasso_cve[,4], error=spc_lasso_cve[,2], 
                          lower = spc_lasso_cve[,2] - spc_lasso_se[,2], 
                          upper = spc_lasso_cve[,2] + spc_lasso_se[,2])

spc_lasso_cv_df <- spc_lasso_cv_df[!apply(spc_lasso_cv_df,1,anyNA),]

spc_lasso_y <- 
  ggplot() + 
  geom_errorbar(data=spc_lasso_cv_df, mapping=aes(x=lasso, ymin=upper, ymax=lower), width=0.1, size=0.5, color="black") + 
  geom_point(data=spc_lasso_cv_df, mapping=aes(x=lasso, y=error), size = 1.5, shape=1, fill="white") +
  theme_bw() +
  scale_x_continuous(name ="Lasso penalty", breaks = c(0:10)) +
  ylab("CV error") 

spc_lasso_cv_df[spc_lasso_cv_df[,2] < spc_lasso_cv_df[which.min(spc_lasso_cv_df[,2]),4],]

spc_lasso_scr <- lasso_range[14]

spc_model <- sscovr2(X = Xy[,-1], Y = Xy[,1], 
                  blockcols = c(14, 36), R = 3, 
                  lambda1 =  rep(spc_lasso_scr, 3), 
                  lambda2 = alpha_ridge_chosen$ridge, 
                  cdstructure = matrix(rep(0,3), nrow = 1), 
                  alpha = alpha_ridge_chosen$alpha, inits = "rational", 
                  nrstart = 1, lambda_y = 1, 
                  MAXITER = 10000, stop_value = 1e-7, 
                  ridge_in_between = T, ridge_end = T)

spc_model$cd_results[[1]]$W
spc_model$cd_results[[1]]$Py

spc_fitted <- Xy[,-1] %*% spc_model$cd_results[[1]]$W %*% spc_model$cd_results[[1]]$P[1,] # fitted y value

(spcovr_rsq_y <- 1 - sum((Xy[,1] - spc_fitted)^2) / sum(Xy[,1]^2))

# final cross validation error #
set.seed(456)
final_seed <- sample(1:1000, 1)

cv_i <- sscovr_cv_general(X = Xy[,-1], Y = Xy[,1], 
                          ridge = alpha_ridge_chosen$ridge, 
                          lasso = rep(spc_lasso_scr, 3), nrFolds = 10,
                          R = 3, cdstructure = matrix(rep(0,3), nrow = 1), 
                          alpha = alpha_ridge_chosen$alpha, 
                          inits = "rational", 
                          nrstart = 1, blockcols = c(14, 36), 
                          seed = final_seed, lambda_y = 1, 
                          MAXITER = 10000, stop_value = 1e-7)

(spc_cv_final <- cv_i)$cve




# final cross validation error plot #
final_cv_df <- data.frame(method = factor(c("SSCovR", "SPCovR", "PCR-SCaDS", "SGCCA"), levels = c("SSCovR", "SPCovR", "PCR-SCaDS",  "SGCCA")), 
                          error = c(ssc_cv_final$cve[2], spc_cv_final$cve[2],
                                    pcr_cv_final$cve[2], pls_cv_final$cve),
                          lower = c(ssc_cv_final$cve[2] - ssc_cv_final$se[2], 
                                    spc_cv_final$cve[2] - spc_cv_final$se[2], 
                                    pcr_cv_final$cve[2] - pcr_cv_final$se[2],
                                    pls_cv_final$cve - pls_cv_final$se), 
                          upper = c(ssc_cv_final$cve[2] + ssc_cv_final$se[2], 
                                    spc_cv_final$cve[2] + spc_cv_final$se[2], 
                                    pcr_cv_final$cve[2] + pcr_cv_final$se[2],
                                    pls_cv_final$cve + pls_cv_final$se))




final_cv_plot <- 
  ggplot() + 
  geom_errorbar(data=final_cv_df, mapping=aes(x=method, ymin=upper, ymax=lower), width=0.2, size=1, color="black") + 
  geom_point(data=final_cv_df, mapping=aes(x=method, y=error), size=2, colour = "red") +
  theme_bw() +
  ylab("CV error") +
  xlab("Method")

# save(list = c("y_fit",
#               "losses", "ml_screeplot", "CD_y", "lasso_y", "ssc_cv_final", "sscovr_rsq_y",
#               "lasso_scr", "alpha_ridge_chosen", "cd_scr",
#               
#               "pls_rsq", "pls_cv_final", "pls_l1", "pls_l1_value",
#               
#               "pcr_ridge", "pcr_CD","pcr_lasso", "pcr_rsq", "pcr_cv_final",
#               "pcr_ridge_scr", "pcr_lasso_scr", "pcr_cd_scr",
#               
#               "spc_lasso_y", "spc_rsq", "spc_cv_final",
#               "spc_lasso_scr", 
#               
#               "final_cv_plot", "final_cv_df"),
#      file = "C:\\Users\\park\\Desktop\\Rstudio Github - Copy\\Project2_SPCovR\\paper2_writing\\empirical_results.Rdata")
