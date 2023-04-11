# Functions for TMLE simulations

# TMLE vs OTMLE vs ROTMLE
library(tmle)
library(SuperLearner)
library(CVXR)
library(dbarts)
library(Rforestry)
library(LaplacesDemon)

library(tibble)
# generate complex data
generate_data <- function(n){ 
  W1 <- (rbinom(n, size=1, prob=0.2)) # binary confounder
  W2 <- (rbinom(n, size=1, prob=0.5)) # binary confounder
  W3 <- round(runif(n, min=2, max=7)) # continuous confounder
  W4 <- round(runif(n, min=0, max=4)) # continuous confounder
  A  <- rbinom(n, size=1, prob= plogis(-0.2 + 0.2*W2 + log(0.1*W3) + 0.3*W4 + 0.2*W1*W4)) # binary treatment depends on confounders
  Y <- rbinom(n, size=1, prob= plogis(-1 + A - 0.1*W1 + 0.2*W2 + 0.3*W3 - 0.1*W4 + sin(0.1*W2*W4))) # binary outcome depends on confounders
  return(tibble(Y, W1, W2, W3, W4, A))
}

# generate simple data
generate_data_simple <- function(n){
  X <- (sample(1:16, size = n, replace = T))
  A <- rbinom(n, size=1, prob = 0.1 + 1/20*X)
  Y <- rbinom(n, size=1, prob = plogis(0.25*A+0.25*(X-8.5)+0.5*(sin(X)*A)))
  X <- as.factor(X)
  return(tibble(Y,A,X))
}

generate_data_prev <- function(n,S){
  X <- sample(1:S, n, replace = T)
  treatment_probs <- 0.5*(X%%4==1) + 0.5*(X%%4==2) + 0.5*(X%%4==3) + 0.05*(X%%4==0)
  A <- rbinom(n, size=1, prob = treatment_probs)
  success_prob <- 0.4*(X%%4==1) + 0.5*(X%%4==2 & A==0) + 
    0.2*(X%%4==3 & A==0) + 0.1*(X%%4==0 & A==0) + 0.2*(X%%4==2 & A==1) + 
    0.5*(X%%4==3 & A==1) + 0.9*(X%%4==0 & A==1)
  #print(success_prob)
  Y <- rbinom(n, size=1, prob = success_prob)
  return(data.frame(Y=Y, A=A, X=as.factor(X)))
}


# influence function direction (only changing Q)
# Q_est, g_est should be the estimated Q and g values from our function
phi.Y <- function(Y,A,X,Q_est,g_est){
  weights <- A/(g_est) + (1-A)/(1-g_est)
  return(weights*(Y-Q_est))
}

# plug_in-bias 
plug_in_bias <- function(Y,A,X,Q_est,g_est){
  return(mean(phi.Y(Y,A,X,Q_est,g_est)))
}

# estimate variance
variance_est <- function(Y,A,X,
                         Q_est, Q1, Q0, g_est, n){
  
  p1 <- mean(phi.Y(Y,A,X,Q_est,g_est)^2)
  p2 <- mean( (Q1-mean(Q1))^2 )
  p3 <- mean( (Q0-mean(Q0))^2 )
  
  return((p1+p2+p3)/n)
}

# update Q
update_Q <- function(Y,A,X,p.y,Q_est){
  probs <- data.frame(unique(cbind(X,A,Y,p.y)[Y==1, ]))
  for (i in 1:nrow(probs)){
    x <- probs$X[i]
    a <- probs$A[i]
    Q_est[X==x & A==a] <- probs$p.y[i]
  }
  return(Q_est)
}

getQs <- function(Y,A,X,Q_est,n, Q1_est, Q0_est){
  Q1 <- Q1_est
  Q0 <- Q0_est
  probs <- data.frame(unique(cbind(X,A,Y,Q_est)[Y==1,]))
  for (x in unique(X)){
    q1x <- probs$Q_est[probs$X==x & probs$A==1]
    q0x <- probs$Q_est[probs$X==x & probs$A==0]
    if (length(q1x) != 0){
      Q1[X==x] <- q1x
    }
    if (length(q0x) != 0){
      Q0[X==x] <- q0x
    }
  }
  return(data.frame(Q1,Q0))
}

# get ate, variance, confidence intervals (95%)
get_values <- function(Y,A,X,Q_est,g_est,n,k,Q1,Q0){
  # get Q1, Q0
  Qs <- getQs(Y,A,X,Q_est,n,Q1,Q0)
  ATE <- mean(Qs$Q1-Qs$Q0)
  var <- variance_est(Y,A,X,Q_est, Qs$Q1, Qs$Q0, g_est, n)
  
  # confidence intervals:
  se <- sqrt(var)
  lower <- ATE - qnorm(0.975)*se
  upper <- ATE + qnorm(0.975)*se
  
  return(list(ATE = ATE,
              Variance = var,
              CI = c(lower, upper),
              iters = k))
}



# TMLE Base Form: ONLY FOR ATE Estimation, need to specify a super-learned q and g
# first.
#' @Y = outcome variable (binary)
#' @A = treatment indicator (binary)
#' @W = control variables
tmle_v0 <- function(Y, A, X){
  
  # To initialize, start with the super-learned distributions
  
  SL.library <- c("SL.glm", "SL.randomForest", "SL.mean")
  
  Q <- SuperLearner(Y, data.frame(A,X), 
                    SL.library = SL.library, 
                    family = binomial())
  g <- SuperLearner(A, data.frame(X), 
                    SL.library = SL.library,
                    family = binomial())
  
  n <- length(Y)
  
  Q_est <- as.vector(predict(Q, data.frame(A,X), onlySL = T)$pred)
  Q1 <- as.vector(predict(Q, data.frame(A=1,X), onlySL = T)$pred)
  Q0 <- as.vector(predict(Q, data.frame(A=0,X), onlySL = T)$pred)
  g_est <- as.vector(predict(g, data.frame(X), onlySL = T)$pred)
  
  # initialize eps and plgb
  plgb <- plug_in_bias(Y,A,X,Q_est,g_est)
  eps <- .Machine$double.xmax
  term_cond <- 1/(sqrt(n)*log(n))
  k <- 0
  
  # Store parameter values
  log_data <- data.frame(eps = NA,
                         plug_in_bias = NA)
  
  eps_term <- NA
  plgb_term <- NA
  eps_first <- 0
  plgb_first <- 0
  
  # Define the fluctuation parameters:
  while(eps^2 > 10^(-20) ){
    
    # increment counter
    k <- k+1
    
    # Solve log-likelihood problem
    eps_dummy <- Variable(1)
    h <- phi.Y(Y,A,X,Q_est,g_est)
    p.y <- Q_est*(Y==1) + (1-Q_est)*(Y==0) 
    obj <-  mean(log( (1+eps_dummy*h)*p.y)) 
    
    
    # log likelihood loss
    prob <- Problem(Maximize(obj))
    result <- solve(prob)
    
    # update Q vectors 
    eps <- result$getValue(eps_dummy)
    p.y <- (1+eps*h)*p.y
    Q_est <- update_Q(Y, A, X, p.y , Q_est)
    
    # check to see if termination conditions have been met:
    plgb <- plug_in_bias(Y,A,X,Q_est,g_est)
    if (abs(plgb) <= term_cond && (plgb_first == 0) ){
      plgb_term <- get_values(Y,A,X,Q_est,g_est,n, k, Q1,Q0)
      plgb_first <- 1
    }
    
    if (abs(eps) <= term_cond && (eps_first == 0) ){
      eps_term <- get_values(Y,A,X,Q_est,g_est,n, k,Q1,Q0)
      eps_first <- 1
    }
    
    # store relevant values:
    log_data <- rbind(log_data,
                      c(eps, plgb))
  }
  
  log_data <- log_data[-1,]
  
  # estimate ATE using plug-in distribution
  final_term <- get_values(Y,A,X,Q_est,g_est,n, k, Q1,Q0)
  
  return(list(log_data = log_data,
              plug_in_term = plgb_term,
              eps_term = eps_term,
              final_termination = final_term))
}

