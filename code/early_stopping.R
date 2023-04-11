# compare stopping criterion

source("helper_2.R")

log_data <- list()
n_list <- c(200, 400, 800, 1600)
sims <- 100
coverage <- data.frame(n = NA,
                       eps_term = NA,
                       plgb_term = NA,
                       final_term = NA)
avg_iters <- data.frame(n = NA,
                    eps_term = NA,
                    plgb_term = NA,
                    final_term = NA)

x <- 1:16
#ate <- mean(plogis(0.25*1+0.25*(x-8.5)+0.1*x) - plogis(0.25*(x-8.5)))
ate <- mean(plogis(0.25+0.25*(x-8.5)+0.5*(sin(x)))-plogis(0.25*(x-8.5)))


# simple data generation process
for (n in n_list){
  n_log_data <- list()
  covered <- c(n,rep(0, 3))
  iters <- c(n,rep(0, 3))
  for (i in 1:sims){
    data <- generate_data_simple(n)
    Y <- data$Y
    X <- data$X
    A <- data$A
    results <- tmle_v0(Y,A,X)
    
    for (j in 2:4){
      result <- results[[j]]
      covered[j] <- covered[j]+((result$CI[1] <= ate) && (result$CI[2] >= ate))
      iters[j] <- iters[j]+result$iters 
    }
    n_log_data <- append(n_log_data, results[[1]])
    
  }
  covered[2:4] <- covered[2:4]/sims
  iters[2:4] <- iters[2:4]/sims
  
  coverage <- rbind(coverage, covered)
  avg_iters <- rbind(avg_iters, iters)
  log_data <- append(log_data, list(n_log_data))
}


# get tables
print(xtable(coverage[-1,]), include.rownames=F)
print(xtable(avg_iters[-1,]), include.rownames=F)



