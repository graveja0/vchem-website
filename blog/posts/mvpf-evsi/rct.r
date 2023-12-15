# library(tidyverse)
# library(MASS)
# library(glue)
# library(janitor)
# 
# 
# C = 1000               # Total number of clusters in population
# n_ = NULL              # Sample size (can be set to NULL if  cluster-level randomization)
# c_ = 100               # Cluster size (can be set to NULL if unit-level randomization)
# p_treat = 0.5          # Probability of treatment
# tau = 65               # Treatment effect

sim_rct <- function(C = 1000, n_ = 100, c_ = NULL, p_treat = 0.5, tau = 65) {

    ###################################################
    # 1. Sample Population and Assign Treatment Status
    ###################################################

    cluster_ <- 
        tibble(cluster_id = 1:C, 
               gamma_j = rnorm(C, mean = 15, sd = 3),
               cluster_size = rpois(C, lambda = 4)) 
    
    N <- sum(cluster_$cluster_size)  # Population sample size 

    if (!is.null(c_)) {
        sampled <- sample(1:C,c_, replace = FALSE)
        cluster_ <- 
            cluster_[sampled, ]
    
        num_treated <-  c_ * p_treat
        W_ <- runif(c_)
        cluster_$W <- as.integer(order(W_) <= num_treated)
        
        df_ <- 
            cluster_ %>% 
            uncount(cluster_size) 
        
    } else if (!is.null(n_)) {
        individ_ <- 
            cluster_ %>% 
            uncount(cluster_size)
        
        sampled <- sample(1:nrow(individ_), n_, replace=FALSE)
        
        individ_ <- 
            individ_[sampled, ]
        
        num_treated <- n_ * p_treat
        W_ <- runif(n_)
        individ_$W <- as.integer(order(W_) <= num_treated)
        
        df_ <- 
            individ_ 
    }

    # Sample size
    n <- nrow(df_)  

    ###################################################
    # 2. Copula-based sampling of potential outcomes
    ###################################################
     
    # Correlation between potential outcomes is encoded in sigma   
    sigma <- matrix(c(1,-.6,-.6,1),byrow=TRUE,nrow=2,ncol=2)
    
    # Sample from multivariate normal with Sigma = sigma
    z <- mvrnorm(n, mu = rep(0,2), Sigma = sigma, empirical = TRUE)
    u <- pnorm(z)
    
    Y_1 <- qlnorm(u[,1], meanlog = log(1000) - 0.5 * 0.2, 0.2 ) + df_$gamma_j + rgamma(n, shape = tau / 2, scale = 2)
    Y_0 <- qlnorm(u[,2], meanlog = log(1000) - 0.5 * 0.2, 0.2 ) + df_$gamma_j 
    df_$Y <- ifelse(df_$W==1, Y_1, Y_0)
    
    df <- df_ %>% 
        mutate(id = 1:nrow(df_)) %>% 
        dplyr::select(id, cluster_id, W, Y)
     
    return(df_)
    
}

