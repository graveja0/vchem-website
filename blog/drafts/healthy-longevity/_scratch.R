# Healthy Longevity 

library(tidyverse)
library(MASS)
library(expm)
library(knitr)
library(kableExtra)
options(scipen = 5) 
transpose <- purrr::transpose
select <- dplyr::select

vec <-  # a simple function to return the vec of an array
    function(x) {
        y <- c(x)
        return(y)
    }

vecperm <- 
    # vecperm
    # function to calculate the vec permutation matrix K of index m,n
    # let X be a m x n matrix, and X' the transpose of X
    # then K satisfies 
    # vec(X') = K*vec(X)
    function(m, n) {
        K <- matrix(0, m * n, m * n)
        a <- matrix(0, m, n)
        
        for (i in 1:m) {
            for (j in 1:n) {
                e <- a
                e[i, j] <- 1
                K <- K + kronecker(e, t(e))
            }
        }
        
        return(K)
    }


gen_wcc <- function (n_cycles, method = c("Simpson1/3", "half-cycle", "none")) 
{
    if (n_cycles <= 0) {
        stop("Number of cycles should be positive")
    }
    method <- match.arg(method)
    n_cycles <- as.integer(n_cycles)
    if (method == "Simpson1/3") {
        v_cycles <- seq(1, n_cycles + 1)
        v_wcc <- ((v_cycles%%2) == 0) * (2/3) + ((v_cycles%%2) != 
                                                     0) * (4/3)
        v_wcc[1] <- v_wcc[n_cycles + 1] <- 1/3
    }
    if (method == "half-cycle") {
        v_wcc <- rep(1, n_cycles + 1)
        v_wcc[1] <- v_wcc[n_cycles + 1] <- 0.5
    }
    if (method == "none") {
        v_wcc <- rep(1, n_cycles + 1)
    }
    return(v_wcc)
}

# Parameterize
params <- list(
    tx_names = c("natural_history"), 
    n_tx = 1,
    
    tr_names = c("Healthy","Sick","Sicker"),
    ab_names = c("DeadBG","DeadDisease"),
    
    omega = 250,
    r_ann = 0.03,
    Delta_t = 1,
    age0 = 20,
    s0T = c(1,0,0,0,0),
    
    rHS1 = 0.05,
    rD = 0.01,
    rS1S2 = 0.12,
    rS1D = 0.02,
    rS2D = 0.03,
    
    uS1 = 0.95,
    uS2 = 0.8,
    uD = 0,
    uH = 1,
    
    dwS1 = 0.1,
    dwS2 = 0.2,
    
    ExR = 
        # Global Burden of Disease Collaborative Network. Global Burden of Disease Study 
        # 2019 (GBD 2019) Reference Life Table. Seattle, United States of America: 
        # Institute for Health Metrics and Evaluation (IHME), 2020
        tibble::tribble(
            ~Age, ~Life.Expectancy,
            0L,       88.8718951,
            1L,      88.00051053,
            5L,      84.03008056,
            10L,      79.04633476,
            15L,       74.0665492,
            20L,      69.10756792,
            25L,      64.14930031,
            30L,       59.1962771,
            35L,      54.25261364,
            40L,      49.31739311,
            45L,      44.43332057,
            50L,      39.63473787,
            55L,      34.91488095,
            60L,      30.25343822,
            65L,      25.68089534,
            70L,      21.28820012,
            75L,      17.10351469,
            80L,      13.23872477,
            85L,      9.990181244,
            90L,      7.617724915,
            95L,      5.922359078
        ) 
)

params <- with(params,modifyList(params,list(
    alpha = length(ab_names),
    tau = length(tr_names), 
    s = length(tr_names)*omega + length(ab_names), #total number of states;s=τω+α
    r_Delta_t = r_ann * Delta_t ,
    cycles = omega/Delta_t,
    ages = (0:(omega/Delta_t))*Delta_t + age0,
    fExR <- function(x) pmax(0,unname(Hmisc::approxExtrap(ExR$Age, ExR$Life.Expectancy,xout = x)$y)),
    # Age transition Matrix
    D = {
        # Create diagonal age advancement matrix
        D <- matrix(0, omega, omega)
        vec <- rep(1, omega)
        D[row(D) == col(D)] <- vec
        D
    }
)))

# note that because we have a time homogeneous model here, we do not need to do this separately
# by cycle. but in principle, we could create a list object with each of these matrices defined
# for each time cycle. 

mR_nh <- 
    with(params,{
        matrix(c(
        -(rHS1 + rD), rHS1, 0, rD, 0 ,
        0, -(rS1S2 + rS1D + rD), rS1S2, (rS1D + rD), rS1D,
        0, 0, -(rS2D + rD), (rS2D + rD), rS2D,
        0,0,0,0,0,
        0,0,0,0,0),
        nrow = 5, 
        ncol = 5,
        byrow=TRUE, 
        dimnames = list(c(tr_names,ab_names),
                        c(tr_names,ab_names)
        ))
    })

mR <- with(params,{
    array(as.vector(mR_nh), dim = c(length(tr_names)+ length(ab_names),length(tr_names)+ length(ab_names),length(tx_names)),
          dimnames = list(c(tr_names,ab_names),c(tr_names,ab_names),tx_names))
    }) %>% 
    apply(.,3,function(x) x, simplify=FALSE)

mV <- 
    with(params,{
        mR %>% map(~({
            m <- .x[tr_names,tr_names] 
            diag(m) <- 0
            t(m)
        }))
    })

# mM <- 
#     with(params,{
#         mR %>% map(~({
#             m <- .x[tr_names,ab_names] 
#             diag(m) <- 0
#             t(m)
#         }))
#     })

mQ <- with(params,{
    mR %>% map(~({
        V = t(.x[tr_names,tr_names])
        M = t(.x[tr_names,ab_names])
        zero_ <- matrix(0, nrow = length(tr_names)+length(ab_names), ncol = length(ab_names))
        tmp <- cbind(rbind(V,M),zero_)
        dimnames(tmp) <- list(c(tr_names,ab_names),c(tr_names,ab_names))
        tmp
    }))
})

mP <- with(params,{
    mQ %>% map(~({
        expm(.x * Delta_t)
    }))
})

mU <- with(params,{
    mP %>% map(~({
        .x[tr_names,tr_names]
    }))
})
   
mM <- with(params,{
    mP %>% map(~({
        .x[ab_names,tr_names]
    }))
}) 

# Block diagonal age transition matrix
bbD_ <- with(params,kronecker(diag(tau), D))

# Block diagonal transient matrix
# note that because we have a time homogeneous model here, we do not need to do this separately
# by cycle. but in principle, we could create a list object with each of these matrices defined
# for each time cycle. 

bbU_ <-
    with(params, {
        (1:omega %>% map(~ (mU))) %>%
            transpose() %>%
            map(~ (bdiag(.x)))
    })

K <- with(params,{vecperm(tau, omega)})

Utilde <- bbU_ %>% map(~({
    t(K) %*% bbD_ %*% K %*% .x
}))

# note that because we have a time homogeneous model here, we do not need to do this separately
# by cycle. but in principle, we could create a list object with each of these matrices defined
# for each time cycle. 

Mtilde <-
    with(params, {
        (1:omega %>% map( ~ (mM)))  %>%
            transpose()
    }) %>% 
    map(~({
        do.call(cbind,.x) 
    }))

# shoudl be alpha x tau*omega   
#c(params$alpha,params$tau * params$omega) == dim(Mtilde$natural_history)

# Age-stage multistate Markov chain transition matrix
Ptilde <-
    with(params, {
        map2(Utilde, Mtilde,  ~ ({
            rbind(cbind(.x, matrix(0, tau * omega, alpha)) ,
                  cbind(.y, diag(alpha)))
        }))
    })







# Fundamental matrix of the multistate Markov chain
# The (i, j) entry of Ñtilde is the mean time spent in state i, before eventual absorption, conditional on starting
# in state j, with i and j ranging over all the age×stage combinations.

Ntilde <-  # eq 16.
    Utilde %>% map( ~ ({
        num_states <- params$tau * (params$omega / params$Delta_t)
        solve(diag(num_states) - .x)
    }))  %>%
    map( ~ as.matrix(.x))


mZ <- # eq. 17
    with(params, {
        cbind(diag(tau * omega), matrix(0, nrow = tau * omega, ncol = alpha))
    })

#The vec operator stacks the columns of a m × n matrix into a mn × 1 column vector.

# To create reward matrices, we begin with an array. 
# rows of H correspond to health stages and the columns to ages

mH_LE <- with(params,{
    matrix(1,nrow = length(tr_names), ncol=omega)
})





# Begin function code
mH <- mH_LE
h <- vec(mH) %>% as.matrix()
mH_c = 1 - mH
h_c = vec(mH_c) %>% as.matrix()
a_ = 0.5 # fractional occupancy experienced by an individual that enters or leaves during the time step.

mBtilde1 = h %*% t(h) + (a_ * (h_c %*% t(h))) + ((1 - a_) * (h %*% t(h_c))) # EQ. 29
#mBtilde1 = with(params,a_ * (rep(1,tau*omega) %*% t(h) + h %*% t(rep(1,tau*omega))))  # Eq. 30
mCtilde1 = with(params,a_ * (rep(1,alpha) %*% t(h))) # Eq. 31
mRtilde1 = # Eq. 12
    with(params,{
        rbind(cbind(mBtilde1, matrix(0, tau * omega, alpha)) ,
          cbind(mCtilde1, diag(alpha)))
    })
mRtilde2 <- 
    mRtilde1 * mRtilde1
mRtilde3 <- 
    mRtilde1 * mRtilde1 * mRtilde1

#rho1 = ((Ntilde$natural_history %*% mZ) %*% t(Ptilde$natural_history * mRtilde1)) %*% rep(1,params$s)

mNtilde = Ntilde$natural_history
mPtilde = Ptilde$natural_history
mUtilde = Utilde$natural_history

rho1 <- ((t(mNtilde) %*% mZ)  %*% (t(as.matrix(mPtilde) * mRtilde1))) %*% rep(1,params$s)
rho2 <- t(mNtilde) %*% (mZ %*% t(as.matrix(mPtilde) * mRtilde2) %*% rep(1,params$s) + 2 * t(mUtilde * mBtilde1) %*% rho1)
rho3 <- NULL # not defined yes as need mBtilde2

trace = 
    0:(params$omega) %>% 
    map_df(~({
        as.data.frame(params$s0T %*% (t(mP$natural_history) %^% .x))
    })) %>% 
    as.matrix()

le_ = matrix(c(1,
              1 ,
              1,
              0,
              0),
            dimnames = list(c(
                c(params$tr_names,params$ab_names)
            ), c("LE")))

tr_simp <- sum((trace %*% le_) * gen_wcc(params$omega,method = "Simpson1/3"))
tr_half <- sum((trace %*% le_) * gen_wcc(params$omega,method = "half-cycle"))
tr_unadj <- sum((trace %*% le_) * gen_wcc(params$omega,method = "none"))
m_unadj <- sum(Ntilde$natural_history[,1])
m_reward <- (kronecker(diag(params$omega),t(c(1,0,0))) %*% rho1)[1] # healthy as a function of age
m_reward2 <- (kronecker(t(c(1,rep(0,params$omega-1))) ,diag(params$tau)) %*% rho1)[1,1] # healthy (from initial_age but no others)

list(m_unaj = m_unadj,
     tr_unadj = tr_unadj,
     tr_simp = tr_simp,
     tr_half = tr_half,
     m_reward =m_reward,
     m_reward2 = m_reward2) %>% cbind()
