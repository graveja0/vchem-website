#| echo: true
#| warning: false
#| code-fold: true
#| code-summary: Setup R session
library(tidyverse)
library(here)
library(glue)
library(copula)
library(knitr)
library(kableExtra)
library(ggsci)
library(MASS)
library(psych)
library(patchwork)
library(directlabels)
library(progress)
library(latex2exp)
select <- dplyr::select
options("scipen"=100, "digits"=6)
theme_set(hrbrthemes::theme_ipsum())

library(tufte)
# invalidate cache when the tufte version changes
knitr::opts_chunk$set(cache.extra = packageVersion('tufte'))
options(htmltools.dir.version = FALSE)

# knitr::purl(input = "./partial-id-intro.qmd",
#             output = "./partial-id-intro.r")



#| echo: true
#| warning: false
#| message: false
#| code-fold: true
#| code-summary: "Generate the perfect doctor data"

df <- # basic perfect doctor data (from Rubin)
  data.frame(patient = LETTERS[1:10],
            Y_1 = c(7,5,5,7,4,10,1,5,3,9),
             Y_0 = c(5,3,8,8,3,1,6,8,8,4))  %>% 
  mutate(Y_0 = c(1,6,1,8,2,1,10,6,7,8)) %>% 
  mutate(delta = Y_1 - Y_0) %>% 
  mutate(D = c(1,1,1,1,1,0,0,0,0,0)) %>% 
  mutate(Y_obs = case_when(D==1 ~ Y_1, TRUE ~ Y_0)) %>% 
  mutate(sim = 0) %>% 
  mutate(X = rnorm(10,mean = 0, sd = 1)) 

# Simulate additional data from a copula that defines the correlation b/t potential outcomes, and
# between a single covariate X and the potential outcomes. 
m <- 3 # number of columns
n <- 10000-10 # sample size
corXY <- 0.8 # corelation between potential outcomes and X
sigma <- matrix(c(1,0.6,corXY,0.6,1,corXY,corXY,corXY,1),nrow = 3) # correlation matrix

# Sample the joint CDF quantiles from a multivariate normal centered at 0
set.seed(100)
z <- mvrnorm(n,mu=rep(0, m),Sigma=sigma,empirical=T) 
u <- pnorm(z) # get the quantiles associated with each value

# Generate the variables
delta <-  # Treatment effect varies, but has mean 2 and SD 1. 
  rnorm(n,mean = 2, sd = 1)

sY_1 <- # simulated potential outcome under Tx
  qpois(u[,1],lambda = mean(df[df$D==1,]$Y_1)) + delta
sY_0 <- # simulated potential outcome under non-Tx
  qpois(u[,2],lambda=mean(df[df$D==1,]$Y_1)) + rnorm(n=n)

X <-  # a single covariate that is predictive of the outcome. 
  qnorm(u[,3],mean = 0, sd = 1)

df_sim <- # Construct the final data. 
  tibble(patient = glue("sim{1:n}")) %>% 
  mutate(Y_1 = sY_1,
         Y_0 = sY_0,
         X = X) %>% 
  mutate(delta = Y_1 - Y_0) %>% 
  mutate(D = as.integer(runif(nrow(.))<.5)) %>% 
  mutate(Y_obs = case_when(D==1 ~ Y_1, TRUE ~ Y_0)) %>% 
  mutate(sim = 1)

df <- # Bind it with the general example rows (10 rows)
  df %>% 
  bind_rows(df_sim)



#| echo: false
#| warning: false
#| message: false
#| tbl-cap: Example Data
#| label: tbl-perfdoc1
df %>% 
  filter(sim==0) %>% 
  dplyr::select(patient, Y_0, Y_1,D) %>% 
  kable() %>% 
  kableExtra::kable_styling()


#| echo: false
#| warning: false
#| message: false
#| tbl-cap: Example Data With Unit-Level Treatment Effects
#| label: tbl-perfdoc1b
df %>% 
  filter(sim==0) %>% 
  dplyr::select(patient, Y_0, Y_1,D,delta) %>% 
  kable() %>% 
  kableExtra::kable_styling()


#| fig-cap: Distribution of the Treatment Effect (delta)
#| label: fig-test
#| fig-height: 6
#| fig-width: 7
df %>% 
  filter(D==1) %>% 
  ggplot(aes(x = delta)) + geom_density() +
  labs(x = TeX("$\\Delta$"), y = "Density")+
  geom_vline(aes(xintercept = 2),lty=2,lwd=1.1, col = "darkred") + 
  annotate("text", x = 2, y = .1, hjust=-.1, label = "Avg. Treatment Effect", col = "darkred") + 
  scale_x_continuous(limits = c(min(df$delta),max(df$delta)),breaks = seq(-8,10,2)) + 
  geom_vline(aes(xintercept = 0),lty=3) 



frac_harmed <- 
  df %>% 
  filter(D==1) %>% 
  mutate(harmed = as.integer(delta<0))  %>% 
  summarise(harmed = mean(harmed)) %>% 
  pull(harmed)

y1 <- df %>% filter(D==1) %>% pull(delta)
frac_harmed <- ecdf(y1)(0)


#| fig-cap: Cumulative Distribution Function of the Treatment Effect (delta)
#| label: fig-distdelta
#| fig-height: 6
#| fig-width: 7

delta_ <- sort(df %>% filter(D==1) %>% pull(delta))
df_ <- data.frame(x = delta_, y = seq_along(delta_)/length(delta_), region = delta_ < 0 )

br0_ <- unique(sort(c(seq(0,1,0.25),ecdf(y1)(0))))
br0_l <- paste0(round(br0_,2))

figdelta <- 
  ggplot(df_, aes(x, y)) +
  geom_area(aes(fill = region)) +
  geom_line() +
  labs(x = "x = Treatment Effect", y = TeX("$Pr(\\Delta$<=x)"))+
  scale_fill_manual(values = c('lightgrey', '#C14E4295'), guide = "none") +
  annotate('text', x = 0, y = 0.08, label = 'Harmed', col = "black",size = 4, hjust=1) +
   geom_vline(aes(xintercept = 2),lty=2,lwd=1.1, col = "darkred") + 
  annotate("text", x = 1.5, y = 0.8, hjust=1, label = "Avg. Treatment Effect", col = "darkred", size = 4) +
  scale_y_continuous(breaks = br0_, labels = br0_l) +
  geom_point(data = tibble(x = 0, y = ecdf(y1)(0)), aes(x = x, y = y),size=5, pch=10) 


figdelta


#| echo: false
#| warning: false
#| message: false
#| tbl-cap: Perfect Doctor Data
#| label: tbl-perfdoc2
df %>% 
  filter(sim == 0) %>% 
  select(patient, Y_0,Y_1,D,delta) %>% 
  mutate(Y_obs = case_when(D==1 ~ Y_1, TRUE ~ Y_0)) %>% 
  mutate(Y_0 = case_when(D==1~"?",TRUE ~ paste0(Y_0))) %>% 
  mutate(Y_1 = case_when(D==0~"?",TRUE ~ paste0(Y_1))) %>% 
  mutate(delta = "?") %>% 
  select(patient,D,Y_obs,Y_0,Y_1,delta) %>% 
  mutate_all(function(x) paste0(x)) %>% 
  kable() %>% 
  kable_styling()


#| echo: false
#| warning: false
#| message: false
#| tbl-cap: Mean Outocomes by Treatment Status, and Average Treatment Effect Estimate
#| label: tbl-att

df_p <- 
  df %>% 
  group_by(D) %>% 
  summarise(Y_obs = mean(Y_obs)) %>% 
  mutate(tmp_ = 1) %>% 
  mutate(D = factor(D,levels = c(0,1), labels= c("Untreated","Treated"))) %>% 
  spread(D,Y_obs) %>% 
  mutate(hat_delta = Treated - Untreated) %>% 
  select(-tmp_)

df_p %>% 
  mutate_all(function(x) paste0(round(x,3))) %>% 
  kable() %>% 
  kable_styling()


#| echo: false
#| warning: false
#| message: false
#| fig-cap: Quantile Treatment Effects
#| label: fig-qtt
#| fig-height: 6
#| fig-width: 7

F1 <- df %>% filter(D==1) %>% pull(Y_obs)
F0 <- df %>% filter(D==0) %>% pull(Y_obs)
taus <- seq(0.05,.95,0.01)
Q1 <- quantile(F1,taus)
Q0 <- quantile(F0,taus)
QTT <- Q1 - Q0

p_qtt <- 
  tibble(tau = taus, qtt = QTT) %>% 
  ggplot(aes(x = tau, y = qtt)) + 
  geom_point() + 
  scale_y_continuous(limits = c(-2,4)) +
  labs(x = TeX("$\\tau$"),y = "QTT")

p_qtt


#| echo: false
#| message: false
#| warning: false
#|
y1 <- # Vector of treated outcomes
    df %>% filter(D==1) %>% pull(Y_obs)
F1 <- # Empirical CDF of treated outcome vector
    ecdf(y1)

y0 <- # Vector of untreated outcomes
    df %>% filter(D==0) %>% pull(Y_obs)
F0 <- # Empirical CDF of untreated outcome vector
    ecdf(y0)

delta_ <- # Vector of hypothetical treatment effects
    seq(min(df$delta), max(df$delta), length.out =1000)
y_ <- # Vector of the support of the outcome, with 100 evenly-spaced values along it. 
    seq(floor(min(y0,y1))-1,ceiling(max(y0,y1))+1, length.out = 1000)

tau <-  # Quantiles of the treatment effect distribution to obtain 
    seq(.01,.99,.01)

X1 <- # X values of treated observations
  df %>%
  filter(D==1) %>%
  pull(X)

X0 <- # X values of untreated observations
  df %>% 
  filter(D==0) %>% 
  pull(X)


#| echo: false
#| message: false
#| warning: false

F_l_wd_ <- # Lower bound
  delta_ %>% map(~{ # Outer loop is over possible values of the treatment effect
    delta <- .x
    y0 %>% # Inner loop is over y0
      map_dbl(~{
        y <- .x
        (F1(y) - F0(y-delta)) %>% max(.,0)
      })
  }) %>% 
  map_dbl(~(max(.))) # sup in formula

# Turn it into an eCDF object
  F.wd.l <- approxfun(delta_,  F_l_wd_, method = "constant", yleft = 0, yright = 1, 
            f = 0, ties = "ordered")
  class(F.wd.l) <- c("ecdf", "stepfun", class(F.wd.l))
  assign("nobs", length(delta_), envir = environment(F.wd.l))

F_u_wd_ <- # Upper bound
  delta_ %>% map(~{ # Outer loop is over possible values of the treatment effect
    delta <- .x
    y0 %>% # Inner loop is over y0
      map_dbl(~{
        y <- .x
        (F1(y) - F0(y-delta)) %>% min(.,0)
      })
  }) %>% 
  map(~(min(.))) %>% # inf in formula
  map_dbl(~(1+.x)) # add one per formula

# Turn it into an eCDF object
  F.wd.u <- approxfun(delta_,  F_u_wd_, method = "constant", yleft = 0, yright = 1, 
                      f = 0, ties = "ordered")
  class(F.wd.u) <- c("ecdf", "stepfun", class(F.wd.u))
  assign("nobs", length(delta_), envir = environment(F.wd.u))
  
df_bounds <- 
  tibble(x = delta_, ymax = F.wd.u(delta_), ymin = F.wd.l(delta_)) %>% 
            mutate(region = x<0)

qwdu <- quantile(F.wd.l, tau, type=1)
qwdl <- quantile(F.wd.u, tau, type=1)

df_wd <- 
  data.frame(tau = c(tau,tau), bound = c(qwdu,qwdl), type=c(rep("upper",length(qwdu)),rep("lower",length(qwdl)))) %>% 
  mutate(method = "Worst-Case\n[Williamson-Downs (1990)]") %>% 
  mutate(label = ifelse(type=="upper","Worst-Case\n[Williamson-Downs (1990)]","")) 



#| echo: false
#| warning: false
#| message: false
#| fig-cap: Bounds on the Cumulative Distribution Function of the Treatment Effect (delta)
#| label: fig-bounddistdelta
#| fig-height: 6
#| fig-width: 7

br_ <- sort(c(seq(0,1,0.25),F.wd.u(0),F.wd.l(0)))
br_l <- paste0(round(br_,2))

ggplot() +
  geom_line(data = df_, aes(x, y),lwd=1.5, col = "darkred") +
  annotate("text",x=4.5, y = .3, label = "True treatment effect CDF\nshown as dark red line.",hjust=0,size = 4) +
  labs(x = "x = Treatment Effect", y = TeX("$Pr(\\Delta$<=x)"))+
  scale_fill_manual(values = c('lightgrey', '#C14E4295'), guide = "none") +
  #geom_line(data = tibble(x = delta_, y = F.wd.l(delta_))) +
  #geom_line(data = tibble(x = delta_, y = F.wd.u(delta_))) +
  geom_ribbon(data = df_bounds , 
              aes(x = x, ymin = ymin, ymax=ymax, fill = region), alpha = 0.5)  +
  geom_line(data = df_bounds, aes(x = x, y = ymin), alpha = 1,lwd=1.25) +
  geom_line(data = df_bounds, aes(x = x, y = ymax), alpha = 1,lwd=1.25) + 
  scale_y_continuous(breaks = br_, labels = br_l) +
  geom_point(data = tibble(x = 0, y = F.wd.u(0)), aes(x = x, y = y),size=5, pch=10) +
  geom_point(data = tibble(x = 0, y = F.wd.l(0)), aes(x = x, y = y),size=5, pch=10)
  


#| echo: true
#| message: false
#| warning: false
#| code-fold: true
#| label: fig-bounds-wd
#| fig-cap: Worst-Case Treatment Effect Distribution Bounds
#| fig-width: 10
#| fig-height: 8

#col_scheme <- c("#FF5733" ) # Williamson and Downs
col_scheme <- c("black","#FF5733")

df_wd  %>% 
  mutate(method = factor(method, levels = c("Frandsen & Lefgrens (2021)"     ,       "Frandsen & Lefgrens (2021) Covariates", "Worst-Case\n[Williamson-Downs (1990)]"))) %>% 
  mutate(type = paste0(type,method)) %>% 
  ggplot(aes(x = tau, y = bound, group = type, colour = method)) + 
  geom_line(lwd=1.25) + 
  #geom_hline(aes(yintercept = 2), colour = "darkred",lwd=1.25) + 
  scale_color_manual(values=col_scheme) + 
  geom_dl(method = list("last.points",hjust=1),aes(label = label)) +
  scale_x_continuous(breaks = seq(0,1,0.1)) +
  theme(legend.position = "none")  +
  geom_point(data =  tibble(tau = taus, qtt = QTT, type = "True Values", method = "True Values"), aes(x = tau, y = qtt )) 
  #annotate("text",x = 0.5, y = 2, label = "True Treatment Effect",vjust=-1,colour = "darkred") 




i1 <- # Index of treated observations
  df %>% 
  mutate(i = row_number()) %>% 
  filter(D==1) %>% 
  pull(i)

i0 <- # index of untreated observations
  df %>% 
  mutate(i = row_number()) %>% 
  filter(D==0) %>% 
  pull(i)


F0_fl <-# ecdf for untreated defined across the support of y
  y_ %>% map_dbl(~({
    .y = .x
    mean(df$Y_obs[i0] <= .y)
  })) %>% BMisc::makeDist(y_, .)


F1_fl <-# ecdf for treated defined across the support of y
  y_ %>% map_dbl(~({
    .y = .x
    mean(df$Y_obs[i1] <= .y )
  })) %>% BMisc::makeDist(y_, .)

# Construct bounds

pb <- progress_bar$new(total = length(delta_))
F_l_fl_ <- # Lower bound
  delta_ %>% map(~{ # outer loop is over possible values of the treatment effect
    delt <- .x
    pb$tick()
    y0 %>% # inner loop is over observed outcome values among untreated
      map_dbl(~({
        y0_ = .x 
        (F1_fl(y0_ + delt)-F0_fl(y0_))/(1-F0_fl(y0_))
      }))  %>% 
      map_dbl(~(max(0,.x,na.rm=TRUE))) %>% 
      mean(.)
  }) %>% 
  unlist()

pb <- progress_bar$new(total = length(delta_))
F_u_fl_ <-  # Upper bound
  delta_ %>% map(~{ # outer loop is over possible values of the treatment effect
    delt <- .x
    pb$tick()
    y0 %>% # inner loop is over observed outcome values among untreated 
      map_dbl(~({
        y0_ = .x 
        (F1_fl(y0_+delt))/(F0_fl(y0_))
      }))  %>% 
      map_dbl(~(min(1,.x,na.rm=TRUE))) %>% 
      mean(.)
  }) %>% 
  unlist()


# Convert into an ECDF
F.fl.l <- approxfun(delta_,  F_l_fl_, method = "constant", yleft = 0, yright = 1, 
                    f = 0, ties = "ordered")
class(F.fl.l) <- c("ecdf", "stepfun", class(F.fl.l))
assign("nobs", length(delta_), envir = environment(F.fl.l))

F.fl.u <- approxfun(delta_,  F_u_fl_, method = "constant", yleft = 0, yright = 1, 
                    f = 0, ties = "ordered")
class(F.fl.u) <- c("ecdf", "stepfun", class(F.fl.u))
assign("nobs", length(delta_), envir = environment(F.fl.u))

df_bounds_fl <- 
  tibble(x = delta_, ymax = F.fl.u(delta_), ymin = F.fl.l(delta_)) %>% 
            mutate(region = x<0)

qflu <- quantile(F.fl.l, tau, type=1)
qfll <- quantile(F.fl.u, tau, type=1)

df_fl <- 
  data.frame(tau = c(tau,tau), bound = c(qflu,qfll), type=c(rep("upper",length(qflu)),rep("lower",length(qfll)))) %>% 
  mutate(method = "Frandsen & Lefgrens (2021)") %>% 
  mutate(label = ifelse(type=="upper","Frandsen & Lefgrens (2021)","")) 



#| echo: false
#| warning: false
#| message: false
#| fig-cap: Frandsen-Lefgren (2021) Bounds on the Cumulative Distribution Function of the Treatment Effect (delta)
#| label: fig-flbounddistdelta
#| fig-height: 6
#| fig-width: 7


br_wd <- sort(c(seq(0,1,0.25),F.wd.u(0),F.wd.l(0),F.fl.u(0),F.fl.l(0)))
br_wd_l <- paste0(round(br_wd,2))

ggplot() +
  geom_line(data = df_, aes(x, y),lwd=1.5, col = "darkred") +
  annotate("text",x=4.5, y = .3, label = "True treatment effect CDF\nshown as dark red line.",hjust=0,size = 4) +
   labs(x = "x = Treatment Effect", y = TeX("$Pr(\\Delta$<=x)"))+
  scale_fill_manual(values = c('lightgrey', '#C14E4295'), guide = "none") +
  #geom_line(data = tibble(x = delta_, y = F.wd.l(delta_))) +
  #geom_line(data = tibble(x = delta_, y = F.wd.u(delta_))) +
  geom_ribbon(data = df_bounds_fl , 
              aes(x = x, ymin = ymin, ymax=ymax, fill = region), alpha = 0.5) +
  hrbrthemes::theme_ipsum() + 
  geom_line(data = df_bounds, aes(x = x, y = ymin), alpha = 0.5) +
  geom_line(data = df_bounds, aes(x = x, y = ymax), alpha = 0.5) +
  geom_line(data = df_bounds_fl, aes(x = x, y = ymin), alpha = 1,lwd=1.25) +
  geom_line(data = df_bounds_fl, aes(x = x, y = ymax), alpha = 1,lwd=1.25) +
    scale_y_continuous(breaks = br_wd, labels = br_wd_l) + 
  geom_point(data = tibble(x = 0, y = F.fl.u(0)), aes(x = x, y = y),size=5, pch=10) +
  geom_point(data = tibble(x = 0, y = F.fl.l(0)), aes(x = x, y = y),size=5, pch=10)
  



#| echo: true
#| message: false
#| warning: false
#| code-fold: true
#| label: fig-bounds-fl
#| fig-cap: Fransden-Lefgren Treatment Effect Distribution Bounds
#| fig-width: 10
#| fig-height: 8

#col_scheme <- c("#FF5733" ) # Williamson and Downs
col_schemefl <- c("#FF5733","darkred","#7DC4CC")

df_wd  %>% 
  bind_rows(df_fl) %>% 
  mutate(method = factor(method, levels = c("Frandsen & Lefgrens (2021)"     ,       "Frandsen & Lefgrens (2021) Covariates", "Worst-Case\n[Williamson-Downs (1990)]"))) %>% 
  mutate(type = paste0(type,method)) %>% 
  ggplot(aes(x = tau, y = bound, group = type, colour = method)) + 
  geom_line(lwd=1.25) + 
  #geom_hline(aes(yintercept = 2), colour = "darkred",lwd=1.25) + 
  scale_color_manual(values=col_schemefl) + 
  geom_dl(method = list("last.points",hjust=1),aes(label = label)) +
  scale_x_continuous(breaks = seq(0,1,0.1)) +
  theme(legend.position = "none")  +
  geom_point(data =  tibble(tau = taus, qtt = QTT, type = "True Values", method = "True Values"), aes(x = tau, y = qtt )) 
  #annotate("text",x = 0.5, y = 2, label = "True Treatment Effect",vjust=-1,colour = "darkred") 



#| echo: false
#| message: false
#| warning: false
#| code-fold: false

i1 <-  # index of treated observations
  df %>% 
  mutate(i = row_number()) %>% 
  filter(D==1) %>% 
  pull(i)

i0 <-  # index of untreated observations
  df %>% 
  mutate(i = row_number()) %>% 
  filter(D==0) %>% 
  pull(i)

X1 <- # X values of treated observations
  df %>%
  filter(D==1) %>%
  pull(X)

X0 <- # X values of untreated observations
  df %>% 
  filter(D==0) %>% 
  pull(X)


taus0 <- sort(ecdf(y0)(y0))
quantreg_fit <- quantreg::rq(Y_obs ~ X, tau = taus0, data = df %>% filter(D==0))
pvals <- predict(quantreg_fit,newdata = data.frame(X=X0), stepfun=TRUE)
hatF0X <- pvals %>% map(~(.x(taus0))) %>% 
  map(~(BMisc::makeDist(.x,taus0)))


taus1 <- sort(ecdf(y1)(y1))
quantreg_fit <- quantreg::rq(Y_obs ~ X, tau = taus0, data = df %>% filter(D==1))
pvals <- predict(quantreg_fit,newdata = data.frame(X=X0), stepfun=TRUE)
hatF1X <- pvals %>% map(~(.x(taus0))) %>% 
  map(~(BMisc::makeDist(.x,taus0)))


pb <- progress_bar$new(total = length(delta_))
F_l_flX_ <- 
  delta_ %>% map(~{ # outer loop is over possible values of the treatment effect
    delt <- .x
    pb$tick()
      map2_dbl(y0,X0,~({
        y0_ = .x 
        xx_ = .y
        i = which(X0==xx_); i
        F1_fl_tmp <- hatF1X[[i]]
        F0_fl_tmp <- hatF0X[[i]]
        (F1_fl_tmp(y0_ + delt)-F0_fl_tmp(y0_))/(1-F0_fl_tmp(y0_))
      }))  %>% 
      map_dbl(~(max(0,.x,na.rm=TRUE))) %>% 
      mean(.)
  }) %>% 
  unlist()
# Make it into a proper eCDF object.
F.flX.l <- approxfun(delta_,  F_l_flX_, method = "constant", yleft = 0, yright = 1, 
                    f = 0, ties = "ordered")
class(F.flX.l) <- c("ecdf", "stepfun", class(F.flX.l))
assign("nobs", length(delta_), envir = environment(F.flX.l))

pb <- progress_bar$new(total = length(delta_))
F_u_flX_ <- 
  delta_ %>% 
  map(~{ # outer loop is over possible values of the treatment effect
    delt <- .x
    pb$tick()
    
      map2_dbl(y0,X0,~({
        y0_ = .x 
        xx_ = .y
        i = which(X0==xx_); i
        F1_fl_tmp <- hatF1X[[i]]
        F0_fl_tmp <- hatF0X[[i]]
        
        (F1_fl_tmp(y0_+delt))/(F0_fl_tmp(y0_))
      }))  %>% 
      map_dbl(~(min(1,.x,na.rm=TRUE))) %>% 
      mean(.)
  }) %>% 
  unlist()
# Make it into a proper eCDF object. 
F.flX.u <- approxfun(delta_,  F_u_flX_, method = "constant", yleft = 0, yright = 1, 
                    f = 0, ties = "ordered")
class(F.flX.u) <- c("ecdf", "stepfun", class(F.flX.u))
assign("nobs", length(delta_), envir = environment(F.flX.u))

df_bounds_flX <- 
  tibble(x = delta_, ymax = F.flX.u(delta_), ymin = F.flX.l(delta_)) %>% 
            mutate(region = x<0)

qflXu <- quantile(F.flX.l, tau, type=1)
qflXl <- quantile(F.flX.u, tau, type=1)

df_flX <- 
  data.frame(tau = c(tau,tau), bound = c(qflXu,qflXl), type=c(rep("upper",length(qflXu)),rep("lower",length(qflXl)))) %>% 
  mutate(method = "Frandsen & Lefgrens (2021)\nWith Covariates") %>% 
  mutate(label = ifelse(type=="upper","Frandsen & Lefgrens (2021)\nWith Covariates","")) 



#| echo: false
#| warning: false
#| message: false
#| fig-cap: Frandsen-Lefgren (2021) Bounds on the Cumulative Distribution Function of the Treatment Effect, With Covariates
#| label: fig-flboundXdistdelta
#| fig-height: 6
#| fig-width: 7


br_wdX <- sort(c(seq(0,1,0.25),F.wd.u(0),F.wd.l(0),F.fl.u(0),F.fl.l(0),F.flX.u(0),F.flX.l(0)))
br_wdX_l <- paste0(round(br_wdX,2))

ggplot() +
  geom_line(data = df_, aes(x, y),lwd=1.5, col = "darkred") +
  annotate("text",x=4.5, y = .3, label = "True treatment effect CDF\nshown as dark red line.",hjust=0,size = 4) +
   labs(x = "x = Treatment Effect", y = TeX("$Pr(\\Delta$<=x)"))+
  scale_fill_manual(values = c('lightgrey', '#C14E4295'), guide = "none") +
  #geom_line(data = tibble(x = delta_, y = F.wd.l(delta_))) +
  #geom_line(data = tibble(x = delta_, y = F.wd.u(delta_))) +
  geom_ribbon(data = df_bounds_flX , 
              aes(x = x, ymin = ymin, ymax=ymax, fill = region), alpha = 0.5) +
  geom_line(data = df_bounds, aes(x = x, y = ymin), alpha = 0.5) +
  geom_line(data = df_bounds, aes(x = x, y = ymax), alpha = 0.5) +
  geom_line(data = df_bounds_fl, aes(x = x, y = ymin), alpha = 0.5) +
  geom_line(data = df_bounds_fl, aes(x = x, y = ymax), alpha = 0.5) +
  geom_line(data = df_bounds_flX, aes(x = x, y = ymin), alpha = 1,lwd=1.25) +
  geom_line(data = df_bounds_flX, aes(x = x, y = ymax), alpha = 1,lwd=1.25) +
  geom_point(data = tibble(x = 0, y = F.flX.u(0)), aes(x = x, y = y),size=5, pch=10) +
  geom_point(data = tibble(x = 0, y = F.flX.l(0)), aes(x = x, y = y),size=5, pch=10) + 
  scale_y_continuous(breaks = br_wdX, labels = br_wdX_l) 
  



#| echo: true
#| message: false
#| warning: false
#| code-fold: true
#| label: fig-bounds-flX
#| fig-cap: Fransden-Lefgren Treatment Effect Distribution Bounds With Covariates
#| fig-width: 10
#| fig-height: 8

#col_scheme <- c("#FF5733" ) # Williamson and Downs
col_schemeflX <- c("#FF5733","darkred","#7DC4CC","#0A2E36")

df_wd  %>% 
  bind_rows(df_fl) %>% 
  bind_rows(df_flX) %>% 
  mutate(method = factor(method, levels = c("Frandsen & Lefgrens (2021)"     ,       "Frandsen & Lefgrens (2021) Covariates", "Worst-Case\n[Williamson-Downs (1990)]"))) %>% 
  mutate(type = paste0(type,method)) %>% 
  ggplot(aes(x = tau, y = bound, group = type, colour = method)) + 
  geom_line(lwd=1.25) + 
  #geom_hline(aes(yintercept = 2), colour = "darkred",lwd=1.25) + 
  scale_color_manual(values=col_schemeflX) + 
  geom_dl(method = list("last.points",hjust=1),aes(label = label)) +
  scale_x_continuous(breaks = seq(0,1,0.1)) +
  theme(legend.position = "none")  +
  geom_point(data =  tibble(tau = taus, qtt = QTT, type = "True Values", method = "True Values"), aes(x = tau, y = qtt )) 
  #annotate("text",x = 0.5, y = 2, label = "True Treatment Effect",vjust=-1,colour = "darkred") 



#| echo: false
#| warning: false
#| message: false
#| fig-cap: Decision Thresholds Based on the Maximum Allowable Fraction of the Population Harmed by Treatment
#| label: fig-fllambda1
#| fig-height: 6
#| fig-width: 7

lambda1 <- 0.4
lambda2 <- 0.15
ggplot() +
 labs(x = "x = Treatment Effect", y = TeX("$Pr(\\Delta$<=x)"))+
  geom_ribbon(data = df_bounds_flX , 
              aes(x = x, ymin = ymin, ymax=ymax, fill = region), alpha = 0.5) +
  scale_fill_manual(values = c('lightgrey', '#C14E4295'), guide = "none") +
  geom_line(data = df_bounds_flX, aes(x = x, y = ymin), alpha = 1) +
  geom_line(data = df_bounds_flX, aes(x = x, y = ymax), alpha = 1) +
  #scale_y_continuous() + 
  geom_segment(aes(x = min(df$delta), xend = quantile(F.flX.u,lambda1),y=lambda1,yend=lambda1),lwd=1.25) +
  geom_segment(aes(x = min(df$delta), xend = quantile(F.flX.u,lambda2),y=lambda2,yend=lambda2),lwd=1.25) +
  geom_segment(aes(x = min(df$delta), xend = 0, y = F.flX.u(0), yend = F.flX.u(0)),lwd=1.25, col = "darkblue") +
  scale_y_continuous(breaks = c(lambda1,lambda2,F.flX.u(0)), labels = c(TeX("$\\lambda_1 = 0.4$"),TeX("$\\lambda_2 = 0.15$"),TeX("$F_{\\Delta}^U(0)$"))) +
  theme(axis.text.y = element_text(size=12),panel.grid.minor.y  = element_blank())




#| echo: false
#| warning: false
#| message: false
#| fig-cap: Decision Thresholds Based on the Minimum Treatment Effect
#| label: fig-fllambda2
#| fig-height: 6
#| fig-width: 7


lambda3 <- 2.5
lambda4 <- 5
ggplot() +
  labs(x = "x = Treatment Effect", y = TeX("$Pr(\\Delta$<=x)"))+
  geom_ribbon(data = df_bounds_flX , 
              aes(x = x, ymin = ymin, ymax=ymax, fill = region), alpha = 0.5) +
  scale_fill_manual(values = c('lightgrey', '#C14E4295'), guide = "none") +
  geom_line(data = df_bounds_flX, aes(x = x, y = ymin), alpha = 1) +
  geom_line(data = df_bounds_flX, aes(x = x, y = ymax), alpha = 1) +
  geom_segment(aes(x = lambda3, xend = lambda3, y = 0, yend = F.flX.l(lambda3)),lwd=1.25) +
  geom_segment(aes(x = lambda4, xend = lambda4, y = 0, yend = F.flX.l(lambda4)),lwd=1.25) +
  geom_segment(aes(x = quantile(F.flX.l,0.5),xend = quantile(F.flX.l,0.5), y = 0, yend = F.flX.l(quantile(F.flX.l,0.5))),lwd=1.25, col = "darkblue") +
  geom_segment(aes(x = min(df$delta), xend = quantile(F.flX.l,0.5), y = 0.5, yend = 0.5),lwd=1.25, col = "darkblue")  +
  scale_x_continuous(breaks = c(lambda3,lambda4,quantile(F.flX.l,0.5)), labels = c(TeX("$\\lambda_3 = 2.5$"),TeX("$\\lambda_4 = 5$"),TeX("${F^{-1U}_{\\Delta}(0.5)}$")),guide = guide_axis(n.dodge = 3)) +
  scale_y_continuous(breaks = c(F.flX.l(lambda3),F.flX.l(lambda4),0.5), labels = round(c(F.flX.l(lambda3),F.flX.l(lambda4),0.5),3)) +
  theme(axis.text.y = element_text(size=12),axis.text.x = element_text(size=12),panel.grid.minor.y  = element_blank()) 
 





## knitr::purl(input = here("blog/drafts/partial-id-intro/partial-id-intro.qmd"),
##             output = here("blog/drafts/partial-id-intro/partial-id-intro.r"))

