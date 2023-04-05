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



params <- 
  list(
    N = 2e3,
    sigma_sq_X = 1.0,
    sigma_sq_epsilon = 0.3,
    delta = 1,
    sigma_sq_delta = 1
  )
params <- with(params,
               modifyList(params,list(
                 r_squared = 1 - sigma_sq_epsilon,
                 beta = sqrt(1 - sigma_sq_epsilon),
                 alpha_1 = -0.5 * (sigma_sq_delta / sigma_sq_epsilon),
                 Sigma = matrix(c(sigma_sq_X,0,0,sigma_sq_epsilon),
                                byrow=TRUE, nrow = 2, ncol = 2))))
params <- with(params,
               modifyList(params,list(
                 sigma_sq_eta = sigma_sq_delta - (alpha_1^2)*sigma_sq_epsilon
               )))

gen_data <- function(params) {
  
  with(params, 
       mvrnorm(n = N, mu = c(0,0), Sigma = Sigma)) %>% 
    data.frame() %>% 
    as_tibble() %>% 
    set_names(c("X","epsilon")) %>% 
    mutate(eta_i = rnorm(nrow(.),mean = 0,sd = params$sigma_sq_eta)) %>% 
    mutate(alpha_0 = params$delta) %>% 
    mutate(delta_i = alpha_0 + params$alpha_1 * epsilon + eta_i) %>% 
    mutate(Y_i0 = params$beta * X + epsilon) %>% 
    mutate(Y_i1 = Y_i0 + delta_i) %>% 
    mutate(random = runif(nrow(.))) %>% 
    mutate(D = as.integer(row_number()<=(params$N)/2)) %>% 
    arrange(random) %>% 
    select(-random) %>% 
    mutate(Y = D * Y_i1 + (1 - D) * Y_i0) %>% 
    select(Y,Y_0 = Y_i0,Y_1 = Y_i1,D,X,delta_i) %>% 
    mutate(delta = Y_1 - Y_0)
}

set.seed(123)
df <- params %>% gen_data() %>% 
  mutate(sim = ifelse(row_number() %in% 1:10, 0,1)) %>% 
  mutate(Y_obs = Y) %>% 
  mutate(patient = ifelse(sim==0,LETTERS[row_number()],NA))  

hist(df$delta)

df %>% head(n=10) %>% kable() %>% 
  kable_styling()





#| echo: false
#| warning: false
#| message: false
#| tbl-cap: Example Data
#| label: tbl-perfdoc1
df %>% 
  filter(sim==0) %>% 
  head(n = 10) %>% 
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

quantile(df$delta,0.9)


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

col_scheme <- c("#7DC4CC","#0A2E36","#FF5733")



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
  geom_ribbon(data = df_bounds , 
              aes(x = x, ymin = ymin, ymax=ymax, fill = region), alpha = 0.5)  +
  geom_line(data = df_bounds, aes(x = x, y = ymin), alpha = 1,lwd=1.25) +
  geom_line(data = df_bounds, aes(x = x, y = ymax), alpha = 1,lwd=1.25) + 
  scale_y_continuous(breaks = br_, labels = br_l) +
  geom_point(data = tibble(x = 0, y = F.wd.u(0)), aes(x = x, y = y),size=5, pch=10) +
  geom_point(data = tibble(x = 0, y = F.wd.l(0)), aes(x = x, y = y),size=5, pch=10)
  


#| echo: false
#| warning: false
#| message: false
#| fig-cap: Worst-Case Bounds on the Average Treatment Effect on the Treated
#| label: fig-boundwdatt
#| fig-height: 3
#| fig-width: 7
#|

## attwd <-
##   tibble(delta = delta_,
##        FU = F.wd.l(delta_),
##        FL = F.wd.u(delta_)) %>%
##   mutate(diffU = c(0,diff(FU)),
##          diffL = c(0,diff(FL))) %>%
##   filter((FU>0|FL>0)&(FU<1|FL<1)) %>%
##   mutate(diffdelta = c(0,diff(delta))) %>%
##   arrange(delta) %>%
##   summarise(upper = sum(delta*diffU),
##             lower = sum(delta*diffL)) %>%
##   mutate(method = "Worst-Case\n[Williamson-Downs (1990)]" ) %>%
##   mutate(size =1) %>%
##   mutate(method = factor(method, levels = c("Frandsen & Lefgrens (2021)","Frandsen & Lefgrens (2021) Covariates", "Worst-Case\n[Williamson-Downs (1990)]")))
## 
## 
## attwd %>%
##   mutate(y = 0) %>%
##   mutate(method = factor(method, levels = c("Frandsen & Lefgrens (2021)","Frandsen & Lefgrens (2021) Covariates", "Worst-Case\n[Williamson-Downs (1990)]")))  %>%
##   ggplot() +
##   geom_errorbar(aes(y = 0, xmin = lower, xmax = upper,colour = method),width = 0.8,lwd = 2) +
##   scale_y_continuous(limits = c(-1,1), breaks = NULL) +
##   geom_point(aes(x = mean(df$delta), y = 0), colour = "darkred", size=8, alpha = 0.5) +
##   geom_vline(aes(xintercept = 0), lty=3, col = "darkgrey") +
##   labs(x = "ATT", y = "") +
##   scale_colour_manual(values = col_scheme[3],name="") +
##   theme(legend.position = "bottom")
## 
## 


#| echo: true
#| message: false
#| warning: false
#| code-fold: true
#| label: fig-bounds-wd
#| fig-cap: Worst-Case Treatment Effect Distribution Bounds
#| fig-width: 10
#| fig-height: 8

#col_scheme <- c("#FF5733" ) # Williamson and Downs

df_wd  %>% 
  mutate(method = factor(method, levels = c("Frandsen & Lefgrens (2021)"     ,       "Frandsen & Lefgrens (2021) Covariates", "Worst-Case\n[Williamson-Downs (1990)]"))) %>% 
  mutate(type = paste0(type,method)) %>% 
  ggplot(aes(x = tau, y = bound, group = type, colour = method)) + 
  geom_line(lwd=1.5,colour = "#FF5733") + 
  geom_dl(method = list("last.points",hjust=1),aes(label = label)) +
  scale_x_continuous(breaks = seq(0,1,0.1)) +
  theme(legend.position = "none")  +
  scale_colour_manual(values = col_scheme[3]) +
  labs(x = TeX("$$\\tau$$"),y=TeX("QoTT($$\\tau$$)"))
  # geom_point(data =  tibble(tau = taus, qtt = QTT, type = "True Values", method = "True Values"), aes(x = tau, y = qtt ),
  #            col="darkred") 





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
  



#| echo: false
#| warning: false
#| message: false
#| fig-cap: Worst-Case Bounds on the Average Treatment Effect on the Treated, by Method
#| label: fig-boundwdattfl
#| fig-height: 3
#| fig-width: 7


attfl <- tibble(delta = delta_,
       FU = F.fl.l(delta_),
       FL = F.fl.u(delta_)) %>% 
  mutate(diffU = c(0,diff(FU)),
         diffL = c(0,diff(FL))) %>% 
  summarise(upper = sum(delta*diffU),
            lower = sum(delta*diffL)) %>% 
  mutate(method = "Frandsen & Lefgrens (2021)" ) %>% 
  mutate(size = 2) 

attwd <- tibble(delta = delta_,
       FU = F.wd.l(delta_),
       FL = F.wd.u(delta_)) %>% 
  mutate(diffU = c(0,diff(FU)),
         diffL = c(0,diff(FL))) %>% 
  summarise(upper = sum(delta*diffU),
            lower = sum(delta*diffL)) %>% 
  mutate(method = "Worst-Case\n[Williamson-Downs (1990)]" ) %>% 
  mutate(size =1)

attfl %>% 
  #bind_rows(attwd) %>% 
  mutate(y = 0) %>% 
  mutate(method = factor(method, levels = c("Frandsen & Lefgrens (2021)","Frandsen & Lefgrens (2021) Covariates", "Worst-Case\n[Williamson-Downs (1990)]")))  %>% 
  ggplot() + 
  geom_errorbar(aes(y = 0, xmin = lower, xmax = upper,colour = method),width = 0.8,lwd = 2) + 
  scale_y_continuous(limits = c(-1,1), breaks = NULL) + 
  geom_point(aes(x = mean(df$delta[df$D==1]), y = 0), colour = "darkred", size=8, alpha = 0.5) + 
  geom_vline(aes(xintercept = 0), lty=3, col = "darkgrey") + 
  labs(x = "ATT", y = "") +
  scale_colour_manual(values = col_scheme[c(1,3)],name="") +
  theme(legend.position = "bottom")



#| echo: true
#| message: false
#| warning: false
#| code-fold: true
#| label: fig-bounds-fl
#| fig-cap: Fransden-Lefgren Treatment Effect Distribution Bounds
#| fig-width: 10
#| fig-height: 8

#col_scheme <- c("#FF5733" ) # Williamson and Downs

df_wd  %>% 
  bind_rows(df_fl) %>% 
  mutate(method = factor(method, levels = c("Frandsen & Lefgrens (2021)"     ,       "Frandsen & Lefgrens (2021) Covariates", "Worst-Case\n[Williamson-Downs (1990)]"))) %>% 
  mutate(type = paste0(type,method)) %>% 
  ggplot(aes(x = tau, y = bound, group = type)) + 
  geom_line(lwd=1.25,aes(colour = method)) + 
  scale_color_manual(values=col_scheme[c(1,3)]) + 
  geom_dl(method = list("last.points",hjust=1),aes(label = label)) +
  scale_x_continuous(breaks = seq(0,1,0.1)) +
  theme(legend.position = "none")  +
   labs(x = TeX("$$\\tau$$"),y=TeX("QoTT($$\\tau$$)"))



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
  



#| echo: false
#| warning: false
#| message: false
#| fig-cap: Worst-Case Bounds on the Average Treatment Effect on the Treated, by Method
#| label: fig-boundwdattflX
#| fig-height: 3
#| fig-width: 7


attfl <- tibble(delta = delta_,
       FU = F.fl.l(delta_),
       FL = F.fl.u(delta_)) %>% 
  mutate(diffU = c(0,diff(FU)),
         diffL = c(0,diff(FL))) %>% 
  summarise(upper = sum(delta*diffU),
            lower = sum(delta*diffL)) %>% 
  mutate(method = "Frandsen & Lefgrens (2021)" ) %>% 
  mutate(size = 2) 

attflX <- tibble(delta = delta_,
       FU = F.flX.l(delta_),
       FL = F.flX.u(delta_)) %>% 
  mutate(diffU = c(0,diff(FU)),
         diffL = c(0,diff(FL))) %>% 
  summarise(upper = sum(delta*diffU),
            lower = sum(delta*diffL)) %>% 
  mutate(method = "Frandsen & Lefgrens (2021) Covariates" ) %>% 
  mutate(size = 2) 

attwd <- tibble(delta = delta_,
       FU = F.wd.l(delta_),
       FL = F.wd.u(delta_)) %>% 
  mutate(diffU = c(0,diff(FU)),
         diffL = c(0,diff(FL))) %>% 
  summarise(upper = sum(delta*diffU),
            lower = sum(delta*diffL)) %>% 
  mutate(method = "Worst-Case\n[Williamson-Downs (1990)]" ) %>% 
  mutate(size =1)

attfl %>% 
  bind_rows(attflX) %>% 
  bind_rows(attwd) %>% 
  mutate(y = 0) %>% 
  mutate(method = factor(method, levels = c("Frandsen & Lefgrens (2021)","Frandsen & Lefgrens (2021) Covariates", "Worst-Case\n[Williamson-Downs (1990)]")))  %>% 
  ggplot() + 
  geom_errorbar(aes(y = 0, xmin = lower, xmax = upper,colour = method),width = 0.8,lwd = 2) + 
  scale_y_continuous(limits = c(-1,1), breaks = NULL) + 
  geom_point(aes(x = 2, y = 0), colour = "darkred", size=8, alpha = 0.5) + 
  geom_vline(aes(xintercept = 0), lty=3, col = "darkgrey") + 
  labs(x = "ATT", y = "") +
  scale_colour_manual(values = col_scheme,name="") +
  theme(legend.position = "bottom")



#| echo: true
#| message: false
#| warning: false
#| code-fold: true
#| label: fig-bounds-flX
#| fig-cap: Fransden-Lefgren Treatment Effect Distribution Bounds With Covariates
#| fig-width: 10
#| fig-height: 8

df_wd  %>% 
  bind_rows(df_fl) %>% 
  bind_rows(df_flX) %>% 
  mutate(method = factor(method, levels = c("Frandsen & Lefgrens (2021)"     ,       "Frandsen & Lefgrens (2021)\nWith Covariates", "Worst-Case\n[Williamson-Downs (1990)]"))) %>% 
  mutate(type = paste0(type,method)) %>% 
  ggplot(aes(x = tau, y = bound, group = type)) + 
  geom_line(lwd=1.25,aes(colour = method)) + 
  scale_color_manual(values=col_scheme) + 
  geom_dl(method = list("last.points",hjust=1),aes(label = label,colour = method)) +
  scale_x_continuous(breaks = seq(0,1,0.1)) +
  theme(legend.position = "none")  +
  labs(x = TeX("$$\\tau$$"),y=TeX("QoTT($$\\tau$$)"))
  #annotate("text",x = 0.5, y = 2, label = "True Treatment Effect",vjust=-1,colour = "darkred") 




# Manski Bounds

bounds <- list(
    WC = list(
        EY1 = list(
            Z0_u = "maxY",
            Z0_l = "minY",
            Z1_u = "Y_1",
            Z1_l = "Y_1"
        ),
        EY0 = list(
            Z0_u = "Y_0",
            Z0_l = "Y_0",
            Z1_u = "maxY",
            Z1_l = "minY"
        )
    ),
    MTS = list(
        EY1 = list(
            Z0_u = "maxY",
            Z0_l = "Y_1",# Had the untreated group received the treatment (lower values better) they would have no better outcomes than the treated group
            Z1_u = "Y_1",
            Z1_l = "Y_1"
        ),
        EY0 = list(
            Z0_u = "Y_0",
            Z0_l = "Y_0",
            Z1_u = "Y_0", # Had the untreated group received the treatment (lower values better) they would have no worse outcomes than in the untreated group
            Z1_l = "minY"
        )
    ),
    MTR = list(
        EY1 = list(
            Z0_u = "maxY",
            Z0_l = "Y_0", # control units cant' get any worse if they get treated
            Z1_u = "Y_1",
            Z1_l = "Y_1"
        ),
        EY0 = list(
            Z0_u = "Y_0",
            Z0_l = "Y_0",
            Z1_u = "Y_1",  # treated units cant get any better if they get control
            Z1_l = "minY"
        )
    )
)


df_bounds_ <-
    df %>%
    mutate(W = ifelse(Y_1>Y_0,1,0)) %>%
    mutate(Y = ifelse(W==1,Y_1,Y_0)) %>%
    mutate(Z= D) %>%
    mutate(maxY = max(Y),
           minY = min(Y)) %>%
    group_by(W) %>%
    summarise_at(vars(Y,maxY,minY),mean)  %>%
    mutate(PX = 0.5) %>%
    #mutate(id = "WC") %>%
    gather(measure,value,-W) %>%
    unite(measure,measure,W) %>%
    spread(measure,value) %>%
    rename(maxY = maxY_0,
           minY = minY_0) %>%
    select(-minY_1,-maxY_1) %>%
    mutate(test = 3) %>%
    uncount(test) %>%
    mutate(assumption = c("WC","MTR","MTS"))

bounds_lut <-
    bounds %>% map_dfr(~(.x %>% data.frame())) %>%
    mutate(assumption = names(bounds)) %>%
    gather(measure,quantity,-assumption) %>%
    separate(measure, into = c("potential_outcome","bound"),sep ="\\.")

df_bounds_manski <-
    bounds_lut %>%
    left_join(
        df_bounds_ %>%
            gather(quantity,value,-PX_0,-PX_1,-assumption), c("assumption","quantity")
    ) %>%
    select(-quantity) %>%
    unite(bound,potential_outcome,bound) %>%
    spread(bound,value) %>%
    mutate(LB = (PX_0 * EY1_Z0_l + PX_1 * EY1_Z1_l) - (PX_0 * EY0_Z0_u + PX_1 * EY0_Z1_u),
           UB = (PX_0 * EY1_Z0_u + PX_1 * EY1_Z1_u) - (PX_0 * EY0_Z0_l + PX_1 * EY0_Z1_l)) %>%
    mutate(assumption = factor(assumption, levels = c("WC","MTR","MTS"))) %>% 
    arrange(assumption)



#| echo: false
#| warning: false
#| message: false
#| fig-cap: Worst-Case Bounds on the Average Treatment Effect on the Treated, by Method
#| label: fig-boundmanski
#| fig-height: 6
#| fig-width: 7


attfl <- tibble(delta = delta_,
       FU = F.fl.l(delta_),
       FL = F.fl.u(delta_)) %>% 
  mutate(diffU = c(0,diff(FU)),
         diffL = c(0,diff(FL))) %>% 
  summarise(upper = sum(delta*diffU),
            lower = sum(delta*diffL)) %>% 
  mutate(method = "Frandsen & Lefgrens (2021)" ) %>% 
  mutate(size = 2) 

attflX <- tibble(delta = delta_,
       FU = F.flX.l(delta_),
       FL = F.flX.u(delta_)) %>% 
  mutate(diffU = c(0,diff(FU)),
         diffL = c(0,diff(FL))) %>% 
  summarise(upper = sum(delta*diffU),
            lower = sum(delta*diffL)) %>% 
  mutate(method = "Frandsen & Lefgrens (2021) Covariates" ) %>% 
  mutate(size = 2) 

attwd <- tibble(delta = delta_,
       FU = F.wd.l(delta_),
       FL = F.wd.u(delta_)) %>% 
  mutate(diffU = c(0,diff(FU)),
         diffL = c(0,diff(FL))) %>% 
  summarise(upper = sum(delta*diffU),
            lower = sum(delta*diffL)) %>% 
  mutate(method = "Worst-Case\n[Williamson-Downs (1990)]" ) %>% 
  mutate(size =1)

att_manski_wc <- 
  df_bounds_manski %>% 
  filter(assumption == "WC") %>% 
  select(upper = UB, 
         lower = LB) %>% 
  mutate(method = "Worst-Case\n[Manski]")

att_manski_mtr <- 
  df_bounds_manski %>% 
  filter(assumption == "MTR") %>% 
  select(upper = UB, 
         lower = LB) %>% 
  mutate(method = "Monotone Treatment Response\n[Manski]")

att_manski_mts <- 
  df_bounds_manski %>% 
  filter(assumption == "MTS") %>% 
  select(upper = UB, 
         lower = LB) %>% 
  mutate(method = "Monotone Treatment Selection\n[Manski]")

p <- 
  attfl %>% 
  bind_rows(attflX) %>% 
  bind_rows(attwd) %>% 
  bind_rows(att_manski_wc) %>% 
  bind_rows(att_manski_mtr) %>% 
  bind_rows(att_manski_mts) %>% 
  mutate(y = 0) %>% 
  mutate(method = factor(method, levels = c("Frandsen & Lefgrens (2021)","Frandsen & Lefgrens (2021) Covariates", "Worst-Case\n[Williamson-Downs (1990)]",
      "Worst-Case\n[Manski]", "Monotone Treatment Response\n[Manski]",
      "Monotone Treatment Selection\n[Manski]")))  %>% 
  ggplot() + 
  geom_errorbar(aes(y = y, x = lower, xmin = lower, xmax = upper,colour = method),width = 0.8,lwd = 2) + 
  facet_grid(method~.,switch="y") +
  scale_y_continuous(limits = c(-1,1), breaks = NULL) + 
  geom_point(aes(x = 2, y = 0), colour = "darkred", size=8, alpha = 0.5) + 
  geom_vline(aes(xintercept = 0), lty=3, col = "darkgrey") + 
  labs(x = "ATT", y = "") +
  scale_colour_aaas() +
  theme(legend.position = "none",strip.text.y.left = element_text(angle=0,hjust=1),
        strip.placement = "outside")
p


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
 





## knitr::purl(input = here::here("blog/drafts/partial-id-intro/partial-id-intro.qmd"),
##             output = here::here("blog/drafts/partial-id-intro/partial-id-intro.r"))

