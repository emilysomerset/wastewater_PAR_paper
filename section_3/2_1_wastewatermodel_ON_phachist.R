rm(list=ls())
library(dplyr)
library(splines)
# library(OSplines)
library(TMB)
library(aghq)
library(bayesplot)
library(lemon)
library(openxlsx)
library(lubridate)
library(magrittr)
library(reshape2)

load("~/Wastewater/reproduction_number/work_d.RData")
load("~/Wastewater/reproduction_number/work_d_toronto_2024_08_08.RData")

compile(file="~/Wastewater/reproduction_number/cpp/model_ospline_fixedeffects_daily_singleCOV_AR2_transformpaper_censored.cpp")
try(dyn.unload(dynlib("~/Wastewater/reproduction_number/cpp/model_ospline_fixedeffects_daily_singleCOV_AR2_transformpaper_censored")),silent = TRUE)
dyn.load(dynlib("~/Wastewater/reproduction_number/cpp/model_ospline_fixedeffects_daily_singleCOV_AR2_transformpaper_censored"))


compute_weights_precision <- function(x){
  d <- diff(x)
  Precweights <- diag(d)
  Precweights
}
local_poly_helper <- function (knots, refined_x, p = 2) {
  if (min(knots) >= 0) {
    dif <- diff(knots)
    nn <- length(refined_x)
    n <- length(knots)
    D <- matrix(0, nrow = nn, ncol = n - 1)
    for (j in 1:nn) {
      for (i in 1:(n - 1)) {
        if (refined_x[j] <= knots[i]) {
          D[j, i] <- 0
        }
        else if (refined_x[j] <= knots[i + 1] & refined_x[j] >= 
                 knots[i]) {
          D[j, i] <- (1/factorial(p)) * (refined_x[j] - 
                                           knots[i])^p
        }
        else {
          k <- 1:p
          D[j, i] <- sum((dif[i]^k) * ((refined_x[j] - 
                                          knots[i + 1])^(p - k))/(factorial(k) * factorial(p - 
                                                                                             k)))
        }
      }
    }
  }
  else if (max(knots) <= 0) {
    refined_x_neg <- refined_x
    refined_x_neg <- ifelse(refined_x < 0, -refined_x, 0)
    knots_neg <- knots
    knots_neg <- unique(sort(ifelse(knots < 0, -knots, 0)))
    dif <- diff(knots_neg)
    nn <- length(refined_x_neg)
    n <- length(knots_neg)
    D <- matrix(0, nrow = nn, ncol = n - 1)
    for (j in 1:nn) {
      for (i in 1:(n - 1)) {
        if (refined_x_neg[j] <= knots_neg[i]) {
          D[j, i] <- 0
        }
        else if (refined_x_neg[j] <= knots_neg[i + 1] & 
                 refined_x_neg[j] >= knots_neg[i]) {
          D[j, i] <- (1/factorial(p)) * (refined_x_neg[j] - 
                                           knots_neg[i])^p
        }
        else {
          k <- 1:p
          D[j, i] <- sum((dif[i]^k) * ((refined_x_neg[j] - 
                                          knots_neg[i + 1])^(p - k))/(factorial(k) * 
                                                                        factorial(p - k)))
        }
      }
    }
  }
  else {
    refined_x_neg <- refined_x
    refined_x_neg <- ifelse(refined_x < 0, -refined_x, 0)
    knots_neg <- knots
    knots_neg <- unique(sort(ifelse(knots < 0, -knots, 0)))
    dif <- diff(knots_neg)
    nn <- length(refined_x_neg)
    n <- length(knots_neg)
    D1 <- matrix(0, nrow = nn, ncol = n - 1)
    for (j in 1:nn) {
      for (i in 1:(n - 1)) {
        if (refined_x_neg[j] <= knots_neg[i]) {
          D1[j, i] <- 0
        }
        else if (refined_x_neg[j] <= knots_neg[i + 1] & 
                 refined_x_neg[j] >= knots_neg[i]) {
          D1[j, i] <- (1/factorial(p)) * (refined_x_neg[j] - 
                                            knots_neg[i])^p
        }
        else {
          k <- 1:p
          D1[j, i] <- sum((dif[i]^k) * ((refined_x_neg[j] - 
                                           knots_neg[i + 1])^(p - k))/(factorial(k) * 
                                                                         factorial(p - k)))
        }
      }
    }
    refined_x_pos <- refined_x
    refined_x_pos <- ifelse(refined_x > 0, refined_x, 0)
    knots_pos <- knots
    knots_pos <- unique(sort(ifelse(knots > 0, knots, 0)))
    dif <- diff(knots_pos)
    nn <- length(refined_x_pos)
    n <- length(knots_pos)
    D2 <- matrix(0, nrow = nn, ncol = n - 1)
    for (j in 1:nn) {
      for (i in 1:(n - 1)) {
        if (refined_x_pos[j] <= knots_pos[i]) {
          D2[j, i] <- 0
        }
        else if (refined_x_pos[j] <= knots_pos[i + 1] & 
                 refined_x_pos[j] >= knots_pos[i]) {
          D2[j, i] <- (1/factorial(p)) * (refined_x_pos[j] - 
                                            knots_pos[i])^p
        }
        else {
          k <- 1:p
          D2[j, i] <- sum((dif[i]^k) * ((refined_x_pos[j] - 
                                           knots_pos[i + 1])^(p - k))/(factorial(k) * 
                                                                         factorial(p - k)))
        }
      }
    }
    D <- cbind(D1, D2)
  }
  D
}
global_poly_helper <- function (x, p = 2) {
  result <- NULL
  for (i in 1:p) {
    result <- cbind(result, x^(i - 1))
  }
  result
}
prior_conversion <- function (d, prior, p) {
  Cp <- (d^((2 * p) - 1))/(((2 * p) - 1) * (factorial(p - 1)^2))
  prior_q <- list(a = prior$a, u = (prior$u * (1/sqrt(Cp))))
  prior_q
}


work_d <- work_d %>%
  filter(substr(site_id,1,1) %in% "T")



## Run all before. 
source("~/Wastewater/reproduction_number/functions_general/prep_data_original.R")

data_foranalysis <- prep_data(outcome_column_name = "value",
                              denom_column_name = "denom",
                              site_id = "site_id",
                              sample_date = "sample_date", 
                              data = work_d %>% mutate(denom = 1, 
                                                       censored_y = FALSE), 
                              polyOrder = 3,
                              pred_also = FALSE)

data_foranalysis$df_full$site_id %>% unique() #4 stations

ggplot(data_foranalysis$df, aes(sample_date, y))+
  geom_point()+
  facet_wrap(~site_id)+ 
  scale_y_continuous(trans = "log")

tmbdat <- data_foranalysis$tmbdat


prior_IWP <- prior_conversion(d=20, prior =list(u=log(2),alpha = 0.5),p=3)

# Set other priors
tmbdat$u1 = prior_IWP$u
tmbdat$alpha1 = 0.5
tmbdat$u2 = 0.5
tmbdat$alpha2 = 0.5
tmbdat$betaprec = 0.01
tmbdat$lambda_phi = -log(0.5)/50
tmbdat$lambda_tau = -log(0.5)/0.5
tmbdat$lambda_cov = -log(0.5)/0.5



set.seed(2)
init_daily <- rnorm(ncol(tmbdat$daily), 0, 0.1);
init_W <- rnorm(ncol(tmbdat$obs),0,0.1)
tmbparams <- list(
  W = c(rep(0, (ncol(tmbdat$X)+ ncol(tmbdat$B))), init_daily, init_W, 0, diff(init_W)), # W = c(U,beta,Z); 
  theta1 = 10, # -2log(sigma)
  theta2 = 0,
  cov_log = 0,
  theta3 = 0,
  theta4 = 0
)


ff <- TMB::MakeADFun(
  data = tmbdat,
  parameters = tmbparams,
  random = "W",
  DLL = "model_ospline_fixedeffects_daily_singleCOV_AR2_transformpaper_censored",
  silent = TRUE
)

aghq_k = 3

mdl1 <- aghq::marginal_laplace_tmb(ff,k=aghq_k,startingvalue = c(10,0,0,0,0))
samps1 <- aghq::sample_marginal(mdl1, M = 3000) 
marginals <- mdl1$marginals  

df_full =data_foranalysis$df_full

save(file="~/Wastewater/reproduction_number/wastewater_model_ON_phachist.RData", list = c("samps1","marginals","tmbdat","df_full"))


load("~/Wastewater/reproduction_number/wastewater_model_ON_phachist.RData")
source("~/Wastewater/reproduction_number/functions_general/process_results_posteriorsamps.R")

results <- process_results(df_full =df_full,
                           tmbdat = tmbdat, 
                           samps1 = samps1, 
                           polyOrder=3, 
                           AR=TRUE)


save(file="~/Wastewater/reproduction_number/results_model_ON_phachist.RData", list = "results")

# load("~/Wastewater/reproduction_number/wastewater_model_ON_phachist.RData")
# source("~/Wastewater/reproduction_number/functions_general/process_results_original.R")
# results <- process_results(df_full =df_full,
#                            tmbdat = tmbdat, 
#                            samps1 = samps1, 
#                            polyOrder=3,
#                            id_group = 0)
# library(ggplot2)
# ggplot(results$df_full, aes(sample_date, exp_v))+
#   # geom_point(aes(sample_date, y), size = 0.1)+
#   geom_line()+
#   theme_bw()+
#   facet_wrap(~site_id)
# 
# ggplot(results$df_full, aes(sample_date, exp_v_u_fixed))+
#   geom_point(aes(sample_date, y), size = 0.1)+
#     geom_line()+
#   geom_ribbon(aes(ymin = exp_v_u_fixed_lwr, ymax = exp_v_u_fixed_upr), alpha = 0.1)+
#   theme_bw()+
#   facet_wrap(~site_id)
# 
# ggplot(results$df_full, aes(sample_date, v_u))+
#   # geom_point(aes(sample_date, y), size = 0.1)+
#   geom_line()+
#   geom_line(aes(sample_date, v), col = "red")+
#   geom_ribbon(aes(ymin = exp_v_lwr, ymax = v_upr), alpha = 0.1, col = "red")+
#   geom_ribbon(aes(ymin = exp_v_u_lwr, ymax = exp_v_u_upr), alpha = 0.1)+
#   theme_bw()+
#   facet_wrap(~site_id)
# 
# 
# ggplot(results$df_full %>% mutate(resid1 = y - exp_v_u_fixed,
#                                   resid2 = y - exp_v), aes(sample_date, resid1))+
#   geom_point(size = 0.2)+
#   theme_bw()+
#   facet_wrap(~site_id)
# 
# 
# 
# ggplot(results$df_full %>% mutate(resid1 = log(y) - v_u_fixed,
#                                   resid2 = log(y) - v_fixed) %>% filter(!is.na(y)), 
#   aes(sample_date, resid1))+
#   geom_line(size = 0.2)+
#   theme_bw()+
#   facet_wrap(~site_id)
# 
# ggplot(results$df_full %>% mutate(resid1 = log(y) - v_u_fixed,
#                                   resid2 = log(y) - v_fixed) %>% filter(!is.na(y)), 
#        aes(sample_date, resid2))+
#   geom_line(size = 0.2)+
#   theme_bw()+
#   facet_wrap(~site_id)
# 
# results$df_full %>% mutate(resid1 = log(y) - v_u_fixed,
#                            resid2 = log(y) - v_fixed) %>% 
#   filter(site_id=="THU") %>% 
#   filter(!is.na(y))%$% resid2 %>% acf()
# 
# 
# 
# 
# 
# 
# 
# load("~/Wastewater/reproduction_number/wastewater_model_ON_phachist_sepAR.RData")
# source("~/Wastewater/reproduction_number/functions_general/process_results_original_nodaily.R")
# results <- process_results(df_full =df_full,
#                            tmbdat = tmbdat, 
#                            samps1 = samps1, 
#                            polyOrder=3,
#                            id_group = 0)
# library(ggplot2)
# ggplot(results$df_full, aes(sample_date, exp_v_fixed))+
#   geom_point(aes(sample_date, y), size = 0.1)+
#   geom_line()+
#   theme_bw()+
#   facet_wrap(~site_id)
# 
# ggplot(results$df_full, aes(sample_date, exp_v_u_fixed))+
#   geom_point(aes(sample_date, y), size = 0.1)+
#   geom_line()+
#   geom_ribbon(aes(ymin = exp_v_u_fixed_lwr, ymax = exp_v_u_fixed_upr), alpha = 0.1)+
#   theme_bw()+
#   facet_wrap(~site_id)
# 
# ggplot(results$df_full, aes(sample_date, exp_v_u))+
#   # geom_point(aes(sample_date, y), size = 0.1)+
#   geom_line()+
#   geom_line(aes(sample_date, exp_v), col = "red")+
#   geom_ribbon(aes(ymin = exp_v_lwr, ymax = exp_v_upr), alpha = 0.1, col = "red")+
#   geom_ribbon(aes(ymin = exp_v_u_lwr, ymax = exp_v_u_upr), alpha = 0.1)+
#   theme_bw()+
#   facet_wrap(~site_id)
# 
# 
# ggplot(results$df_full %>% mutate(resid1 = y - exp_v_u_fixed,
#                                   resid2 = y - exp_v), aes(sample_date, resid1))+
#   geom_point(size = 0.2)+
#   theme_bw()+
#   facet_wrap(~site_id)
# 
# 
# 
# ggplot(results$df_full %>% mutate(resid1 = log(y) - v_u_fixed,
#                                   resid2 = log(y) - v_fixed) %>% filter(!is.na(y)), 
#        aes(sample_date, resid1))+
#   geom_line(size = 0.2)+
#   theme_bw()+
#   facet_wrap(~site_id)
# 
# ggplot(results$df_full %>% mutate(resid1 = log(y) - v_u_fixed,
#                                   resid2 = log(y) - v_fixed) %>% filter(!is.na(y)), 
#        aes(sample_date, resid2))+
#   geom_line(size = 0.2)+
#   theme_bw()+
#   facet_wrap(~site_id)
# 
# results$df_full %>% mutate(resid1 = log(y) - v_u_fixed,
#                            resid2 = log(y) - v_fixed) %>% 
#   filter(site_id=="TAB") %>% 
#   filter(!is.na(y))%$% resid2 %>% acf()
# 
# 
# results$df_full %>% mutate(resid1 = log(y) - v_u_fixed,
#                            resid2 = log(y) - v_fixed) %>% 
#   filter(site_id=="THC") %>% 
#   filter(!is.na(y))%$% resid2 %>% acf()
# 
# 
# 
# 
# 
# load("~/Wastewater/reproduction_number/wastewater_model_ON_phachist_sepAR2.RData")
# source("~/Wastewater/reproduction_number/functions_general/process_results_original.R")
# results <- process_results(df_full =df_full,
#                            tmbdat = tmbdat, 
#                            samps1 = samps1, 
#                            polyOrder=3,
#                            id_group = 0)
# library(ggplot2)
# ggplot(results$df_full, aes(sample_date, exp_v_fixed))+
#   geom_point(aes(sample_date, y), size = 0.1)+
#   geom_line()+
#   theme_bw()+
#   facet_wrap(~site_id)
# 
# ggplot(results$df_full, aes(sample_date, exp_v_u_fixed))+
#   geom_point(aes(sample_date, y), size = 0.1)+
#   geom_line()+
#   geom_ribbon(aes(ymin = exp_v_u_fixed_lwr, ymax = exp_v_u_fixed_upr), alpha = 0.1)+
#   theme_bw()+
#   facet_wrap(~site_id)
# 
# ggplot(results$df_full, aes(sample_date, exp_v_u))+
#   # geom_point(aes(sample_date, y), size = 0.1)+
#   geom_line()+
#   geom_line(aes(sample_date, exp_v), col = "red")+
#   geom_ribbon(aes(ymin = exp_v_lwr, ymax = exp_v_upr), alpha = 0.1, col = "red")+
#   geom_ribbon(aes(ymin = exp_v_u_lwr, ymax = exp_v_u_upr), alpha = 0.1)+
#   theme_bw()+
#   facet_wrap(~site_id)
# 
# 
# ggplot(results$df_full %>% mutate(resid1 = y - exp_v_u_fixed,
#                                   resid2 = y - exp_v), aes(sample_date, resid1))+
#   geom_point(size = 0.2)+
#   theme_bw()+
#   facet_wrap(~site_id)
# 
# 
# 
# ggplot(results$df_full %>% mutate(resid1 = log(y) - v_u_fixed,
#                                   resid2 = log(y) - v_fixed) %>% filter(!is.na(y)), 
#        aes(sample_date, resid1))+
#   geom_line(size = 0.2)+
#   theme_bw()+
#   facet_wrap(~site_id)
# 
# ggplot(results$df_full %>% mutate(resid1 = log(y) - v_u_fixed,
#                                   resid2 = log(y) - v_fixed) %>% filter(!is.na(y)), 
#        aes(sample_date, resid2))+
#   geom_line(size = 0.2)+
#   theme_bw()+
#   facet_wrap(~site_id)
# 
# results$df_full %>% mutate(resid1 = log(y) - v_u_fixed,
#                            resid2 = log(y) - v_fixed) %>% 
#   filter(site_id=="TAB") %>% 
#   filter(!is.na(y))%$% resid2 %>% acf()
# 
# 
# results$df_full %>% mutate(resid1 = log(y) - v_u_fixed,
#                            resid2 = log(y) - v_fixed) %>% 
#   filter(site_id=="THC") %>% 
#   filter(!is.na(y))%$% resid2 %>% acf()
# 
# 
# 
# 
# 
# 
# 
# 
# 
# #####
# load("~/Wastewater/reproduction_number/results_model_ON_phachist.RData")
# tmp <- results$v %>% group_by(sample_date) %>% summarise(med = median(exp(value)),
#                                                          lwr = quantile(exp(value), 0.025),
#                                                          upr = quantile(exp(value), 0.975))
# 
# tmp2 <- results$v_u %>% group_by(sample_date,site_id) %>% summarise(med = median(exp(value)),
#                                                                     lwr = quantile(exp(value), 0.025),
#                                                                     upr = quantile(exp(value), 0.975))
# 
# tmp3 <- results$v_u %>% 
#   mutate(weights = 1,
#          weights = ifelse(site_id == "TAB", 0.499/(0.177+0.232+0.064+0.499),weights),
#          weights = ifelse(site_id == "THC",0.177/(0.177+0.232+0.064+0.499),weights),
#          weights = ifelse(site_id == "THU",0.232/(0.177+0.232+0.064+0.499),weights),
#          weights = ifelse(site_id == "TNT",0.064/(0.177+0.232+0.064+0.499),weights)) %>% 
#   mutate(value = exp(value)) %>% 
#   mutate(value = value*weights) %>% 
#   group_by(variable, sample_date) %>% 
#   summarise(value = sum(value)) %>% 
#   mutate(value = log(value)) %>% 
#   group_by(sample_date) %>% summarise(med = median(exp(value)),
#                                       lwr = quantile(exp(value), 0.025),
#                                       upr = quantile(exp(value), 0.975))
# 
# 
# 
# ggplot(tmp, aes(sample_date, med), col = "red")+ 
#   geom_line(col = "red")+
#   geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.1)+ 
#   geom_line(data = results$df_full, aes(sample_date, exp_v))
# 
# ggplot(tmp2, aes(sample_date, med), col = "red")+ 
#   geom_line(col = "red")+
#   geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.1)+ 
#   geom_line(data = results$df_full, aes(sample_date, exp_v_u))+
#   facet_wrap(~site_id)
# 
# ggplot(tmp, aes(sample_date, med))+ 
#   geom_line(col = "blue")+
#   geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.1, fill = "blue") +
#   geom_line(data = tmp2, aes(sample_date, med),col = "red")+
#   geom_ribbon(data = tmp2,aes(ymin = lwr, ymax = upr), alpha = 0.1, fill = "red") +
#   geom_line(data = tmp3, aes(sample_date, med),col = "green")+
#   geom_ribbon(data = tmp3,aes(ymin = lwr, ymax = upr), alpha = 0.1, fill = "green") +
#   scale_y_continuous(limits = c(0,500))+
#   facet_wrap(~site_id)
# 
# 
# W(t)~ Gamma( )
# 
# log(eta_i(t)) = V(t) + U_i(t)
# 
# # signal: exp(V(t)) -> compute the reproduction
# # signal: weighted sum of exp(V(t)+U_i(t))
# 
# 
# ###################
# 
# 
# load("~/Wastewater/reproduction_number/results_model_ON.RData")
# 
# tmp <- results$v %>% group_by(sample_date) %>% summarise(med = median(exp(value)),
#                                                          lwr = quantile(exp(value), 0.025),
#                                                          upr = quantile(exp(value), 0.975))
# 
# tmp2 <- results$v_u %>% group_by(sample_date,site_id) %>% summarise(med = median(exp(value)),
#                                                                     lwr = quantile(exp(value), 0.025),
#                                                                     upr = quantile(exp(value), 0.975))
# 
# tmp3 <- results$v_u %>% 
#   mutate(weights = 1,
#          weights = ifelse(site_id == "Toronto Ashbridges Bay", 0.499/(0.177+0.232+0.064+0.499),weights),
#          weights = ifelse(site_id == "Toronto Highland Creek",0.177/(0.177+0.232+0.064+0.499),weights),
#          weights = ifelse(site_id == "Toronto Humber",0.232/(0.177+0.232+0.064+0.499),weights),
#          weights = ifelse(site_id == "Toronto North Toronto",0.064/(0.177+0.232+0.064+0.499),weights)) %>% 
#   mutate(value = exp(value)) %>% 
#   mutate(value = value*weights) %>% 
#   group_by(variable, sample_date) %>% 
#   summarise(value = sum(value)) %>% 
#   mutate(value = log(value)) %>% 
#   group_by(sample_date) %>% summarise(med = median(exp(value)),
#                                       lwr = quantile(exp(value), 0.025),
#                                       upr = quantile(exp(value), 0.975))
# 
# ggplot(tmp, aes(sample_date, med))+ 
#   geom_line(col = "red")+
#   geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.1, fill = "red") +
#   geom_line(data = tmp2, aes(sample_date, med),col = "blue")+
#   geom_ribbon(data = tmp2,aes(ymin = lwr, ymax = upr), alpha = 0.1, fill = "blue") +
#   geom_line(data = tmp3, aes(sample_date, med),col = "green")+
#   geom_ribbon(data = tmp3,aes(ymin = lwr, ymax = upr), alpha = 0.1, fill = "green") +
#   facet_wrap(~site_id)
# 
# 
# load("~/Wastewater/reproduction_number/results_model2.RData")
# 
# tmp <- results$v %>% group_by(sample_date) %>% summarise(med = median(exp(value)),
#                                                          lwr = quantile(exp(value), 0.025),
#                                                          upr = quantile(exp(value), 0.975))
# 
# load("~/Wastewater/reproduction_number/results_model_ON_phachist.RData")
# 
# tmp2 <- results$v %>% group_by(sample_date) %>% summarise(med = median(exp(value)),
#                                                          lwr = quantile(exp(value), 0.025),
#                                                          upr = quantile(exp(value), 0.975))
# 
# ggplot(tmp, aes(sample_date, med))+ 
#   geom_line(col = "red")+
#   geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.1, fill = "red") +
#   geom_line(data = tmp2, aes(sample_date, med),col = "blue")+
#   geom_ribbon(data = tmp2,aes(ymin = lwr, ymax = upr), alpha = 0.1, fill = "blue") 
