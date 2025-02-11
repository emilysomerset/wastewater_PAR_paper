rm(list=ls())

# R version 4.4.2 (2024-10-31)
# Platform: x86_64-pc-linux-gnu
# Running under: Ubuntu 22.04.5 LTS

library(dplyr) # dplyr_1.1.4
library(TMB) # TMB_1.9.16
library(aghq) # aghq_0.4.1
library(magrittr) # magrittr_2.0.3
library(reshape2) # reshape2_1.4.4
library(gridExtra) # gridExtra_2.3
library(cowplot) # cowplot_1.1.3
library(readr) #readr_2.1.5
library(ggplot2) #ggplot2_3.5.1
library(lubridate) #lubridate_1.9.4

load("./section_4/wastewater_model_newzealand_censored_comp_fullAR.RData")
source("./functions_general/process_results_original.R")

## Get site weights 
sites <- read_csv("./section_4/data/sites.csv")
sites <- sites %>% 
  dplyr::select(SampleLocation, Region, Population) %>% 
  rename("site_id"= SampleLocation,
         "weights"=Population)

sites2 <- sites %>% 
  group_by(Region) %>% 
  mutate(region_weight = sum(weights)) %>% 
  mutate(region_weight = ifelse(Region == "Auckland", region_weight-weights[which(site_id == "AU_Mangere")], region_weight)) %>% 
  mutate(weights = ifelse(site_id == "AU_Mangere", 1.716*1000000- region_weight, weights))

# If your Rstudio doesn't have enough memory
# Step 1: Open terminal,
# 
# Step 2:
#   
#   cd ~
#   touch .Renviron
# open .Renviron
# Step 3: Save the following as the first line of .Renviron:
#   
#   R_MAX_VSIZE=100Gb 
# Step 4: Close RStudio and reopen

results <- process_results(df_full =df_full,
                           tmbdat = tmbdat,
                           samps1 = samps1,
                           polyOrder=3,
                           id_group = 1,
                           weights_df = sites2 %>% dplyr::select(site_id, weights) )

save(file="./section_4/results_forplots_model_NZ.RData", list = "results")

##############################################
load(file="./section_4/results_forplots_model_NZ.RData")

gg1 <- results$df_full %>% 
  group_by(sample_date) %>% 
  slice(1) %>%
  ungroup() %>% 
  ggplot(aes(x = sample_date, y)) +
  geom_line(aes(y = exp_v),size=0.2) + 
  geom_ribbon(aes(ymax = exp_v_upr, ymin = exp_v_lwr), alpha = 0.4, size = 0.2, fill = "black") +
  theme_bw()+
  scale_y_continuous(name = "exp(V(t))", breaks = scales::pretty_breaks(n=5))+ 
  scale_x_date(breaks=scales::pretty_breaks(n=10), name = "",date_labels ="%b",
               sec.axis = sec_axis(name = "",trans = ~ .,labels = function(x) year(x)))+
  theme(axis.title.y = element_text(size = 12),
        axis.text.x.top = element_text(vjust = -68),
        axis.ticks.x.top = element_blank())


gg2 <- results$df_full %>% 
  group_by(sample_date) %>%
  slice(1) %>% 
  ungroup() %>% 
  ggplot(aes(x = sample_date,y)) + 
  geom_line(aes(y = exp_v_deriv),size=0.2) + 
  geom_ribbon(aes(ymax = exp_v_deriv_upr, ymin = exp_v_deriv_lwr), alpha = 0.4, size = 0.2, fill = "black") +
  theme_bw()+
  scale_y_continuous(name = "exp(V(t))'", breaks = scales::pretty_breaks(n=5))+ 
  geom_hline(yintercept=0, lty="dashed")+
  scale_x_date(breaks=scales::pretty_breaks(n=10), name = "",date_labels ="%b",
               sec.axis = sec_axis(name = "",trans = ~ .,labels = function(x) year(x)))+
  theme(axis.title.y = element_text(size = 12),
        axis.text.x.top = element_text(vjust = -68),
        axis.ticks.x.top = element_blank())


gg3 <- results$station_ave_df %>% 
  ggplot(aes(x = sample_date, ave_exp_v_u_fixed)) +
  geom_line(size=0.2) + 
  geom_ribbon(aes(ymax = ave_exp_v_u_fixed_upr, ymin = ave_exp_v_u_fixed_lwr), alpha = 0.4, size = 0.2, fill = "black") +
  theme_bw()+
  scale_y_continuous(name = expression(paste(bar(mu),"(t)")), breaks = scales::pretty_breaks(n=5))+ 
  scale_x_date(breaks=scales::pretty_breaks(n=10), name = "",date_labels ="%b",
               sec.axis = sec_axis(name = "",trans = ~ .,labels = function(x) year(x)))+
  theme(axis.title.y = element_text(size = 16),
        axis.text.x.top = element_text(vjust = -68),
        axis.ticks.x.top = element_blank())

gg4 <- results$station_ave_df %>% 
  ggplot(aes(x = sample_date,ave_exp_v_u_fixed_deriv)) + 
  geom_line(size=0.2) + 
  geom_ribbon(aes(ymax = ave_exp_v_u_fixed_deriv_upr, ymin = ave_exp_v_u_fixed_deriv_lwr), alpha = 0.4, size = 0.2, fill = "black") +
  theme_bw()+
  scale_y_continuous(name = expression(paste(bar(mu),"(t)'")), breaks = scales::pretty_breaks(n=5))+ 
  geom_hline(yintercept=0, lty="dashed")+
  scale_x_date(breaks=scales::pretty_breaks(n=10), name = "",date_labels ="%b",
               sec.axis = sec_axis(name = "",trans = ~ .,labels = function(x) year(x)))+
  theme(axis.title.y = element_text(size = 16),
        axis.text.x.top = element_text(vjust = -68),
        axis.ticks.x.top = element_blank())


fest1 = cowplot::plot_grid(add_sub(gg1,"a) Common signal",size = 11),
                           add_sub(gg3,"b) Station average signal",size = 11),
                           add_sub(gg2,"c) Derivative of common signal",size = 11),
                           add_sub(gg4,"d) Derivative of station average signal",size = 11),
                           ncol=2, align="v", byrow = TRUE) 



ggsave(filename = paste0("./section_4/plots/allsignals_NZ_wastewatermodel.pdf"),
       plot = grid.arrange(fest1), 
       device = "pdf",
       width = 8, 
       height = 16/3,
       dpi = 300)

rstudioapi::viewer(paste0("./section_4/plots/allsignals_NZ_wastewatermodel.pdf"))


#### Do the plot with common + station-specific

###############################################
# Posteriors and prior plots
load("./section_4/wastewater_model_newzealand_censored_comp_fullAR.RData")
### 20-day prediction SD 

theta_logprior <- function(theta,prior_alpha = tmbdat$alpha1,prior_u = c(tmbdat$u1)) {
  lambda <- -log(prior_alpha)/prior_u
  log(lambda/2) - lambda * exp(-theta/2) - theta/2
}
priorfunc <- function(x) exp(theta_logprior(x))
priorfuncsigma <- function(x) (2/x) * exp(theta_logprior(-2*log(x)))


### Inference for the smoothing parameter:
prec_marg <- marginals[[1]]
logpostsigma <- compute_pdf_and_cdf(prec_marg,list(totheta = function(x) -2*log(x),fromtheta = function(x) exp(-x/2)),interpolation = 'spline')
smooth_var <- data.frame(SD = logpostsigma$transparam, 
                         density = logpostsigma$pdf_transparam,
                         prior = priorfuncsigma(logpostsigma$transparam))
smooth_var <- rbind(smooth_var, data.frame(SD = 0, density = 0, prior = priorfuncsigma(.Machine$double.eps)))

smooth_var <- data.frame(SD = logpostsigma$transparam, 
                         density = logpostsigma$pdf_transparam)

orig <- log(2)
c <- orig/tmbdat$u1
smooth_var_prior <- data.frame(SD = seq(0, 3/400,0.0001)) %>% 
  mutate(prior = priorfuncsigma(SD),
         prior2 = dexp(c*SD, -log(0.5)/orig))





posterior_sigma_s <- smooth_var_prior %>% 
  ggplot(aes(x = c*SD, y = prior/c)) + 
  geom_area(fill = "orange", alpha = 0.2) +
  theme_classic(base_size = 15) +
  ylab("Density") + 
  scale_x_continuous(name = "",limits = c(0,2), breaks=scales::pretty_breaks(n=6) )+
  geom_line(data= smooth_var, aes(x = c*SD, y = density/c)) 
# +geom_line(aes(x = c*SD, y = prior2)) 

posterior_sigma_s

### Daily RE SD

theta_logprior <- function(theta,prior_alpha = tmbdat$alpha2,prior_u = tmbdat$u2) {
  lambda <- -log(prior_alpha)/prior_u
  log(lambda/2) - lambda * exp(-theta/2) - theta/2
}
priorfunc <- function(x) exp(theta_logprior(x))
priorfuncsigma <- function(x) (2/x) * exp(theta_logprior(-2*log(x)))


prec_marg <- marginals[[2]]
logpostsigma <- compute_pdf_and_cdf(prec_marg,list(totheta = function(x) -2*log(x),fromtheta = function(x) exp(-x/2)),interpolation = 'spline')
smooth_var <- data.frame(SD = logpostsigma$transparam, 
                         density = logpostsigma$pdf_transparam,
                         prior = priorfuncsigma(logpostsigma$transparam))

smooth_var_prior <- data.frame(SD = c(smooth_var$SD), 
                               prior = c(smooth_var$prior))



smooth_var_prior <- rbind(smooth_var_prior, 
                          data.frame(SD = c(0,1), prior = priorfuncsigma(c(.Machine$double.eps,1))))


posterior_sigma_z <- smooth_var_prior %>% ggplot(aes(x = SD, y = prior)) + 
  geom_area(fill = "orange", alpha = 0.2) +
  theme_classic(base_size = 15) +
  geom_line(aes(x = SD, y = density), linetype="solid", data = smooth_var) +
  ylab("Density") +
  xlab("") + theme(legend.position = "none")+
  scale_x_continuous(name = "",limits = c(0,1), breaks=scales::pretty_breaks(n=6) )


posterior_sigma_z



theta_logprior <- function(theta,rate = tmbdat$lambda_cov) {
  lambda <- rate
  log(lambda) - lambda * exp(theta) + theta
}
priorfunc <- function(x) exp(theta_logprior(x))
priorfuncsigma <- function(x) abs(1/x) * exp(theta_logprior(log(x)))


prec_marg <- marginals[[3]]
logpostsigma <- compute_pdf_and_cdf(prec_marg,list(totheta = function(x) log(x),fromtheta = function(x) exp(x)),interpolation = 'spline')
smooth_var <- data.frame(SD = logpostsigma$transparam, 
                         density = logpostsigma$pdf_transparam)

smooth_var_prior <- data.frame(SD = seq(0, 5,0.01)) %>% 
  mutate(prior = priorfuncsigma(SD))



posterior_cov <- smooth_var_prior %>% ggplot(aes(x = SD, y = prior)) + 
  geom_area(fill = "orange", alpha = 0.2) +
  theme_classic(base_size = 15) +
  geom_line(aes(x = SD, y = density), linetype="solid", data = smooth_var) +
  ylab("Density") +
  xlab("") + theme(legend.position = "none")+
  scale_x_continuous(name = "",limits = c(0,1), breaks=scales::pretty_breaks(n=6) )

posterior_cov


theta_logprior <- function(theta,rate = tmbdat$lambda_phi) {
  lambda <- rate
  log(lambda) - lambda * exp(theta) + theta
}
priorfunc <- function(x) exp(theta_logprior(x))
priorfuncsigma <- function(x) abs(1/x) * exp(theta_logprior(log(x)))

smooth_var <- data.frame(tau = seq(0, 100, 0.01)) %>% 
  mutate(prior = priorfuncsigma(tau))

range_posterior_samps = exp((samps1$thetasamples[[4]] - samps1$thetasamples[[5]]/2)/sqrt(2))
drange = density(range_posterior_samps, adjust = 2.5)


posterior_range <- smooth_var %>% ggplot(aes(x = tau, y = prior)) + 
  geom_area(fill = "orange", alpha = 0.2) +
  theme_classic(base_size = 15) +
  geom_line(aes(x = x, y = y), linetype="solid", data = data.frame(x = drange$x, y = drange$y)) +
  ylab("Density")  + theme(legend.position = "none")+
  scale_x_continuous(name = "",limits = c(0,40), breaks=scales::pretty_breaks(n=6) )

posterior_range


theta_logprior <- function(theta,rate = tmbdat$lambda_tau) {
  lambda <- rate
  log(lambda) - lambda * exp(theta) + theta
}
priorfunc <- function(x) exp(theta_logprior(x))
priorfuncsigma <- function(x) abs(1/x) * exp(theta_logprior(log(x)))

smooth_var <- data.frame(tau = seq(0, 2, 0.01)) %>% 
  mutate(prior = priorfuncsigma(tau))

range_posterior_samps = exp((samps1$thetasamples[[4]] + samps1$thetasamples[[5]]/2)/sqrt(2))
drange = density(range_posterior_samps, adjust = 2.5)


posterior_tau <- smooth_var %>% ggplot(aes(x = tau, y = prior)) + 
  geom_area(fill = "orange", alpha = 0.2) +
  theme_classic(base_size = 15) +
  geom_line(aes(x = x, y = y), linetype="solid", data = data.frame(x = drange$x, y = drange$y)) +
  ylab("Density") +
  theme(legend.position = "none")+
  scale_x_continuous(name = "",limits = c(0,2), breaks=scales::pretty_breaks(n=6) )

posterior_tau


gg <- grid.arrange(posterior_sigma_s + ggtitle(expression(sigma[v](20))), 
                   posterior_sigma_z + ggtitle(expression(sigma[z])), 
                   posterior_cov + ggtitle(expression(kappa)), 
                   posterior_range + ggtitle(expression(phi[u])), 
                   posterior_tau+ ggtitle(expression(sigma[u])), nrow = 2)



ggsave(filename = "./section_4/plots/hyperparameters_NZ_wastewatermodel.pdf",
       plot = gg, 
       device = "pdf",
       width = 11, 
       height = 7,
       dpi = 300)

rstudioapi::viewer(paste0("./section_4/plots/hyperparameters_NZ_wastewatermodel.pdf"))

~/Library/CloudStorage/OneDrive-UniversityofToronto/Research/wastewater_PAR_paper_code
