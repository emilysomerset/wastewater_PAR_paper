# > sessionInfo()
# R version 4.4.1 (2024-06-14)
# Platform: aarch64-apple-darwin20
# Running under: macOS Monterey 12.4

# R version 4.4.2 (2024-10-31)
# Platform: x86_64-pc-linux-gnu
# Running under: Ubuntu 22.04.5 LTS

rm(list=ls())
library(dplyr) # dplyr_1.1.4
library(TMB) # TMB_1.9.16
library(magrittr) # magrittr_2.0.3
library(reshape2) # reshape2_1.4.4
library(gridExtra) # gridExtra_2.3
library(cowplot) # cowplot_1.1.3
library(lubridate) # lubridate_1.9.4
library(ggplot2) # ggplot2_3.5.1

## Source the following functions
source('./functions_general/prep_data_covid_with_fittedwastewater.R')
source('./functions_general/process_results_epidemic.R')

## Full model results
load("./section_3/results_pred_popweighted_fixed/results_126.RData")


## Load and prep the data for plots
load("./section_3/data/work_d_toronto_2024_08_08.RData") # Public health ontario covid cases
load("./section_3/data/work_d_pho_testing_2024_10_07.RData") # Public health ontario testing
load(file="./section_3/results_model_ON_phachist.RData") # from the wastewater model

work_d_toronto <- work_d_toronto %>% 
  left_join(work_d_testing %>% 
              dplyr::select(week_end_date, total_number_of_tests) %>% 
              mutate(week_end_date = ymd(week_end_date)), 
            by = c("earliest_week_end_date"="week_end_date"))

# bring in hospital data too

load("./section_3/data/work_d_pho_hosp_cases_2024_10_07.RData")
load("./section_3/data/work_d_toronto_hosp_cases_2024_09_05.RData")

work_d_toronto <- work_d_toronto %>% 
  left_join(work_d_hosp2 %>% 
              filter(outcome == "COVID-19 hospital admissions (up to January 20, 2024)") %>% 
              dplyr::select(week_end_date, number) %>% 
              mutate(week_end_date = ymd(week_end_date)) %>% 
              rename("admissions" = number), 
            by = c("earliest_week_end_date"="week_end_date")) %>% 
  left_join(work_d_hosp2 %>% 
              filter(outcome == "COVID-19 deaths") %>% 
              dplyr::select(week_end_date, number,population) %>% 
              mutate(week_end_date = ymd(week_end_date)) %>% 
              rename("deaths" = number),
            by = c("earliest_week_end_date"="week_end_date")) %>% 
  left_join(work_d_hosp %>% 
              dplyr::select(episode_week, hospitalized_cases) %>% 
              mutate(episode_week = ymd(episode_week)), 
            by = c("week_start_date"="episode_week"))

# These weights are from 
# https://www.toronto.ca/wp-content/uploads/2023/10/8e9f-PublicHealthWastewaterSurveillanceTechNotes.pdf
weights_datadic = data.frame(site_id = c("TAB","THC","THU","TNT"), 
                             weight = c(0.499/(0.177+0.232+0.064+0.499),
                                        0.177/(0.177+0.232+0.064+0.499),
                                        0.232/(0.177+0.232+0.064+0.499),
                                        0.064/(0.177+0.232+0.064+0.499)))

data_foranalysis_full <- prep_data(case_data = work_d_toronto,
                                   y_var = "number_of_cases",
                                   results = results,
                                   AR=TRUE,
                                   weight_ratio = TRUE,
                                   weights_datadic = weights_datadic)

pop <- work_d_toronto %>% 
  filter(earliest_week_end_date < mdy("09-30-2020")) %>% slice(1) %$% population[1]


### Effective reproduction numbers
set.seed(2917)
ss <- sample(1:3000,500)
tmp = Reduce(rbind, lapply(as.list(ss),function(i){data_foranalysis_full$analysis_d[[i]]$version = i; data_foranalysis_full$analysis_d[[i]]}))
tmp <- tmp %>% 
  group_by(earliest_week_end_date,y) %>% 
  summarise(ratio_v_u_fixed_lwr = quantile(ratio_v_u_fixed,0.025),
            ratio_v_u_fixed_upr = quantile(ratio_v_u_fixed,0.975),
            ratio_v_u_fixed_med = quantile(ratio_v_u_fixed,0.5),
            ratio_lwr = quantile(ratio,0.025),
            ratio_upr = quantile(ratio,0.975),
            ratio_med = quantile(ratio,0.5))


tmp <- tmp %>%ungroup() %>%  
  mutate(lagged_case = lag(y)) %>% 
  mutate(ratio_crude = y/lagged_case)


### 
which(tmp$ratio_v_u_fixed_med==(tmp$ratio_v_u_fixed_med %>% max()))
which(tmp$ratio_crude==(tmp$ratio_crude %>% max(na.rm = TRUE)))
tmp[63,]

gg1 = tmp %>% 
  ggplot(aes(earliest_week_end_date, ratio_crude))+
  geom_ribbon(aes(earliest_week_end_date,ymax = ratio_v_u_fixed_upr, ymin = ratio_v_u_fixed_lwr), alpha = 0.3)+  
  geom_line(aes(earliest_week_end_date,ratio_v_u_fixed_med))+  
  # geom_point(col="red")+
  geom_line(col ="red", size = 0.5)+
  theme_bw()+ 
  scale_y_continuous(name = expression(R[j]),
                     breaks = scales::pretty_breaks(n=10))+
  scale_x_date(breaks=scales::pretty_breaks(n=10), name = "",date_labels ="%b",
               sec.axis = sec_axis(name = "",trans = ~ .,labels = function(x) year(x)))+
  theme(axis.title.y = element_text(size = 16),
        axis.text.x.top = element_text(vjust = -68),
        axis.ticks.x.top = element_blank())

### Results
results <- process_results(data_foranalysis = data_foranalysis_full,
                           MI_models = MI_models, 
                           adj = TRUE, 
                           mean_adj = 2.57/100*pop, 
                           sd_adj = 0.47/100*pop)

results %>% 
  filter(earliest_week_end_date %in% c(ymd("2021-12-25"),ymd("2021-12-25")-days(7), ymd("2021-12-25")+days(7)))

### Weekly new infections 
gg2 = ggplot(results, aes(earliest_week_end_date, z_med))+
  geom_ribbon(aes(ymin = z_lwr, ymax = z_upr), alpha = 0.3)+
  geom_line()+
  geom_line(aes(earliest_week_end_date, y), col = "red")+
  theme_bw()+ 
  scale_y_continuous(name = expression(I[j]), breaks = scales::pretty_breaks(n=8),
                     labels = scales::comma)+
  scale_x_date(breaks=scales::pretty_breaks(n=10), name = "",date_labels ="%b",
               sec.axis = sec_axis(name = "",trans = ~ .,labels = function(x) year(x)))+
  theme(axis.title.y = element_text(size = 16),
        axis.text.x.top = element_text(vjust = -68),
        axis.ticks.x.top = element_blank())

### Ascertainment probability
gg3 = ggplot(results ,  aes(earliest_week_end_date, p_med))+
  # geom_point()+
  geom_ribbon(aes(ymin = p_lwr, ymax = p_upr), alpha = 0.3)+
  geom_line()+
  theme_bw()+
  scale_y_continuous(name = expression(pi[j]), breaks = scales::pretty_breaks(n=6))+
  scale_x_date(breaks=scales::pretty_breaks(n=10), name = "",date_labels ="%b",
               sec.axis = sec_axis(name = "",trans = ~ .,labels = function(x) year(x)))+
  theme(axis.title.y = element_text(size = 16),
        axis.text.x.top = element_text(vjust = -68),
        axis.ticks.x.top = element_blank())+
  geom_line(aes(earliest_week_end_date,total_number_of_tests/max(total_number_of_tests)),col = "red")

gg4 = ggplot(results %>% filter(earliest_week_end_date<"2021-12-01"), aes(earliest_week_end_date, z_cumsum_noadj_delta_med/pop*100))+
  geom_ribbon(aes(ymin = z_cumsum_noadj_delta_lwr/pop*100, ymax = z_cumsum_noadj_delta_upr/pop*100), alpha = 0.3)+
  geom_line()+
  theme_bw()+ 
  scale_y_continuous(name = "Cumulative Incidence (%)", breaks = scales::pretty_breaks(n=10))+
  geom_errorbar(data= data.frame(earliest_week_end_date=ymd("2021-03-31"),
                                 lwr = 3.74,
                                 upr = 7.27), 
                aes(x=earliest_week_end_date,ymax =upr, ymin = lwr),inherit.aes = FALSE,size=0.5,width = 10)+
  geom_point(data = data.frame(earliest_week_end_date=ymd("2021-03-31"),y= 5.38 ),
             aes(earliest_week_end_date, y),size= 0.5)+
  geom_line(aes(earliest_week_end_date, number_of_cases_cumsum_delta/pop*100), col = "red")+
  scale_x_date(breaks=scales::pretty_breaks(n=8), name = "")+
  theme(axis.title.y = element_text(size = 10),
        axis.ticks.x.top = element_blank())

fest1 = cowplot::plot_grid(add_sub(gg1+ theme(axis.title.y=element_text(vjust=-3)), 
                                              "a) Effective reproduction numbers"),
                   add_sub(gg2, 
                           "b) Weekly new infections"),
                   add_sub(gg3+ theme(axis.title.y=element_text(vjust=-2)), 
                                      "c) Weekly ascertainment probability"),
                   add_sub(gg4+ theme(axis.title.y=element_text(vjust=-1)) , 
                           "d) Pre-Omicron cumulative incidence"),
                   align = "hv")

ggsave(filename = paste0("./section_3/plots/epidemic_model_toronto.pdf"),
       plot = grid.arrange(fest1), 
       device = "pdf",
       width = 8.3, 
       height = 8/3*2,
       dpi = 300)

### Validation plots
results %>% 
  filter(earliest_week_end_date<= ymd("2021-04-03") & earliest_week_end_date >= ymd("2021-03-27")) %>% 
  dplyr::select(earliest_week_end_date, 
                z_cumsum_noadj_delta_med,
                z_cumsum_noadj_delta_lwr,
                z_cumsum_noadj_delta_upr,
                number_of_cases_cumsum_delta) %>% 
  mutate(z_cumsum_noadj_delta_med=z_cumsum_noadj_delta_med/pop*100,
         z_cumsum_noadj_delta_lwr=z_cumsum_noadj_delta_lwr/pop*100,
         z_cumsum_noadj_delta_upr=z_cumsum_noadj_delta_upr/pop*100,
         number_of_cases_cumsum_delta = number_of_cases_cumsum_delta/pop*100)




load("./section_3/results_pred_popweighted_fixed/results_pred.RData")
rect_left = results_pred %>% 
  group_by(model) %>% 
  slice(c(1))%$% earliest_week_end_date

rect_right = results_pred %>% 
  group_by(model) %>% 
  slice(c(5))%$% earliest_week_end_date


rectangles <- data.frame(
  xmin = rect_left,
  xmax = rect_right,
  ymin = -Inf,
  ymax = Inf
)

dd = scales::pretty_breaks(n=10)
round_custom <- function(x) {
  round(x, -floor(log10(x)) * (x >= 100)) + (x < 100) * (round(x, -1) - x)
}
gg2 = results_pred %>%
  ggplot(aes(earliest_week_end_date, y))+
  theme_bw()+
  geom_rect(data=rectangles, aes(xmin=xmin, xmax=xmax, ymin=0, ymax=Inf), 
            fill='gray80', alpha=0.3,inherit.aes = FALSE)+
  geom_point(size = 0.4,shape = 2)+
  geom_point(aes(earliest_week_end_date, pred_y_med), show.legend = FALSE,size = 0.4)+
  geom_errorbar(aes(earliest_week_end_date, ymin = pred_y_95_lwr, ymax = pred_y_95_upr), 
                show.legend = FALSE,
                linewidth = 0.2)+
  scale_y_continuous(name = expression("Confirmed COVID-19 case counts"), 
                     trans="log10",
                     labels = function(x)format(x,digits=2,big.mark = ","),
                     breaks = round_custom(10^(dd(log10(results_pred$pred_y_med),8))))+
  scale_x_date(breaks=scales::pretty_breaks(n=10), 
               name = "",
               date_labels ="%b",
               sec.axis = sec_axis(name = "",
                                   trans = ~ .,
                                   labels = function(x) {
                                   years <- year(x)
                                   years[duplicated(years)] <- ""  # Remove duplicate year labels
                                   years}),
               expand = c(0.01, 0.01))+
  theme(axis.title.y = element_text(size = 10),
        axis.ticks.x.top = element_blank(),
        axis.text.x.top = element_text(vjust = -91),
        panel.grid.minor.x = element_blank())


ggsave(filename = paste0("./section_3/plots/epidemic_model_validation_toronto.pdf"),
       plot =gg2, 
       device = "pdf",
       width = 7, 
       height = 3,
       dpi = 300)

rstudioapi::viewer(paste0("./section_3/plots/epidemic_model_validation_toronto.pdf"))


### Prior and posterior

rm(list=ls())

library(fanetc)
library(splines)
library(cowplot)
library(OSplines)
library(TMB)
library(aghq)
library(bayesplot)
library(lemon)
library(openxlsx)

# Load Model Object
## Model object should have samps1, df, df_full, tmbdat and marginals
load("./section_3/results_pred_popweighted_fixed/results_126.RData")

### ascertainment probability SD

theta_logprior <- function(theta,prior_alpha = 0.5,prior_u = 1) {
  lambda <- -log(prior_alpha)/prior_u
  log(lambda/2) - lambda * exp(-theta/2) - theta/2
}
priorfunc <- function(x) exp(theta_logprior(x))
priorfuncsigma <- function(x) (2/x) * exp(theta_logprior(-2*log(x)))


posterior_samps <- Reduce(c,lapply(MI_models, function(x)x$thetasamples[[1]]))
dsigma = density(posterior_samps,adjust=2.5)
plot(range)

smooth_var <- data.frame(sigma = seq(0, 5, 0.01)) %>% 
  mutate(prior = priorfuncsigma(sigma))


posterior_sigma_pi = smooth_var %>% ggplot(aes(x = sigma, y = prior)) + 
  geom_area(fill = "orange", alpha = 0.2) +
  theme_classic(base_size = 15) +
  geom_line(aes(x = x, y = y), linetype="solid", data = data.frame(x = dsigma$x, y = dsigma$y)) +
  ylab("Density") +
  theme(legend.position = "none")+
  scale_x_continuous(name = "",limits = c(0,5), breaks=scales::pretty_breaks(n=6) )


library(gridExtra)
gg <- grid.arrange(posterior_sigma_pi + ggtitle(expression(sigma[pi])), nrow = 1)

ggsave(filename = "./section_3/plots/hyperparameters_epidemic.pdf",
       plot = gg, 
       device = "pdf",
       width = 5, 
       height = 4,
       dpi = 300)

rstudioapi::viewer(paste0("./section_3/plots/hyperparameters_epidemic.pdf"))
