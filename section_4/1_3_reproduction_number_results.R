rm(list=ls())

# R version 4.4.2 (2024-10-31)
# Platform: x86_64-pc-linux-gnu
# Running under: Ubuntu 22.04.5 LTS

library(dplyr) # dplyr_1.1.4
library(magrittr) # magrittr_2.0.3
library(reshape2) # reshape2_1.4.4
library(gridExtra) # gridExtra_2.3
library(cowplot) # cowplot_1.1.3
library(readr) #readr_2.1.5
library(ggplot2) #ggplot2_3.5.1
library(lubridate)

load("./section_4/results_forplots_model_NZ.RData")
results_ww <- results

genTimeDist <- read_csv("section_4/data/genTimeDist.csv")

####
dir = "./section_4/data/comparison/"

filename = "final"
st_date = as.Date("2022-01-01")
en_date = as.Date("2023-03-31")

# Load model results
df1 = read.csv(paste0(dir, filename, "_alpha2e9", ".csv")) %>%  mutate(alpha = 2e9)
df2 = read.csv(paste0(dir, filename, ".csv")) %>%  mutate(alpha = 3e9)
df3 = read.csv(paste0(dir, filename, "_alpha4e9", ".csv")) %>%  mutate(alpha = 4e9)

borderWorkerInfections <- read_csv(paste0(dir,"/borderWorkerInfections.csv"))
borderWorkerInfections <- borderWorkerInfections %>% 
  filter(!is.na(dailyInfectionsPer100K)) %>% 
  mutate(cases = dailyInfectionsPer100K*(5.15*1000000/100000)*7) %>% 
  mutate(date = date(dmy_hm(date)))

# Print initial CARs (for interest) and then tidy
CAR0_a = df1 %>% filter(variable=="CARt", t==141) %>% select(mean)
CAR0_b = df2 %>% filter(variable=="CARt", t==141) %>% select(mean)
CAR0_c = df3 %>% filter(variable=="CARt", t==141) %>% select(mean)
CAR0_a = CAR0_a[[1]]
CAR0_b = CAR0_b[[1]]
CAR0_c = CAR0_c[[1]]

# Update with relative CARt
df1_relcar = df1 %>% filter(variable=="CARt") %>% mutate(mean=mean/CAR0_a, lower=lower/CAR0_a, upper=upper/CAR0_a, variable="CARtRelative")
df2_relcar = df2 %>% filter(variable=="CARt") %>% mutate(mean=mean/CAR0_b, lower=lower/CAR0_b, upper=upper/CAR0_b, variable="CARtRelative")
df3_relcar = df3 %>% filter(variable=="CARt") %>% mutate(mean=mean/CAR0_c, lower=lower/CAR0_c, upper=upper/CAR0_c, variable="CARtRelative")

df = rbind(df1, df2, df3, df1_relcar, df2_relcar, df3_relcar) %>%
  mutate(date = as.Date(date), alpha=factor(alpha))
# filter(date >= st_date, date <= en_date)
rm(df1, df2, df3, df1_relcar, df2_relcar, df3_relcar)


# Make plot
df_plt = df %>%
  filter(variable %in% c("It", "CI", "Rt", "CARt", "CARtRelative")) %>%
  mutate(popline = ifelse(variable=="CI", 5.15e6, NA),
         Rline = ifelse(variable=="Rt", 1, NA),
         variable = case_when(variable=="It" ~ "(a) Infection incidence",
                              variable=="CI" ~ "(b) Cumulative infections",
                              variable=="CARt" ~ "(c) Absolute case ascertainment rate",
                              variable=="CARtRelative" ~ "(d) Relative case ascertainment rate",
                              variable=="Rt" ~ "(e) Instantaneous reproduction number"),
         variable = factor(variable, levels = c("(a) Infection incidence",
                                                "(b) Cumulative infections",
                                                "(c) Absolute case ascertainment rate",
                                                "(d) Relative case ascertainment rate",
                                                "(e) Instantaneous reproduction number")))

# cases <- read_csv("./section_4/data/cases.csv")
# ratio_raw = data.frame(date = seq(min(cases$date), max(cases$date),1)) %>% 
#   left_join(cases, by = "date") %>% filter(date>=ymd('2022-01-01')-days(7)) %>% 
#   mutate(ratio = c(NA, nCasesMovingAverage[2:nrow(.)]/nCasesMovingAverage[1:(nrow(.)-1)]),
#          ratio2 = c(NA, nCases[2:nrow(.)]/nCases[1:(nrow(.)-1)])) %>% 
#   dplyr::select(date, ratio,ratio2)

gg1 = df_plt %>% filter(variable == "(e) Instantaneous reproduction number") %>% 
  # left_join(ratio_raw, by = "date") %>%
  group_by(alpha) %>% 
  mutate(cumprod = cumprod(mean)) %>% 
  ggplot(aes(date, mean, col = alpha, fill = alpha))+
  geom_line()+ 
  geom_ribbon(aes(ymin= lower, ymax = upper), alpha = 0.2)+ 
  theme_bw()+ 
  scale_y_continuous(name = expression(R(t)),
                     breaks = scales::pretty_breaks(n=10),
                     limits = c(0,3.2))+
  scale_x_date(breaks=scales::pretty_breaks(n=10), 
               date_labels ="%b",
               limits = c(ymd("2022-01-01"),ymd("2023-03-25")),
               name = "",
               sec.axis = sec_axis(name = "",
                                   trans = ~ .,
                                   labels = function(x) {
                                     years <- year(x)
                                     years[duplicated(years)] <- ""  # Remove duplicate year labels
                                     years}))+
  theme(axis.ticks.x.top = element_blank(),
        axis.text.x.top = element_text(vjust = -66),
        legend.position = c(1, 1),  # Places legend inside the top-right
        legend.justification = c(1, 1),
        legend.direction = "horizontal",
        legend.background = element_blank(),
        legend.key = element_blank())+ 
  labs(col = '',fill = '')

gg2 = df_plt %>% filter(variable == "(e) Instantaneous reproduction number") %>% 
  # left_join(ratio_raw, by = "date") %>%
  group_by(alpha) %>% 
  mutate(cumprod = cumprod(mean)) %>% 
  ggplot(aes(date, cumprod, col = alpha, fill = alpha))+
  # geom_line()+
  geom_line(data=results_ww$station_ave_df, aes(sample_date,prod_inst_repro_med), inherit.aes = FALSE)+
  # geom_line(data=results_ww$station_ave_df %>% mutate(cumprod = cumprod(inst_repro_med)), aes(sample_date,cumprod), inherit.aes = FALSE,col = "red")+
  geom_ribbon(data=results_ww$station_ave_df, aes(sample_date, ymin = prod_inst_repro_lwr, ymax = prod_inst_repro_upr), inherit.aes = FALSE, alpha = 0.5)+
  theme_bw()+ 
  scale_y_continuous(name = expression(tilde(I)(t)),
                     breaks = scales::pretty_breaks(n=8))+
  scale_x_date(breaks=scales::pretty_breaks(n=10), 
               date_labels ="%b",
               limits = c(ymd("2022-01-01"),ymd("2023-03-25")),
               name = "",
               sec.axis = sec_axis(name = "",
                                   trans = ~ .,
                                   labels = function(x) {
                                     years <- year(x)
                                     years[duplicated(years)] <- ""  # Remove duplicate year labels
                                     years}))+
  theme(axis.ticks.x.top = element_blank(),
        axis.text.x.top = element_text(vjust = -66))

R = df_plt %>% filter(variable == "(e) Instantaneous reproduction number") %>% filter(alpha == 2e+09)%$% mean
I = c(1,rep(0,99))
weights = genTimeDist$p[2:101]
weights_matrix = (weights%*% t(R)) 

recursive_weighted_values <- numeric(100)

gg3 = df_plt %>% filter(variable == "(e) Instantaneous reproduction number") %>% 
  group_by(alpha) %>% 
  mutate(cumprod = cumprod(mean)) %>% 
  ggplot(aes(date, cumprod, col = alpha, fill = alpha))+
  geom_line()+
  theme_bw()+ 
  scale_y_continuous(name =  expression(tilde(I)(t)/10^9),
                     breaks = scales::pretty_breaks(n=6),
                     labels = function(x) x / 1000000000)+
  scale_x_date(breaks=scales::pretty_breaks(n=10), 
               date_labels ="%b",
               limits = c(ymd("2022-01-01"),ymd("2023-03-25")),
               name = "",
               sec.axis = sec_axis(name = "",
                                   trans = ~ .,
                                   labels = function(x) {
                                     years <- year(x)
                                     years[duplicated(years)] <- ""  # Remove duplicate year labels
                                     years}))+
  theme(axis.ticks.x.top = element_blank(),
        axis.text.x.top = element_text(vjust = -66),
        legend.position = c(1, 1),  # Places legend inside the top-right
        legend.justification = c(1, 1),
        legend.direction = "horizontal",
        legend.background = element_blank(),
        legend.key = element_blank())+ 
  labs(col = '',fill = '')



gg4 = results_ww$station_ave_df %>% 
  ggplot(aes(sample_date, inst_repro_med))+
  geom_line() +
  geom_ribbon(aes(sample_date, ymin = inst_repro_lwr, ymax = inst_repro_upr),alpha = 0.5)+
  theme_bw()+ 
  scale_y_continuous(name = expression(R(t)),
                     breaks = scales::pretty_breaks(n=10),
                     limits = c(0,3.2))+
  scale_x_date(breaks=scales::pretty_breaks(n=10), 
               date_labels ="%b",
               limits = c(ymd("2022-01-01"),ymd("2023-03-25")),
               name = "",
               sec.axis = sec_axis(name = "",
                                   trans = ~ .,
                                   labels = function(x) {
                                     years <- year(x)
                                     years[duplicated(years)] <- ""  # Remove duplicate year labels
                                     years}))+
  theme(axis.ticks.x.top = element_blank(),
        axis.text.x.top = element_text(vjust = -66),
        legend.position = c(1, 1),  # Places legend inside the top-right
        legend.justification = c(1, 1),
        legend.direction = "horizontal",
        legend.background = element_blank(),
        legend.key = element_blank())+ 
  labs(col = '',fill = '')


fest1 = cowplot::plot_grid(add_sub(gg1+ theme(axis.title.y=element_text(vjust=0)), 
                                   "a) Instantaneous reproduction numbers"),
                           add_sub(gg4+ theme(axis.title.y=element_text(vjust=0)), 
                                   "b) Reproduction numbers"),
                           add_sub(gg3, 
                                   "c) Cumulative product of R(t)"),
                           add_sub(gg2+ theme(axis.title.y=element_text(vjust=0)) , 
                                   "d) Cumulative product of R(t)"),
                           align = "v")

ggsave(filename = paste0("./section_4/plots/reproduction_numbers_comparisons.pdf"),
       plot = grid.arrange(fest1), 
       device = "pdf",
       width = 8.3, 
       height = 8/3*2,
       dpi = 300)

rstudioapi::viewer(paste0("./section_4/plots/reproduction_numbers_comparisons.pdf"))
