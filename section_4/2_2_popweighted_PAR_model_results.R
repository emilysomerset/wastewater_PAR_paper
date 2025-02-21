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

load(file="./section_4/results_popweighted/results_processed_NZ_v_u_fixed.RData")


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

ratio_raw = data.frame(date = seq(min(cases$date), max(cases$date),1)) %>% 
  left_join(cases, by = "date") %>% filter(date>=ymd('2022-01-01')-days(7)) %>% 
  mutate(ratio = c(NA, nCasesMovingAverage[2:nrow(.)]/nCasesMovingAverage[1:(nrow(.)-1)])) %>% 
  dplyr::select(date, ratio)

df_plt %>% filter(variable == "(e) Instantaneous reproduction number") %>% 
  left_join(ratio_raw, by = "date") %>% 
  ggplot(aes(date, mean, col = alpha, fill = alpha))+
  geom_line()+ 
  geom_ribbon(aes(ymin= lower, ymax = upper), alpha = 0.2)+ 
  geom_line(data=results_ww$station_ave_df, aes(sample_date, inst_repro), inherit.aes = FALSE)+
  geom_ribbon(data=results_ww$station_ave_df, aes(sample_date, ymin = inst_repro_lwr, ymax = inst_repro_upr), inherit.aes = FALSE, alpha = 0.5)+
  geom_line(aes(date, ratio),inherit.aes = FALSE, col = "red")+ 
  scale_x_date(breaks=scales::pretty_breaks(n=10), 
               date_labels ="%b",
               limits = c(ymd("2022-01-01"),ymd("2023-03-25")),
               name = "",
               sec.axis = sec_axis(name = "",
                                   trans = ~ .,
                                   labels = function(x) {
                                     years <- year(x)
                                     years[duplicated(years)] <- ""  # Remove duplicate year labels
                                     years}))

df_plt_weekly <- df_plt %>% 
  filter(variable == "(a) Infection incidence") %>% 
  mutate(dow = wday(date,label = TRUE),
         dow_num = wday(date)) %>% 
  group_by(alpha) %>% 
  mutate(starting_sun = date[which(dow=='Sun')[1]]) %>% 
  filter(date >= starting_sun) %>% 
  mutate(mean = c(zoo::rollapply(mean,7,sum,na.rm=TRUE),rep(NA,6))) %>%
  filter(!is.na(mean)) %>% 
  filter(date %in% seq(min(date),max(date),7)) %>% 
  dplyr::select(mean, t, date,alpha) %>% 
  rename('earliest_week_start_date' = date) %>% 
  mutate(earliest_week_end_date = earliest_week_start_date + 6) 

####
### Effective reproduction numbers
set.seed(2917)
ss <- sample(1:3000,500)
tmp = Reduce(rbind, lapply(as.list(ss),function(i){data_foranalysis$analysis_d[[i]]$version = i; data_foranalysis$analysis_d[[i]]}))
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

gg1 = tmp %>% 
  ggplot(aes(earliest_week_end_date, ratio_crude))+
  geom_ribbon(aes(earliest_week_end_date,ymax = ratio_v_u_fixed_upr, ymin = ratio_v_u_fixed_lwr), alpha = 0.3)+  
  geom_line(aes(earliest_week_end_date,ratio_v_u_fixed_med))+  
  # geom_point(col="red")+
  geom_line(col ="red", size = 0.5)+
  theme_bw()+ 
  scale_y_continuous(name = expression(R[j]),
                     breaks = scales::pretty_breaks(n=10))+
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

##

gg2=ggplot(results , aes(earliest_week_end_date, z_cumsum_noadj_med))+
  geom_line()+
  geom_line(aes(earliest_week_end_date, number_of_cases_cumsum), linetype = "dashed")+
  geom_ribbon(aes(ymin = z_cumsum_noadj_lwr, ymax = z_cumsum_noadj_upr), alpha = 0.3)+
  theme_bw()+
  scale_y_continuous(name = expression(C[j]/1000000), breaks = scales::pretty_breaks(n=8),
                     labels = function(x) x / 1000000)+
  # # geom_hline(yintercept = 5.15*1000000)+ 
  geom_line(data = df_plt %>% filter(variable == "(b) Cumulative infections"),
            aes(date, mean, group = alpha, col=alpha))+
  geom_ribbon(data = df_plt %>% filter(variable == "(b) Cumulative infections") ,
              aes(date, ymin = lower, ymax = upper, group = alpha, fill=alpha),inherit.aes = FALSE, alpha = 0.2)+
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
        legend.position = c(0, 1),  # Places legend inside the top-right
        legend.justification = c(0, 1),
        legend.direction = "horizontal",
        legend.background = element_blank(),
        legend.key = element_blank())+
  geom_point(data=borderWorkerInfections,aes(date, cumsum(cases)),size=0.2) + 
  labs(col = "",
       fill = "")
  



## 
gg3 = ggplot(results, aes(earliest_week_end_date, z_med))+
  geom_line()+
  # geom_point()+
  geom_line(aes(earliest_week_end_date, y), linetype = "dashed")+
  # geom_point(aes(earliest_week_end_date, y), col = "red")+
  geom_ribbon(aes(ymin = z_lwr, ymax = z_upr), alpha = 0.3)+
  theme_bw()+ 
  scale_y_continuous(name = expression(I[j]/100000), breaks = scales::pretty_breaks(n=8),
                     labels = function(x) x / 100000)+
  geom_line(data = df_plt_weekly,
            aes(earliest_week_end_date, mean, group = alpha, col=alpha))+
  # geom_point(data = df_plt_weekly,
  #            aes(earliest_week_end_date, mean, group = alpha, col=alpha))+
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
  labs(col = '')



gg4=ggplot(results,  aes(earliest_week_end_date, p_med))+
  geom_line()+
  # geom_point()+
  geom_ribbon(aes(ymin = p_lwr, ymax = p_upr), alpha = 0.3)+
  theme_bw()+
  scale_y_continuous(name = expression(pi[j]), 
                     breaks = scales::pretty_breaks(n=10))+
  geom_line(data = df_plt %>% filter(variable == "(c) Absolute case ascertainment rate"),
            aes(date, mean, group = alpha, col=alpha))+
  # geom_ribbon(data = df_plt %>% filter(variable == "(c) Absolute case ascertainment rate"),
  #             aes(date, ymin = lower, ymax = upper, group = alpha, fill=alpha),inherit.aes = FALSE, alpha = 0.4)+
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
  labs(col = '')

fest1 = cowplot::plot_grid(add_sub(gg1+ theme(axis.title.y=element_text(vjust=0)), 
                                   "a) Reproduction numbers"),
                           add_sub(gg3, 
                                   "b) Infection counts"),
                           add_sub(gg4+ theme(axis.title.y=element_text(vjust=0)), 
                                   "c) Ascertainment probability"),
                           add_sub(gg2+ theme(axis.title.y=element_text(vjust=0)) , 
                                   "d) Cumulative infections"),
                           align = "v")

ggsave(filename = paste0("./section_4/plots/epidemic_model_NZ.pdf"),
       plot = grid.arrange(fest1), 
       device = "pdf",
       width = 8, 
       height = 8/3*2,
       dpi = 300)

rstudioapi::viewer(paste0("./section_4/plots/epidemic_model_NZ.pdf"))







####

df_plt%>% filter(variable == "(b) Cumulative infections") %>% 
  dplyr::select(lower, upper, date,alpha) %>% 
  filter(alpha == 2*10^9) %>% 
  left_join(borderWorkerInfections, by = "date") %>% 
  filter(!is.na(cases)) %>% 
  mutate(cases = cumsum(cases)) %>% 
  mutate(in_interval = cases >= lower & cases <= upper)%$% in_interval %>% table()

df_plt%>% filter(variable == "(b) Cumulative infections") %>% 
  dplyr::select(lower, upper, date,alpha) %>% 
  filter(alpha == 3*10^9) %>% 
  left_join(borderWorkerInfections, by = "date") %>% 
  filter(!is.na(cases)) %>% 
  mutate(cases = cumsum(cases)) %>% 
  mutate(in_interval = cases >= lower & cases <= upper)%$% in_interval %>% table()

df_plt%>% filter(variable == "(b) Cumulative infections") %>% 
  dplyr::select(lower, upper, date,alpha) %>% 
  filter(alpha == 4*10^9) %>% 
  left_join(borderWorkerInfections, by = "date") %>% 
  filter(!is.na(cases)) %>% 
  mutate(cases = cumsum(cases)) %>% 
  mutate(in_interval = cases >= lower & cases <= upper)%$% in_interval %>% table()


data.frame(date = seq(min(results$earliest_week_end_date), max(results$earliest_week_end_date),1)) %>% 
  left_join(results, by = c("date"="earliest_week_end_date")) %>% 
  dplyr::select(date, z_cumsum_noadj_lwr,  z_cumsum_noadj_upr) %>% 
  mutate_all(zoo::na.locf) %>% 
  left_join(borderWorkerInfections, by = "date") %>% 
  filter(!is.na(cases)) %>% 
  mutate(cases = cumsum(cases)) %>% 
  mutate(in_interval = cases >= z_cumsum_noadj_lwr & cases <= z_cumsum_noadj_upr)%$% in_interval %>% table()

gg2=ggplot(results , aes(earliest_week_end_date, z_cumsum_noadj_med))+
  # geom_line()+
  # geom_line(aes(earliest_week_end_date, number_of_cases_cumsum), col = "red")+
  # geom_ribbon(aes(ymin = z_cumsum_noadj_lwr, ymax = z_cumsum_noadj_upr), alpha = 0.3)+
  theme_bw()+
  scale_y_continuous(name = "Cumulative Infections", breaks = scales::pretty_breaks(n=8),trans = "log")+
  # # geom_hline(yintercept = 5.15*1000000)+ 
  # geom_line(data = df_plt %>% filter(variable == "(b) Cumulative infections"),
  #           aes(date, mean, group = alpha, col=alpha))+
  geom_ribbon(data = df_plt %>% filter(variable == "(b) Cumulative infections") %>% filter(alpha == 4*10^9),
              aes(date, ymin = lower, ymax = upper, group = alpha, fill=alpha),inherit.aes = FALSE, alpha = 0.2)+
  scale_x_date(breaks=scales::pretty_breaks(n=10), 
               date_labels ="%b",
               limits = c(ymd("2022-02-01"),ymd("2022-07-24")),
               name = "",
               sec.axis = sec_axis(name = "",
                                   trans = ~ .,
                                   labels = function(x) {
                                     years <- year(x)
                                     years[duplicated(years)] <- ""  # Remove duplicate year labels
                                     years}))+
  theme(axis.ticks.x.top = element_blank(),
        axis.text.x.top = element_text(vjust = -66))+
  geom_point(data=borderWorkerInfections,aes(date, cumsum(cases)),size=0.2) 
