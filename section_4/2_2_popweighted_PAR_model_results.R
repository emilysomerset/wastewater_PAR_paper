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

# df_border = read.csv(paste0(dir,"data/borderWorkerInfections.csv")) %>%
#   filter(!is.nan(dailyInfectionsPer100K)) %>%
#   mutate(incidence=dailyInfectionsPer100K * 7 * 51.5,
#          date = as.Date(substr(date, 1, 10), "%d/%m/%Y"),
#          car = NaN) %>%
#   arrange(date) %>%
#   mutate(cumulative=cumsum(incidence))

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

gg4=ggplot(results , aes(earliest_week_start_date, z_cumsum_noadj_med))+
  geom_line()+
  geom_line(aes(earliest_week_start_date, number_of_cases_cumsum), col = "red")+
  geom_ribbon(aes(ymin = z_cumsum_noadj_lwr, ymax = z_cumsum_noadj_upr), alpha = 0.3)+
  theme_bw()+
  scale_y_continuous(name = "Cumulative Infections", breaks = scales::pretty_breaks(n=10))+
  # # geom_hline(yintercept = 5.15*1000000)+ 
  # geom_point(data=borderWorkerInfections,aes(date, cumsum(cases)),size=1) +
  geom_line(data = df_plt %>% filter(variable == "(b) Cumulative infections"),
            aes(date, mean, group = alpha, col=alpha))+
  geom_ribbon(data = df_plt %>% filter(variable == "(b) Cumulative infections") ,
              aes(date, ymin = lower, ymax = upper, group = alpha, fill=alpha),inherit.aes = FALSE, alpha = 0.2)
