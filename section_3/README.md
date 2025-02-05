
<!-- README.md is generated from the README.Rmd. Please edit that file -->

# Ontario epidemic curve reconstruction

This is the code to replicate the results in Section 3 of the
manuscript.

Start by running the script **1_1_wastewatermodel_ON.R** to fit the
wastewater model. For more information on the wastewater model, please
see <https://doi.org/10.1093/jrsssc/qlae073> and
<https://github.com/emilysomerset/wastewater_paper_code>.

For the remainder of the analysis, you will need the following
libraries. The R version used for this analysis and the package versions
are listed below.

``` r
# > sessionInfo()
# R version 4.4.1 (2024-06-14)
# Platform: aarch64-apple-darwin20
# Running under: macOS Monterey 12.4

# R version 4.4.2 (2024-10-31)
# Platform: x86_64-pc-linux-gnu
# Running under: Ubuntu 22.04.5 LTS

library(dplyr) # dplyr_1.1.4
library(TMB) # TMB_1.9.16
library(magrittr) # magrittr_2.0.3
library(reshape2) # reshape2_1.4.4
library(gridExtra) # gridExtra_2.3
library(cowplot) # cowplot_1.1.3
library(lubridate) # lubridate_1.9.4
library(ggplot2) # ggplot2_3.5.1
```

You will need to source the following functions, provided in this
folder.

``` r
source('../functions_general/prep_data_covid_with_fittedwastewater.R')
source('../functions_general/process_results_epidemic.R')
```

Compile and load the following C++ file

``` r
compile(file="../section_3/cpp/model4.cpp")
```

    ## [1] 0

``` r
try(dyn.unload(dynlib("../section_3/cpp/model4")),silent = TRUE)
dyn.load(dynlib("../section_3/cpp/model4"))
```

## Load case data and prep data

``` r
load("../section_3/data/work_d_toronto_2024_08_08.RData") # Public health ontario covid cases
load("../section_3/data/work_d_pho_testing_2024_10_07.RData") # Public health ontario testing
load(file="../section_3/results_model_ON_phachist.RData") # from the wastewater model

work_d_toronto <- work_d_toronto %>% 
  left_join(work_d_testing %>% 
              dplyr::select(week_end_date, total_number_of_tests) %>% 
              mutate(week_end_date = ymd(week_end_date)), 
            by = c("earliest_week_end_date"="week_end_date"))

# bring in hospital data too

load("../section_3/data/work_d_pho_hosp_cases_2024_10_07.RData")
load("../section_3/data/work_d_toronto_hosp_cases_2024_09_05.RData")

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
```

    ## `summarise()` has grouped output by 'earliest_week_end_date'. You can override
    ## using the `.groups` argument.

![](README_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

![](README_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

## Fit model

``` r
# This model is fit using the code in 2_1_popweighted_fixed_pred.R
# We load the model that was fit using the full dataset (j = 126)

load("../section_3/results_pred_popweighted_fixed/results_126.RData")
```

## Process the results

``` r
results <- process_results(data_foranalysis = data_foranalysis_full,
                MI_models = MI_models,
                mean_adj = 2.57/100*pop, 
                sd_adj = 0.47/100*pop)
```

## Plot the results

### Ascertainment probability

    ## # A tibble: 186 × 36
    ##    week_start_date     y public_health_unit    earliest_week_start_date
    ##    <date>          <dbl> <chr>                 <date>                  
    ##  1 2020-10-04       1762 Toronto Public Health 2020-10-04              
    ##  2 2020-10-11       1881 Toronto Public Health 2020-10-11              
    ##  3 2020-10-18       2213 Toronto Public Health 2020-10-18              
    ##  4 2020-10-25       2243 Toronto Public Health 2020-10-25              
    ##  5 2020-11-01       2562 Toronto Public Health 2020-11-01              
    ##  6 2020-11-08       3281 Toronto Public Health 2020-11-08              
    ##  7 2020-11-15       3189 Toronto Public Health 2020-11-15              
    ##  8 2020-11-22       3460 Toronto Public Health 2020-11-22              
    ##  9 2020-11-29       3846 Toronto Public Health 2020-11-29              
    ## 10 2020-12-06       3739 Toronto Public Health 2020-12-06              
    ## # ℹ 176 more rows
    ## # ℹ 32 more variables: earliest_week_end_date <date>,
    ## #   total_number_of_tests <dbl>, admissions <dbl>, deaths <dbl>,
    ## #   population <dbl>, hospitalized_cases <dbl>, ratio_cases <dbl>,
    ## #   index_start <int>, Y0 <dbl>, ratio_v_u <dbl>, ratio_v_u_fixed <dbl>,
    ## #   row_id <int>, p_med <dbl>, p_upr <dbl>, p_lwr <dbl>, z_med <dbl>,
    ## #   z_upr <dbl>, z_lwr <dbl>, y_med <dbl>, y_upr <dbl>, y_lwr <dbl>, …

<img src="README_files/figure-gfm/unnamed-chunk-9-1.png" height="50%" />

### Weekly new infection counts

![](README_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->

### Cumulative infections

    ## # A tibble: 2 × 5
    ##   earliest_week_end_date z_cumsum_noadj_delta_med z_cumsum_noadj_delta_lwr
    ##   <date>                                    <dbl>                    <dbl>
    ## 1 2021-03-27                                 6.85                     5.42
    ## 2 2021-04-03                                 7.48                     5.92
    ## # ℹ 2 more variables: z_cumsum_noadj_delta_upr <dbl>,
    ## #   number_of_cases_cumsum_delta <dbl>

![](README_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->
