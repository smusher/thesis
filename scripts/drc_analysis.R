library(tidyverse)
library(broom)
library(nls.multstart)

concresp <-
    read_csv("data/combined_drc_data.csv") %>%
    mutate(log_concentration = log10(concentration)) %>%
    filter(log_concentration > -Inf) %>%
    group_by(unique_experiment_id) %>%
    mutate(observations = n()) %>%
    filter(observations >= 3) %>%
    ungroup() %>%
    drop_na()

currents <-
    concresp %>%
    filter(measure == "current")

fluorescence <-
    concresp %>%
    filter(measure == "fluorescence")

hill_equation <- function(log_concentration, ec50, hill, floor){
    floor + ((1 - floor) / (1 + 10^((ec50 - log_concentration) * -hill)))
}

#population fits
unadjusted_fits_free_floor_pop <-
    concresp %>%
    group_by(construct, method, measure, nucleotide) %>%
    nest() %>%
    mutate(fit = map(data, ~ nls_multstart(
            response ~ hill_equation(log_concentration, ec50, hill, floor),
            data = data.frame(.),
            iter = c(10,10,10),
            lower = c(-Inf, 0, 0),
            start_lower = c(ec50 = -7, hill = 0, floor= 0),
            start_upper = c(ec50 = -1, hill = 2, floor = 1),
            na.action = na.omit
              )
    )
    )

adjusted_fits_free_floor_pop <-
    fluorescence %>%
    group_by(construct, method, measure, nucleotide) %>%
    nest() %>%
    mutate(fit = map(data, ~ nls_multstart(
            log2(response + 1) ~ hill_equation(log_concentration, ec50, hill, floor),
            data = data.frame(.),
            iter = c(10, 10, 10),
            lower = c(-Inf, 0, 0),
            start_lower = c(ec50 = -7, hill = 0, floor= 0),
            start_upper = c(ec50 = -1, hill = 2, floor = 1),
            na.action = na.omit
              )
    )
    )

adjusted_fits_fixed_floor_pop <-
    fluorescence %>%
    group_by(construct, method, measure, nucleotide) %>%
    nest() %>%
    mutate(fit = map(data, ~ nls_multstart(
            log2(response + 1) ~ 0.1 + ((1 - 0.1) / (1 + 10^((ec50 - log_concentration) * -hill))),
            data = data.frame(.),
            iter = c(15, 15),
            lower = c(-Inf, 0),
            start_lower = c(ec50 = -7, hill = 0),
            start_upper = c(ec50 = -2, hill = 2),
            na.action = na.omit
              )
    )
    )

saveRDS(unadjusted_fits_free_floor_pop, "data/hill_fits/ua_free_fit_pop.rds")
saveRDS(adjusted_fits_free_floor_pop, "data/hill_fits/a_free_fit_pop.rds")
saveRDS(adjusted_fits_fixed_floor_pop, "data/hill_fits/a_fixed_fit_pop.rds")

unadjusted_fits_free_floor_pop %>%
mutate(tidied = map(fit, tidy)) %>%
unnest(tidied) %>%
select(construct, method, measure, nucleotide, term, estimate) %>%
pivot_wider(names_from = term, values_from = estimate) %>%
rename(pop_ec50 = ec50) %>%
right_join(concresp) -> ua_free_tidy_pop

adjusted_fits_free_floor_pop %>%
mutate(tidied = map(fit, tidy)) %>%
unnest(tidied) %>%
select(construct, method, measure, nucleotide, term, estimate) %>%
pivot_wider(names_from = term, values_from = estimate) %>%
rename(pop_ec50 = ec50) %>%
right_join(fluorescence) -> a_free_tidy_pop

adjusted_fits_fixed_floor_pop %>%
mutate(tidied = map(fit, tidy)) %>%
unnest(tidied) %>%
select(construct, method, measure, nucleotide, term, estimate) %>%
pivot_wider(names_from = term, values_from = estimate) %>%
rename(pop_ec50 = ec50) %>%
mutate(floor = 0.1) %>%
right_join(fluorescence) -> a_fixed_tidy_pop

#individual fits
unadjusted_fits_free_floor <-
    ua_free_tidy_pop %>%
    group_by(construct, method, measure, nucleotide, unique_experiment_id) %>%
    nest() %>%
    mutate(fit = map(data, ~ nls_multstart(
            response ~ hill_equation(log_concentration, ec50, hill, floor),
            data = data.frame(.),
            iter = c(100),
            lower = c(-Inf),
            start_lower = c(ec50 = -7),
            start_upper = c(ec50 = -1),
            na.action = na.omit
              )
    )
    )

adjusted_fits_free_floor <-
    a_free_tidy_pop %>%
    group_by(construct, method, measure, nucleotide, unique_experiment_id) %>%
    nest() %>%
    mutate(fit = map(data, ~ nls_multstart(
            log2(response + 1) ~ hill_equation(log_concentration, ec50, hill, floor),
            data = data.frame(.),
            iter = c(100),
            lower = c(-Inf),
            start_lower = c(ec50 = -7),
            start_upper = c(ec50 = -1),
            na.action = na.omit
              )
    )
    )

adjusted_fits_fixed_floor <-
    a_fixed_tidy_pop %>%
    group_by(construct, method, measure, nucleotide, unique_experiment_id) %>%
    nest() %>%
    mutate(fit = map(data, ~ nls_multstart(
            log2(response + 1) ~ hill_equation(log_concentration, ec50, hill, floor),
            data = data.frame(.),
            iter = c(100),
            lower = c(-Inf),
            start_lower = c(ec50 = -7),
            start_upper = c(ec50 = -1),
            na.action = na.omit
              )
    )
    )

saveRDS(unadjusted_fits_free_floor, "data/hill_fits/ua_free_fit.rds")
saveRDS(adjusted_fits_free_floor, "data/hill_fits/a_free_fit.rds")
saveRDS(adjusted_fits_fixed_floor, "data/hill_fits/a_fixed_fit.rds")
