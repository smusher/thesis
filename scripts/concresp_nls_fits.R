library(tidyverse)
library(broom)
library(nls.multstart)

currents_1 <-
    read_csv("data/electrophys_data.csv") %>%
    select(n, method, measure, construct, nucleotide, concentration, response) %>%
    mutate(experimenter = "Sam")

currents_2 <-
    read_csv("data/natascias_data.csv") %>%
    select(n, method, measure, construct, nucleotide, concentration, response) %>%
    mutate(experimenter = "Natascia")

currents <-
    bind_rows(currents_1, currents_2) %>%
    group_by(experimenter, n, method, measure, construct, nucleotide) %>%
    mutate(unique_experiment_id = cur_group_id())

fluorescence <-
    read_csv("data/unroofed_concresp_data.csv") %>%
    filter(dye < 480) %>%
    select(unique_experiment_id, method, measure, construct, concentration, response) %>%
    mutate(nucleotide = "TNP-ATP", experimenter = "Sam")

pcf_data <-
    read_csv("data/pcf_data.csv") %>%
    filter(dye < 480 | is.na(dye)) %>%
    select(unique_experiment_id, method, measure, construct, concentration, response) %>%
    mutate(nucleotide = "TNP-ATP", experimenter = "Sam")

concresp <-
    bind_rows(currents, fluorescence, pcf_data) %>%
    mutate(log_concentration = log10(concentration) %>%
    filter(log_concentration > -Inf) %>%
    group_by(unique_experiment_id, experimenter, method, measure, construct, nucleotide) %>%
    mutate(n = n()) %>%
    filter(n > 3)

hill_equation <- function(log_concentration, ec50, hill, floor){
    floor + ((1 - floor) / (1 + 10^((ec50 - log_concentration) * -hill)))
}

population_hill_fits_free_floor <-
    concresp %>%
    group_by(method, measure, construct, nucleotide) %>%
    nest() %>%
    mutate(fit = map(data, ~ nls_multstart(
            response ~ hill_equation(log_concentration, ec50, hill, floor),
            data = data.frame(.),
            iter = 500,
            lower = c(-Inf, 0, 0),
            start_lower = c(ec50 = -8, hill = 0, floor= 0),
            start_upper = c(ec50 = 0, hill = 2, floor = 1),
            na.action = na.omit
              )
    )
    )

population_hill_fits_fixed_floor <-
    concresp %>%
    group_by(method, measure, construct, nucleotide) %>%
    nest() %>%
    mutate(fit = map(data, ~ nls_multstart(
            response ~ hill_equation(log_concentration, ec50, hill, 0),
            data = data.frame(.),
            iter = 500,
            lower = c(-Inf, 0),
            start_lower = c(ec50 = -8, hill = 0),
            start_upper = c(ec50 = 0, hill = 2),
            na.action = na.omit
              )
    )
    )

individual_hill_fits_free_floor <-
    concresp %>%
    group_by(unique_experiment_id, method, measure, construct, nucleotide) %>%
    nest() %>%
    mutate(fit = map(data, ~ nls_multstart(
            response ~ hill_equation(log_concentration, ec50, hill, floor),
            data = data.frame(.),
            iter = 500,
            lower = c(-Inf, 0, 0),
            start_lower = c(ec50 = -8, hill = 0, floor= 0),
            start_upper = c(ec50 = 0, hill = 2, floor = 1),
            na.action = na.omit
              )
    )
    )

individual_hill_fits_fixed_floor <-
    concresp %>%
    group_by(unique_experiment_id, method, measure, construct, nucleotide) %>%
    nest() %>%
    mutate(fit = map(data, ~ nls_multstart(
            response ~ hill_equation(log_concentration, ec50, hill, 0),
            data = data.frame(.),
            iter = 500,
            lower = c(-Inf, 0),
            start_lower = c(ec50 = -8, hill = 0),
            start_upper = c(ec50 = 0, hill = 2),
            na.action = na.omit
              )
    )
    )

saveRDS(population_hill_fits_free_floor, "data/population_hill_fits_free_floor.rds")
saveRDS(population_hill_fits_fixed_floor, "data/population_hill_fits_fixed_floor.rds")
saveRDS(individual_hill_fits_free_floor, "data/individual_hill_fits_free_floor.rds")
saveRDS(individual_hill_fits_fixed_floor, "data/individual_hill_fits_fixed_floor.rds")
