library(tidyverse)
library(broom)
library(nls.multstart)
library(brms)
library(tidybayes)
library(ggbeeswarm)
library(ggdist)
library(distributional)

concresp_1 <-
    read_csv("data/combined_drc_data.csv") %>%
    mutate(log_concentration = log10(concentration)) %>%
    filter(log_concentration > -Inf) %>%
    group_by(unique_experiment_id) %>%
    mutate(observations = n()) %>%
    filter(observations >= 3) %>%
    ungroup() %>%
    mutate(
        response = case_when(
            response < -0.1 ~ -0.1,
            response > 1.1 ~ 1.1,
            TRUE ~ response
            )
        ) %>%
    drop_na() %>%
    filter(construct == "W311*-GFP+SUR", nucleotide == "TNP-ATP")

concresp_2 <-
    read_csv("data/activation_data.csv") %>%
    mutate(log_concentration = log10(concentration)) %>%
    filter(log_concentration > -Inf)

concresp <- bind_rows(concresp_1, concresp_2)

currents <-
    concresp %>%
    filter(measure == "current")

fluorescence <-
    concresp %>%
    filter(measure == "fluorescence")

ggplot(concresp) +
geom_quasirandom(aes(x = concentration, y = response, colour = measure, shape = nucleotide)) +
scale_x_log10() +
facet_grid(rows = vars(method))