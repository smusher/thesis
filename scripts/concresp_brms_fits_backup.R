library(tidyverse)
library(broom)
library(brms)
library(tidybayes)
library(bayesplot)
library(modelr)
library(patchwork)
library(ggbeeswarm)
library(ggdist)
library(distributional)
brms_iter <- 10000
brms_warmup <- 2000
brms_chains <- 4
brms_cores <- 4
brms_seed <- 2021

concresp <-
    read_csv("data/combined_drc_data.csv") %>%
    mutate(log_concentration = log10(concentration)) %>%
    filter(log_concentration > -Inf) %>%
    group_by(unique_experiment_id) %>%
    mutate(observations = n()) %>%
    filter(observations >= 3)

currents <-
    concresp %>%
    filter(measure == "current", construct %in% c("WT-GFP+SUR", "W311*-GFP+SUR"), nucleotide %in% c("ATP", "TNP-ATP")) %>%
    ungroup() %>%
    mutate(
        construct = fct_relevel(factor(construct), "WT-GFP+SUR"),
        nucleotide = fct_relevel(factor(nucleotide), "ATP")
        )

fluorescence <-
    concresp %>%
    filter(measure == "fluorescence", construct == "W311*-GFP+SUR") %>%
    ungroup() %>%
    mutate(
        construct = fct_relevel(factor(construct), "W311*-GFP+SUR"),
        method = fct_relevel(factor(method), "unroofed")
        )

current_formula_1 <-
    bf(#estimate each experiment as a group level effect - no correlations
        response ~ floor + ((1 - floor) / (1 + 10^((ec50 - log_concentration) * -hill))),
        ec50 + hill + floor ~ 0 + construct:nucleotide + (1 + construct:nucleotide||unique_experiment_id),
        sigma ~ (1||unique_experiment_id),
        nl = TRUE
    )

current_formula_2 <-
    bf(#estimate each experiment as a group level effect only for ec50
        response ~ floor + ((1 - floor) / (1 + 10^((ec50 - log_concentration) * -hill))),
        ec50 ~ 0 + construct:nucleotide + (1 + construct:nucleotide||unique_experiment_id),
        hill + floor ~ 0 + construct:nucleotide,
        sigma ~ (1||unique_experiment_id),
        nl = TRUE
    )

current_formula_3 <-
    bf(#estimate each pool of experiments as a group level effect
        response ~ floor + ((1 - floor) / (1 + 10^((ec50 - log_concentration) * -hill))),
        ec50 + hill + floor ~ 0 + construct:nucleotide + (1||construct:nucleotide),
        sigma ~ (1||construct:nucleotide),
        nl = TRUE
    )

ec50_intercept_prior <- prior(uniform(-8, 0), nlpar = "ec50", class = b, lb = -8, ub = 0)
hill_intercept_prior <- prior(uniform(0, 2), nlpar = "hill", class = b, lb = 0, ub = 2)
floor_intercept_prior <- prior(uniform(0, 1), nlpar = "floor", class = b, lb = 0, ub = 1)
ec50_sd_prior <- prior(cauchy(0, 1), nlpar = "ec50", class = sd)
hill_sd_prior <- prior(cauchy(0, 1), nlpar = "hill", class = sd)
floor_sd_prior <- prior(cauchy(0, 1), nlpar = "floor", class = sd)
sigma_prior <- prior(cauchy(0, 3), dpar = "sigma", class = sd)

brm(
    data = currents,
    formula = current_formula_1,
    prior = c(
        ec50_intercept_prior,
        hill_intercept_prior,
        floor_intercept_prior,
        ec50_sd_prior,
        hill_sd_prior,
        floor_sd_prior,
        sigma_prior
        ),
    family = gaussian(),
    iter = brms_iter,
    warmup = brms_warmup,
    chains = brms_chains,
    cores = brms_cores,
    seed = brms_seed,
    save_pars = save_pars(all=TRUE),
    sample_prior = "yes",
    file = "data/hill_fits/current_fit_1.rds",
    control = list(adapt_delta = 0.99, max_treedepth=15)
    ) %>%
    add_criterion(c("waic", "loo")) -> test_1

brm(
    data = currents,
    formula = current_formula_2,
    prior = c(
        ec50_intercept_prior,
        hill_intercept_prior,
        floor_intercept_prior,
        ec50_sd_prior,
        sigma_prior
        ),
    family = gaussian(),
    iter = brms_iter,
    warmup = brms_warmup,
    chains = brms_chains,
    cores = brms_cores,
    seed = brms_seed,
    save_pars = save_pars(all=TRUE),
    sample_prior = "yes",
    file = "data/hill_fits/current_fit_2.rds",
    control = list(adapt_delta = 0.99, max_treedepth=15)
    ) %>%
    add_criterion(c("waic", "loo")) -> test_2

brm(
    data = currents,
    formula = current_formula_3,
    prior = c(
        ec50_intercept_prior,
        hill_intercept_prior,
        floor_intercept_prior,
        ec50_sd_prior,
        hill_sd_prior,
        floor_sd_prior,
        sigma_prior
        ),
    family = gaussian(),
    iter = brms_iter,
    warmup = brms_warmup,
    chains = brms_chains,
    cores = brms_cores,
    seed = brms_seed,
    save_pars = save_pars(all=TRUE),
    sample_prior = "yes",
    file = "data/hill_fits/current_fit_3.rds",
    control = list(adapt_delta = 0.99, max_treedepth=15)
    ) %>%
    add_criterion(c("waic", "loo")) -> test_3

brm(
    data = cfluorescence,
    formula = current_formula_1,
    prior = c(
        ec50_intercept_prior,
        hill_intercept_prior,
        floor_intercept_prior,
        ec50_sd_prior,
        hill_sd_prior,
        floor_sd_prior,
        sigma_prior
        ),
    family = gaussian(),
    iter = brms_iter,
    warmup = brms_warmup,
    chains = brms_chains,
    cores = brms_cores,
    seed = brms_seed,
    save_pars = save_pars(all=TRUE),
    sample_prior = "yes",
    file = "data/hill_fits/fluorescence_fit_1.rds",
    control = list(adapt_delta = 0.99, max_treedepth=15)
    ) %>%
    add_criterion(c("waic", "loo"), moment_match=TRUE) -> test_4

brm(
    data = cfluorescence,
    formula = current_formula_2,
    prior = c(
        ec50_intercept_prior,
        hill_intercept_prior,
        floor_intercept_prior,
        ec50_sd_prior,
        sigma_prior
        ),
    family = gaussian(),
    iter = brms_iter,
    warmup = brms_warmup,
    chains = brms_chains,
    cores = brms_cores,
    seed = brms_seed,
    save_pars = save_pars(all=TRUE),
    sample_prior = "yes",
    file = "data/hill_fits/fluorescence_fit_2.rds",
    control = list(adapt_delta = 0.99, max_treedepth=15)
    ) %>%
    add_criterion(c("waic", "loo"), moment_match=TRUE) -> test_5

brm(
    data = cfluorescence,
    formula = current_formula_3,
    prior = c(
        ec50_intercept_prior,
        hill_intercept_prior,
        floor_intercept_prior,
        ec50_sd_prior,
        hill_sd_prior,
        floor_sd_prior,
        sigma_prior
        ),
    family = gaussian(),
    iter = brms_iter,
    warmup = brms_warmup,
    chains = brms_chains,
    cores = brms_cores,
    seed = brms_seed,
    save_pars = save_pars(all=TRUE),
    sample_prior = "yes",
    file = "data/hill_fits/fluorescence_fit_3.rds",
    control = list(adapt_delta = 0.99, max_treedepth=15)
    ) %>%
    add_criterion(c("waic", "loo"), moment_match=TRUE) -> test_6
