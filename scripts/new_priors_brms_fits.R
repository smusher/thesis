library(tidyverse)
library(broom)
library(brms)
library(rstan)
library(tidybayes)
library(future)
source("scripts/specify_model_formulae.R")
brms_iter <- 30000
brms_warmup <- 20000
brms_chains <- 4
brms_thin <- 10
brms_seed <- 2019

pcf_data <-
    read_csv("data/pcf_data.csv") %>%
    filter(dye < 480 | is.na(dye), concentration > 0) %>%
    select(unique_experiment_id, construct, measure, method, concentration, response) %>%
    mutate(
        response = case_when(measure == "fluorescence" ~ 1 - response, TRUE ~ response),
        binding_mask = case_when(measure == "fluorescence" ~ 1, TRUE ~ 0)
    ) %>%
    group_by(unique_experiment_id, construct, measure) %>%
    mutate(response_normalised = case_when(
        measure == "fluorescence" ~ (response-min(response))/(max(response)-min(response)),
        measure == "current" ~ response
    )) %>%
    ungroup()

mwc_formula_transformed <-
    mwc_formula %>%
    str_replace_all("Ka", "(10^logKa)") %>%
    str_replace_all("L", "(10^logL)")

mwc_population <-
    bf(
        as.formula(mwc_formula_transformed),
        logL ~ 1,
        D ~ 1,
        logKa ~ 1,
        nl = TRUE
    )

coop_formula_transformed <-
    mwc_formula_with_cooperativity %>%
    str_replace_all("Ka", "(10^logKa)") %>%
    str_replace_all("L", "(10^logL)")

coop_population <-
    bf(
        as.formula(coop_formula_transformed),
        logL ~ 1,
        D ~ 1,
        logKa ~ 1,
        c ~ 1,
        nl = TRUE
    )

mwc_prior_prep <- 
    readRDS("data/mwc_fits/W311*-GFP+SUR_mwc_population.rds") %>%
    spread_draws(`b_.*`, regex = T) %>%
    select(starts_with("b_"))

logL_dist <-
    MASS::fitdistr(mwc_prior_prep$b_logL_Intercept, "normal")

logKa_dist <-
    MASS::fitdistr(mwc_prior_prep$b_logKa_Intercept, "normal")

D_dist <-
    MASS::fitdistr(mwc_prior_prep$b_D_Intercept, "normal")

mwc_priors <-
    c(
        prior(normal(0.056, 0.051), nlpar = "D", lb = 0, ub = 1),
        prior(normal(4.321, 0.056), nlpar = "logKa"),
        prior(normal(-1.076, 0.349), nlpar = "logL")
    )

mwc_broader_priors <-
    c(
        prior(normal(0.05, 0.5), nlpar = "D", lb = 0, ub = 1),
        prior(normal(4.32, 0.5), nlpar = "logKa"),
        prior(normal(-1.08, 3.5), nlpar = "logL")
    )

coop_prior_prep <- 
    readRDS("data/mwc_fits/W311*-GFP+SUR_coop_population.rds") %>%
    spread_draws(`b_.*`, regex = T) %>%
    select(starts_with("b_"))

logL_dist <-
    MASS::fitdistr(coop_prior_prep$b_logL_Intercept, "normal")

logKa_dist <-
    MASS::fitdistr(coop_prior_prep$b_logKa_Intercept, "normal")

D_dist <-
    MASS::fitdistr(coop_prior_prep$b_D_Intercept, "normal")

c_dist <-
    MASS::fitdistr(coop_prior_prep$b_c_Intercept, "normal")

coop_priors <-
    c(
        prior(normal(0.149, 0.073), nlpar = "D", lb = 0, ub = 1),
        prior(normal(0.176, 0.076), nlpar = "c", lb = 0, ub = 1),
        prior(normal(4.839, 0.187), nlpar = "logKa"),
        prior(normal(-0.431, 0.475), nlpar = "logL")
    )

coop_broader_priors <-
    c(
        prior(normal(0.15, 0.7), nlpar = "D", lb = 0, ub = 1),
        prior(normal(0.18, 0.8), nlpar = "c", lb = 0, ub = 1),
        prior(normal(4.84, 1.9), nlpar = "logKa"),
        prior(normal(-0.43, 4.8), nlpar = "logL")
    )

do_fit <- function(construct_name, formula, priors){
    model_name <- deparse(substitute(formula))
    save_name <- paste(construct_name, model_name, sep = "_")
    print(save_name)
    data <- pcf_data %>% filter(construct == construct_name)
    fit <-
        brm(
            formula = formula,
            data = data,
            prior = priors,
            family = gaussian(),
            iter = brms_iter,
            warmup = brms_warmup,
            chains = brms_chains,
            thin = brms_thin,
            seed = brms_seed,
            control = list(adapt_delta = 0.99, max_treedepth = 15),
            file = paste("data/mwc_fits_new_priors/", save_name, sep = ""),
            sample_prior = "yes",
            future = TRUE
        )
    return(fit)
}

mwc_cs_pop_prior <- do_fit("W311*,C166S-GFP+SUR", mwc_population, mwc_priors)
mwc_ka_pop_prior <- do_fit("W311*-GFP+SUR-K205A", mwc_population, mwc_priors)
mwc_ke_pop_prior <- do_fit("W311*-GFP+SUR-K205E", mwc_population, mwc_priors)
coop_cs_pop_prior <- do_fit("W311*,C166S-GFP+SUR", coop_population, coop_priors)
coop_ka_pop_prior <- do_fit("W311*-GFP+SUR-K205A", coop_population, coop_priors)
coop_ke_pop_prior <- do_fit("W311*-GFP+SUR-K205E", coop_population, coop_priors)

do_fit <- function(construct_name, formula, priors){
    model_name <- deparse(substitute(formula))
    save_name <- paste(construct_name, model_name, "broader", sep = "_")
    print(save_name)
    data <- pcf_data %>% filter(construct == construct_name)
    fit <-
        brm(
            formula = formula,
            data = data,
            prior = priors,
            family = gaussian(),
            iter = brms_iter,
            warmup = brms_warmup,
            chains = brms_chains,
            thin = brms_thin,
            seed = brms_seed,
            control = list(adapt_delta = 0.99, max_treedepth = 15),
            file = paste("data/mwc_fits_new_priors/", save_name, sep = ""),
            sample_prior = "yes",
            future = TRUE
        )
    return(fit)
}

mwc_cs_pop_broader_prior <- do_fit("W311*,C166S-GFP+SUR", mwc_population, mwc_broader_priors)
mwc_ka_pop_broader_prior <- do_fit("W311*-GFP+SUR-K205A", mwc_population, mwc_broader_priors)
mwc_ke_pop_broader_prior <- do_fit("W311*-GFP+SUR-K205E", mwc_population, mwc_broader_priors)
coop_cs_pop_broader_prior <- do_fit("W311*,C166S-GFP+SUR", coop_population, coop_broader_priors)
coop_ka_pop_broader_prior <- do_fit("W311*-GFP+SUR-K205A", coop_population, coop_broader_priors)
coop_ke_pop_broader_prior <- do_fit("W311*-GFP+SUR-K205E", coop_population, coop_broader_priors)
