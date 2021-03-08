library(tidyverse)
library(broom)
library(nls.multstart)
library(brms)

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

# hill_equation <- function(log_concentration, ec50, hill, floor){
#     floor + ((1 - floor) / (1 + 10^((ec50 - log_concentration) * -hill)))
# }

# #population fits
# unadjusted_fits_free_floor_pop <-
#     concresp %>%
#     group_by(construct, method, measure, nucleotide) %>%
#     nest() %>%
#     mutate(fit = map(data, ~ nls_multstart(
#             response ~ hill_equation(log_concentration, ec50, hill, floor),
#             data = data.frame(.),
#             iter = c(10,10,10),
#             lower = c(-Inf, 0, 0),
#             start_lower = c(ec50 = -7, hill = 0, floor= 0),
#             start_upper = c(ec50 = -1, hill = 2, floor = 1),
#             na.action = na.omit
#               )
#     )
#     )

# adjusted_fits_free_floor_pop <-
#     fluorescence %>%
#     group_by(construct, method, measure, nucleotide) %>%
#     nest() %>%
#     mutate(fit = map(data, ~ nls_multstart(
#             log2(response + 1) ~ hill_equation(log_concentration, ec50, hill, floor),
#             data = data.frame(.),
#             iter = c(10, 10, 10),
#             lower = c(-Inf, 0, 0),
#             start_lower = c(ec50 = -7, hill = 0, floor= 0),
#             start_upper = c(ec50 = -1, hill = 2, floor = 1),
#             na.action = na.omit
#               )
#     )
#     )

# adjusted_fits_fixed_floor_pop <-
#     fluorescence %>%
#     group_by(construct, method, measure, nucleotide) %>%
#     nest() %>%
#     mutate(fit = map(data, ~ nls_multstart(
#             log2(response + 1) ~ 0.1 + ((1 - 0.1) / (1 + 10^((ec50 - log_concentration) * -hill))),
#             data = data.frame(.),
#             iter = c(15, 15),
#             lower = c(-Inf, 0),
#             start_lower = c(ec50 = -7, hill = 0),
#             start_upper = c(ec50 = -2, hill = 2),
#             na.action = na.omit
#               )
#     )
#     )

# saveRDS(unadjusted_fits_free_floor_pop, "data/hill_fits/ua_free_fit_pop.rds")
# saveRDS(adjusted_fits_free_floor_pop, "data/hill_fits/a_free_fit_pop.rds")
# saveRDS(adjusted_fits_fixed_floor_pop, "data/hill_fits/a_fixed_fit_pop.rds")

# unadjusted_fits_free_floor_pop %>%
# mutate(tidied = map(fit, tidy)) %>%
# unnest(tidied) %>%
# select(construct, method, measure, nucleotide, term, estimate) %>%
# pivot_wider(names_from = term, values_from = estimate) %>%
# rename(pop_ec50 = ec50) %>%
# right_join(concresp) -> ua_free_tidy_pop

# adjusted_fits_free_floor_pop %>%
# mutate(tidied = map(fit, tidy)) %>%
# unnest(tidied) %>%
# select(construct, method, measure, nucleotide, term, estimate) %>%
# pivot_wider(names_from = term, values_from = estimate) %>%
# rename(pop_ec50 = ec50) %>%
# right_join(fluorescence) -> a_free_tidy_pop

# adjusted_fits_fixed_floor_pop %>%
# mutate(tidied = map(fit, tidy)) %>%
# unnest(tidied) %>%
# select(construct, method, measure, nucleotide, term, estimate) %>%
# pivot_wider(names_from = term, values_from = estimate) %>%
# rename(pop_ec50 = ec50) %>%
# mutate(floor = 0.1) %>%
# right_join(fluorescence) -> a_fixed_tidy_pop

# #individual fits
# unadjusted_fits_free_floor <-
#     ua_free_tidy_pop %>%
#     group_by(construct, method, measure, nucleotide, unique_experiment_id) %>%
#     nest() %>%
#     mutate(fit = map(data, ~ nls_multstart(
#             response ~ hill_equation(log_concentration, ec50, hill, floor),
#             data = data.frame(.),
#             iter = c(100),
#             lower = c(-Inf),
#             start_lower = c(ec50 = -7),
#             start_upper = c(ec50 = -1),
#             na.action = na.omit
#               )
#     )
#     )

# adjusted_fits_free_floor <-
#     a_free_tidy_pop %>%
#     group_by(construct, method, measure, nucleotide, unique_experiment_id) %>%
#     nest() %>%
#     mutate(fit = map(data, ~ nls_multstart(
#             log2(response + 1) ~ hill_equation(log_concentration, ec50, hill, floor),
#             data = data.frame(.),
#             iter = c(100),
#             lower = c(-Inf),
#             start_lower = c(ec50 = -7),
#             start_upper = c(ec50 = -1),
#             na.action = na.omit
#               )
#     )
#     )

# adjusted_fits_fixed_floor <-
#     a_fixed_tidy_pop %>%
#     group_by(construct, method, measure, nucleotide, unique_experiment_id) %>%
#     nest() %>%
#     mutate(fit = map(data, ~ nls_multstart(
#             log2(response + 1) ~ hill_equation(log_concentration, ec50, hill, floor),
#             data = data.frame(.),
#             iter = c(100),
#             lower = c(-Inf),
#             start_lower = c(ec50 = -7),
#             start_upper = c(ec50 = -1),
#             na.action = na.omit
#               )
#     )
#     )

# saveRDS(unadjusted_fits_free_floor, "data/hill_fits/ua_free_fit.rds")
# saveRDS(adjusted_fits_free_floor, "data/hill_fits/a_free_fit.rds")
# saveRDS(adjusted_fits_fixed_floor, "data/hill_fits/a_fixed_fit.rds")


free_brms_formula <-
    bf(
        response ~ floor + ((1 - floor) / (1 + 10^((ec50 - log_concentration) * -hill))),
        nl = TRUE,
        floor + ec50 + hill ~ 0 + construct:method:measure:nucleotide + (construct:method:measure:nucleotide||unique_experiment_id),
        sigma ~ (1||method:measure),
        family = gaussian()
        )

free_brms_priors <-
    c(
        #most likely near 0, tails off towards 1
        set_prior("uniform(0, 1)", nlpar = "floor", class = "b", lb = 0, ub = 1),
        #ec50 based on concentration range
        set_prior("normal(log(1e-4), 1)", nlpar = "ec50", class = "b"),
        #hill coefficient between 0 and 2
        set_prior("normal(1, 0.3)", nlpar = "hill", class = "b"),
        #standard cauchy prior for sigmas
        set_prior("cauchy(0, 1)", nlpar = "floor", class = "sd"),
        set_prior("cauchy(0, 1)", nlpar = "ec50", class = "sd"),
        set_prior("cauchy(0, 1)", nlpar = "hill", class = "sd"),
        set_prior("cauchy(0, 1)", dpar = "sigma", class = "Intercept"),
        set_prior("cauchy(0, 1)", dpar = "sigma", class = "sd")
        )

brms_iter <- 2000
brms_warmup <- 1000
brms_chains <- 4
brms_thin <- 1
brms_seed <- 2021

brm(
    formula = free_brms_formula,
    prior = free_brms_priors,
    data = concresp,
    chains = brms_chains,
    iter = brms_iter,
    warmup = brms_warmup,
    thin  = brms_thin,
    seed = brms_seed,
    control = list(adapt_delta = 0.95, max_treedepth = 10),
    cores = getOption("mc.cores", 4),
    file = "data/hill_fits/full_fit_070321_short",
    sample_prior = "yes",
    save_all_pars = TRUE
    ) %>%
    add_criterion("loo", moment_match=TRUE) -> test_run
