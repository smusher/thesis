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
brms_iter <- 5000
brms_warmup <- 2000
brms_chains <- 4
brms_thin <- 1
brms_seed <- 2020

concresp <-
    read_csv("data/combined_drc_data.csv") %>%
    mutate(log_concentration = log10(concentration)) %>%
    filter(log_concentration > -Inf) %>%
    group_by(unique_experiment_id) %>%
    mutate(observations = n()) %>%
    filter(observations >= 3)

currents <-
    concresp %>%
    filter(measure == "current") %>%
    ungroup() %>%
    mutate(
        construct = fct_relevel(factor(construct), "WT-GFP+SUR"),
        nucleotide = fct_relevel(factor(nucleotide), "ATP")
        )

fluorescence <-
    concresp %>%
    filter(measure == "fluorescence") %>%
    ungroup() %>%
    mutate(
        construct = fct_relevel(factor(construct), "W311*-GFP+SUR"),
        method = fct_relevel(factor(method), "unroofed")
        )

hill_priors <-
    c(
        prior(uniform(-8, 0), nlpar = "ec50", class = b, lb = -8, ub = 0),
        prior(uniform(0, 2), nlpar = "hill", class = b, lb = 0, ub = 2),
        prior(uniform(0, 1), nlpar = "floor", class = b, lb = 0, ub = 1),
        prior(cauchy(0, 1), class = sd, nlpar = "ec50"),
        prior(cauchy(0, 1), class = sd, nlpar = "hill"),
        prior(cauchy(0, 1), class = sd, nlpar = "floor"),
        prior(cauchy(0, 3), class = sigma)
    )

hill_priors_2 <-
    c(
        prior(uniform(-8, 0), nlpar = "ec50", class = b, lb = -8, ub = 0),
        prior(uniform(0, 2), nlpar = "hill", lb = 0, ub = 2),
        prior(uniform(0, 1), nlpar = "floor", lb = 0, ub = 1),
        prior(cauchy(0, 1), class = sd, nlpar = "ec50"),
        prior(cauchy(0, 3), class = sigma)
    )

hill_priors_prep <-
    c(
        prior(uniform(-8, 0), nlpar = "ec50", class = b, lb = -8, ub = 0),
        prior(uniform(0, 2), nlpar = "hill", class = b, lb = 0, ub = 2),
        prior(uniform(0, 1), nlpar = "floor", class = b, lb = 0, ub = 1),
        prior(cauchy(0, 3), class = sigma)
    )

hill_priors_3 <-
    c(
        prior(uniform(-8, 0), nlpar = "ec50", class = b, lb = -8, ub = 0),
        prior(uniform(0, 2), nlpar = "hill", class = b, lb = 0, ub = 2),
        prior(cauchy(0, 1), class = sd, nlpar = "ec50"),
        prior(cauchy(0, 1), class = sd, nlpar = "hill"),
        prior(cauchy(0, 3), class = sigma)
    )

hill_priors_4 <-
    c(
        prior(uniform(-8, 0), nlpar = "ec50", class = b, lb = -8, ub = 0),
        prior(uniform(0, 2), nlpar = "hill", class = b, lb = 0, ub = 2),
        prior(cauchy(0, 1), class = sd, nlpar = "ec50"),
        prior(cauchy(0, 3), class = sigma)
    )

brm_hill_eq_1a <-
    bf(#without correlations
        response ~ floor + ((1 - floor) / (1 + 10^((ec50 - log_concentration) * -hill))),
        ec50 + hill + floor ~ 0 + Intercept + (1||construct:nucleotide),
        nl = TRUE
        )

brm_hill_eq_1aa <-
    bf(#without correlations
        response ~ floor + ((1 - floor) / (1 + 10^((ec50 - log_concentration) * -hill))),
        ec50 + hill + floor ~ 0 + construct:nucleotide + (construct:nucleotide|n),
        nl = TRUE
        )

brm_hill_eq_1b <-
    bf(#group effects only for ec50
        response ~ floor + ((1 - floor) / (1 + 10^((ec50 - log_concentration) * -hill))),
        ec50 ~ 0 + Intercept + (1||construct:nucleotide),
        floor + hill ~ 0 + construct:nucleotide,
        nl = TRUE
        )

brm_hill_eq_1c <-
    bf(#hill fixed between constructs
        response ~ floor + ((1 - floor) / (1 + 10^((ec50 - log_concentration) * -hill))),
        ec50 ~ 0 + Intercept + (1||construct:nucleotide),
        floor ~ 0 + construct:nucleotide,
        hill ~ 1,
        nl = TRUE
        )

brm_hill_eq_2a <-
    bf(#without correlations
        response ~ floor + ((1 - floor) / (1 + 10^((ec50 - log_concentration) * -hill))),
        ec50 + hill + floor ~ 0 + Intercept + (1||construct:method),
        nl = TRUE
        )

brm_hill_eq_2b <-
    bf(#group effects only for ec50
        response ~ floor + ((1 - floor) / (1 + 10^((ec50 - log_concentration) * -hill))),
        ec50 ~ 0 + Intercept + (1||construct:method),
        floor + hill ~ 0 + construct:method,
        nl = TRUE
        )

brm_hill_eq_2c <-
    bf(#hill fixed between constructs
        response ~ floor + ((1 - floor) / (1 + 10^((ec50 - log_concentration) * -hill))),
        ec50 ~ 0 + Intercept + (1||construct:method),
        floor ~ 0 + construct:method,
        hill ~ 1,
        nl = TRUE
        )

brm_hill_eq_3_prep <-
    bf(#calculate maximum FRET for control construct in unroofed
        log2(response + 1) ~ floor + ((1 - floor) / (1 + 10^((ec50 - log_concentration) * -hill))),
        ec50 ~ 1,
        floor ~ 1,
        hill ~ 1,
        nl = TRUE
        )

brm_hill_eq_3a <-
    bf(#without correlations
        log2(response + 1) ~ 0.07 + ((1 - 0.07) / (1 + 10^((ec50 - log_concentration) * -hill))),
        ec50 + hill ~ 0 + Intercept + (1||construct:method),
        nl = TRUE
        )

brm_hill_eq_3b <-
    bf(#group effects only for ec50
        log2(response + 1) ~ 0.07 + ((1 - 0.07) / (1 + 10^((ec50 - log_concentration) * -hill))),
        ec50 ~ 0 + Intercept + (1||construct:method),
        hill ~ 0 + construct:method,
        nl = TRUE
        )

brm_hill_eq_3c <-
    bf(#hill fixed between constructs
        log2(response + 1) ~ 0.07 + ((1 - 0.07) / (1 + 10^((ec50 - log_concentration) * -hill))),
        ec50 ~ 0 + Intercept + (1||construct:method),
        hill ~ 1,
        nl = TRUE
        )

hill_brm_fit_1a <-
    brm(
        data = currents,
        family = gaussian(),
        prior = hill_priors,
        formula = brm_hill_eq_1a,
        control = list(adapt_delta = 0.99, max_treedepth=15),
        iter = brms_iter,
        warmup = brms_warmup,
        chains = brms_chains,
        thin = brms_thin,
        seed = brms_seed,
        file = "data/hill_fits/1a_currents.rds",
        sample_prior = "yes",
        save_pars = save_pars(all=TRUE),
        cores = 4
        ) %>%
    add_criterion(c("loo", "waic"))

hill_brm_fit_1b <-
    brm(
        data = currents,
        family = gaussian(),
        prior = hill_priors_2,
        formula = brm_hill_eq_1b,
        control = list(adapt_delta = 0.99, max_treedepth=15),
        iter = brms_iter,
        warmup = brms_warmup,
        chains = brms_chains,
        thin = brms_thin,
        seed = brms_seed,
        file = "data/hill_fits/1b_currents.rds",
        sample_prior = "yes",
        save_pars = save_pars(all=TRUE),
        cores = 4
        ) %>%
    add_criterion(c("loo", "waic"))

hill_brm_fit_1c <-
    brm(
        data = currents,
        family = gaussian(),
        prior = hill_priors_2,
        formula = brm_hill_eq_1c,
        control = list(adapt_delta = 0.99, max_treedepth=15),
        iter = brms_iter,
        warmup = brms_warmup,
        chains = brms_chains,
        thin = brms_thin,
        seed = brms_seed,
        file = "data/hill_fits/1c_currents.rds",
        sample_prior = "yes",
        save_pars = save_pars(all=TRUE),
        cores = 4
        ) %>%
    add_criterion(c("loo", "waic"))

hill_brm_fit_2a <-
    brm(
        data = fluorescence,
        family = gaussian(),
        prior = hill_priors,
        formula = brm_hill_eq_2a,
        control = list(adapt_delta = 0.99, max_treedepth=15),
        iter = brms_iter,
        warmup = brms_warmup,
        chains = brms_chains,
        thin = brms_thin,
        seed = brms_seed,
        file = "data/hill_fits/2a_fluorescence.rds",
        sample_prior = "yes",
        save_pars = save_pars(all=TRUE),
        cores = 4
        ) %>%
    add_criterion(c("loo", "waic"))

hill_brm_fit_2b <-
    brm(
        data = fluorescence,
        family = gaussian(),
        prior = hill_priors_2,
        formula = brm_hill_eq_2b,
        control = list(adapt_delta = 0.99, max_treedepth=15),
        iter = brms_iter,
        warmup = brms_warmup,
        chains = brms_chains,
        thin = brms_thin,
        seed = brms_seed,
        file = "data/hill_fits/2b_fluorescence.rds",
        sample_prior = "yes",
        save_pars = save_pars(all=TRUE),
        cores = 4
        ) %>%
    add_criterion(c("loo", "waic"))

hill_brm_fit_2c <-
    brm(
        data = fluorescence,
        family = gaussian(),
        prior = hill_priors_2,
        formula = brm_hill_eq_2c,
        control = list(adapt_delta = 0.99, max_treedepth=15),
        iter = brms_iter,
        warmup = brms_warmup,
        chains = brms_chains,
        thin = brms_thin,
        seed = brms_seed,
        file = "data/hill_fits/2c_fluorescence.rds",
        sample_prior = "yes",
        save_pars = save_pars(all=TRUE),
        cores = 4
        ) %>%
    add_criterion(c("loo", "waic"))

hill_brm_fit_3_prep <-
    brm(
        data = fluorescence %>% filter(construct == "W311*-GFP+SUR", method == "unroofed"),
        family = gaussian(),
        prior = hill_priors_prep,
        formula = brm_hill_eq_3_prep,
        control = list(adapt_delta = 0.99, max_treedepth=15),
        iter = brms_iter,
        warmup = brms_warmup,
        chains = brms_chains,
        thin = brms_thin,
        seed = brms_seed,
        file = "data/hill_fits/fluorescence_prep.rds",
        sample_prior = "yes",
        save_pars = save_pars(all=TRUE),
        cores = 4
        )

hill_brm_fit_3a <-
    brm(
        data = fluorescence,
        family = gaussian(),
        prior = hill_priors_3,
        formula = brm_hill_eq_3a,
        control = list(adapt_delta = 0.99, max_treedepth=15),
        iter = brms_iter,
        warmup = brms_warmup,
        chains = brms_chains,
        thin = brms_thin,
        seed = brms_seed,
        file = "data/hill_fits/3a_fluorescence.rds",
        sample_prior = "yes",
        save_pars = save_pars(all=TRUE),
        cores = 4
        ) %>%
    add_criterion(c("loo", "waic"))

hill_brm_fit_3b <-
    brm(
        data = fluorescence,
        family = gaussian(),
        prior = hill_priors_4,
        formula = brm_hill_eq_3b,
        control = list(adapt_delta = 0.99, max_treedepth=15),
        iter = brms_iter,
        warmup = brms_warmup,
        chains = brms_chains,
        thin = brms_thin,
        seed = brms_seed,
        file = "data/hill_fits/3b_fluorescence.rds",
        sample_prior = "yes",
        save_pars = save_pars(all=TRUE),
        cores = 4
        ) %>%
    add_criterion(c("loo", "waic"))

hill_brm_fit_3c <-
    brm(
        data = fluorescence,
        family = gaussian(),
        prior = hill_priors_4,
        formula = brm_hill_eq_3c,
        control = list(adapt_delta = 0.99, max_treedepth=15),
        iter = brms_iter,
        warmup = brms_warmup,
        chains = brms_chains,
        thin = brms_thin,
        seed = brms_seed,
        file = "data/hill_fits/3c_fluorescence.rds",
        sample_prior = "yes",
        save_pars = save_pars(all=TRUE),
        cores = 4
        ) %>%
    add_criterion(c("loo", "waic"))

currents_2 <-
    currents %>%
    mutate(anap_present = case_when(stringr::str_detect(construct, "W311*") == TRUE ~ TRUE, TRUE ~ FALSE))

currents_2 %>%
data_grid(
    log_concentration = seq_range(-7:-2, n = 51),
    nesting(construct, nucleotide)
    ) %>%
add_predicted_draws(hill_brm_fit_1b, re_formula = NULL) %>%
group_by(construct, nucleotide, log_concentration) %>%
point_interval(.prediction, .point = median, .interval = qi, .width = .95) %>%
mutate(anap_present = case_when(stringr::str_detect(construct, "W311*") == TRUE ~ TRUE, TRUE ~ FALSE)) -> fitted_draws

w311 <-
    fitted_draws %>%
    filter(construct == "W311*-GFP+SUR") %>%
    ungroup() %>%
    select(-construct)

wt <-
    fitted_draws %>%
    filter(construct == "WT-GFP+SUR") %>%
    ungroup() %>%
    select(-construct)

ggplot() +
geom_ribbon(data = wt, aes(x = log_concentration, ymin = .lower, ymax = .upper), fill = "#636363", alpha = 0.5) +
geom_ribbon(data = fitted_draws %>% filter(anap_present == FALSE), aes(x = log_concentration, ymin = .lower, ymax = .upper), fill = "#3182bd", alpha = 0.5) +
geom_line(data = fitted_draws %>% filter(anap_present == FALSE), aes(x = log_concentration, y = .prediction)) +
geom_quasirandom(data = currents_2 %>% filter(anap_present == FALSE), aes(x = log_concentration, y = response), shape = 21, fill = "white", size = 2.5, width = 0.2) +
facet_grid(rows = vars(nucleotide), cols = vars(construct)) +
cowplot::theme_cowplot(font_family = "IBM Plex Sans Condensed", font_size = 9) +
theme(strip.background  = element_rect(fill = "white"), legend.position = "none") +
scale_y_continuous(breaks = c(0.0, 0.5, 1.0)) +
labs(
    x = "[nucleotide] (M)",
    y = expression(I/I[max])
    )

ggplot() +
geom_ribbon(data = w311, aes(x = log_concentration, ymin = .lower, ymax = .upper), fill = "#636363", alpha = 0.5) +
geom_ribbon(data = fitted_draws %>% filter(anap_present == TRUE), aes(x = log_concentration, ymin = .lower, ymax = .upper), fill = "#3182bd", alpha = 0.5) +
geom_line(data = fitted_draws %>% filter(anap_present == TRUE), aes(x = log_concentration, y = .prediction)) +
geom_quasirandom(data = currents_2 %>% filter(anap_present == TRUE), aes(x = log_concentration, y = response), shape = 21, fill = "white", size = 2.5, width = 0.2) +
facet_grid(rows = vars(nucleotide), cols = vars(construct)) +
cowplot::theme_cowplot(font_family = "IBM Plex Sans Condensed", font_size = 9) +
theme(strip.background  = element_rect(fill = "white"), legend.position = "none") +
scale_y_continuous(breaks = c(0.0, 0.5, 1.0)) +
labs(
    x = "[nucleotide] (M)",
    y = expression(I/I[max])
    )

fluorescence %>%
data_grid(
    log_concentration = seq_range(-7:-2, n = 51),
    nesting(construct, method)
    ) %>%
add_fitted_draws(hill_brm_fit_3c, re_formula = NULL) %>%
group_by(construct, method, log_concentration) %>%
point_interval(.value, .point = median, .interval = qi, .width = .95) -> fitted_draws

fluorescence %>%
data_grid(
    log_concentration = seq_range(-7:-2, n = 51),
    nesting(construct, method)
    ) %>%
add_predicted_draws(hill_brm_fit_3c, re_formula = NULL) %>%
group_by(construct, method, log_concentration) %>%
point_interval(.prediction, .point = median, .interval = qi, .width = .95) -> predicted_draws

w311 <-
    fitted_draws %>%
    filter(construct == "W311*-GFP+SUR") %>%
    ungroup() %>%
    select(-construct)

ggplot() +
geom_ribbon(data = w311, aes(x = log_concentration, ymin = .lower, ymax = .upper), fill = "#636363", alpha = 0.5) +
geom_ribbon(data = fitted_draws, aes(x = log_concentration, ymin = .lower, ymax = .upper), fill = "#3182bd", alpha = 0.5) +
geom_ribbon(data = predicted_draws, aes(x = log_concentration, ymin = .lower, ymax = .upper, colour = "black"), fill = NA, linetype = 2) +
geom_line(data = fitted_draws, aes(x = log_concentration, y = .value)) +
geom_quasirandom(data = fluorescence, aes(x = log_concentration, y = response), shape = 21, fill = "white", size = 2.5, width = 0.2) +
facet_grid(rows = vars(method), cols = vars(construct)) +
cowplot::theme_cowplot(font_family = "IBM Plex Sans Condensed", font_size = 9) +
theme(strip.background  = element_rect(fill = "white"), legend.position = "none") +
scale_y_continuous(breaks = c(0.0, 0.5, 1.0)) +
labs(
    x = "[nucleotide] (M)",
    y = expression(I/I[max])
    )

get_variables(hill_brm_fit_1a)

hill_brm_fit_1a %>%
recover_types(currents) %>%
spread_draws(`b_.*`, `r_construct:nucleotide__.*`[i,j], regex = TRUE) %>%
ungroup() %>%
mutate(i = case_when(i == "W311*"~str_c(i, j, sep = ","), TRUE~i)) %>%
separate(i, into = c("construct", "nucleotide"), sep = "_(?=A)|_(?=M)|_(?=T)") %>%
group_by(construct, nucleotide, .draw) %>%
summarise(
    ec50 = b_ec50_Intercept + `r_construct:nucleotide__ec50`,
    hill = b_hill_Intercept + `r_construct:nucleotide__hill`,
    floor = b_floor_Intercept + `r_construct:nucleotide__floor`
    ) %>%
pivot_longer(
    cols = ec50:floor,
    names_to = "parameter",
    values_to = "value"
    ) -> draws

control_draws <-
    draws %>%
    ungroup() %>%
    filter(
        construct %in% c("WT-GFP+SUR", "W311*-GFP+SUR"),
        nucleotide %in% c("ATP", "TNP-ATP")
        ) %>%
    pivot_wider(
        names_from = c("construct", "nucleotide"),
        values_from = "value"
        )

contrasted_draws <-
    left_join(draws, control_draws) %>%
    mutate(
        construct = factor(construct),
        nucleotide = factor(nucleotide),
        parameter = factor(parameter)
        ) %>%
    group_by(construct, nucleotide, .draw, parameter) %>%
    summarise(across(contains("GFP"), ~value - .x)) %>%
    pivot_longer(
        cols = !c(construct, nucleotide, .draw, parameter),
        names_to = "contrast",
        values_to = "value"
        ) %>%
    mutate(
        contrast = factor(contrast)
        )

contrasted_draws_summary <-
    contrasted_draws %>%
    group_by(construct, nucleotide, parameter, contrast) %>%
    summarise(
        sigma = sd(value),
        mu = mean(value)
        )

ggplot(
    contrasted_draws_summary %>% filter(parameter == "ec50", nucleotide != "MG-ATP"),
    aes(
        y = construct,
        dist = dist_normal(mu, sigma),
        fill = nucleotide,
        fill_ramp =  stat(cut_cdf_qi(cdf, .width = c(.5, .8, .95), labels = scales::percent_format()))
        )
    ) +
stat_dist_slab() +
facet_grid(cols = vars(contrast), scales = "free") +
scale_colour_brewer(palette = "PuOr") +
scale_fill_ramp_discrete(range = c(1, 0.2), na.translate = FALSE) +
labs(fill_ramp = "Interval")

ggplot(
    contrasted_draws_summary %>% filter(parameter == "ec50", nucleotide != "MG-ATP", construct %in% c("WT-GFP+SUR", "W311*-GFP+SUR")),
    aes(
        y = construct,
        dist = dist_normal(mu, sigma),
        fill = nucleotide,
        fill_ramp =  stat(cut_cdf_qi(cdf, .width = c(.5, .8, .95), labels = scales::percent_format()))
        )
    ) +
stat_dist_slab() +
facet_grid(rows = vars(contrast), scales = "free") +
scale_colour_brewer(palette = "PuOr") +
scale_fill_ramp_discrete(range = c(1, 0.2), na.translate = FALSE) +
labs(fill_ramp = "Interval")

get_variables(hill_brm_fit_3c)

hill_brm_fit_3c %>%
recover_types(fluorescence) %>%
spread_draws(`b_.*`, `r_construct:method__ec50`[i,j], regex = TRUE) %>%
ungroup() %>%
mutate(i = case_when(i == "W311*"~str_c(i, j, sep = ","), TRUE~i)) %>%
separate(i, into = c("construct", "method"), sep = "_(?=u)|_(?=p)") %>%
group_by(construct, method, .draw) %>%
summarise(
    ec50 = b_ec50_Intercept + `r_construct:method__ec50`
    ) %>%
pivot_longer(
    cols = ec50,
    names_to = "parameter",
    values_to = "value"
    ) -> draws

control_draws <-
    draws %>%
    ungroup() %>%
    filter(
        construct == "W311*-GFP+SUR"
        ) %>%
    pivot_wider(
        names_from = c("construct", "method"),
        values_from = "value"
        )

contrasted_draws <-
    left_join(draws, control_draws) %>%
    mutate(
        construct = factor(construct),
        method = factor(method),
        parameter = factor(parameter)
        ) %>%
    group_by(construct, method, .draw, parameter) %>%
    summarise(across(contains("GFP"), ~value - .x)) %>%
    pivot_longer(
        cols = !c(construct, method, .draw, parameter),
        names_to = "contrast",
        values_to = "value"
        ) %>%
    mutate(
        contrast = factor(contrast)
        )

contrasted_draws_summary <-
    contrasted_draws %>%
    group_by(construct, method, parameter, contrast) %>%
    summarise(
        sigma = sd(value),
        mu = mean(value)
        )

ggplot(
    contrasted_draws_summary,
    aes(
        y = construct,
        dist = dist_normal(mu, sigma),
        fill = method,
        fill_ramp =  stat(cut_cdf_qi(cdf, .width = c(.5, .8, .95), labels = scales::percent_format()))
        )
    ) +
stat_dist_slab() +
facet_grid(cols = vars(contrast), scales = "free") +
scale_colour_brewer(palette = "PuOr") +
scale_fill_ramp_discrete(range = c(1, 0.2), na.translate = FALSE) +
labs(fill_ramp = "Interval")
