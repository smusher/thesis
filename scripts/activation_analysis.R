library(tidyverse)
library(broom)
library(nls.multstart)
library(brms)
library(tidybayes)
library(ggbeeswarm)
library(ggdist)
library(distributional)
library(bayesplot)

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
    filter(construct == "W311*-GFP+SUR", nucleotide == "TNP-ATP") %>%
    mutate(
        response = case_when(measure == "fluorescence" ~ (1 - (log2(response + 1))) / 0.9, TRUE ~ response),
        inhibition_mask = 1
        )

concresp_2 <-
    read_csv("data/activation_data.csv") %>%
    mutate(log_concentration = log10(concentration)) %>%
    filter(log_concentration > -Inf) %>%
    mutate(
        inhibition_mask = 0,
        unique_experiment_id = within_experiment_id
        )

concresp <-
    bind_rows(concresp_1, concresp_2) %>%
    mutate(
        binding_mask = case_when(measure == "fluorescence" ~ 1, measure == "current" ~ 0),
        method = case_when(method %in% c("pcf", "patch-clamp") ~ "patch", TRUE ~ method)
        )

ggplot(concresp) +
geom_quasirandom(aes(x = concentration, y = response, colour = measure, shape = method)) +
scale_x_log10() +
scale_shape_manual(values = c("unroofed" = 21, "patch" = 19)) +
facet_grid(cols = vars(construct))

inhibition_bound <-
"((Ka * Fa * ((1 + Ka * Fa) ^ 3)) + (L * D * Ka * Fa * ((1 + D * Ka * Fa) ^ 3))) / (((1 + Ka * Fa) ^ 4) + (L * ((1 + D * Ka * Fa) ^ 4)))"

inhibition_open <-
"((L * ((1 + D * Ka * Fa) ^ 4)) / (((1 + Ka * Fa) ^ 4) + (L * ((1 + D * Ka * Fa) ^ 4)))) / (L / (L + 1))"

activation_bound <-
"((Kb * Fb * ((1 + Kb * Fb) ^ 3)) + (L * E * Kb * Fb * ((1 + E * Kb * Fb) ^ 3))) / (((1 + Kb * Fb) ^ 4) + (L * ((1 + E * Kb * Fb) ^ 4)))"

activation_open <-
"(((L * ((1 + E * Kb * Fb) ^ 4)) / (((1 + Kb * Fb) ^ 4) + (L * ((1 + E * Kb * Fb) ^ 4)))) - (L / (L + 1))) / (((L * E^4) / (1 + L * E^4)) - (L / (L + 1)))"

inhibition_formula_str <-
paste(
    "(binding_mask * (",
    inhibition_bound,
    ")) + ((1 - binding_mask) * (",
    inhibition_open,
    "))",
    sep = ""
    ) %>%
    str_replace_all(c(
        "Fa" = "concentration",
        "Ka" = "exp(logKa)",
        "L" = "exp(logL)",
        "D" = "exp(logD)"
        ))

activation_formula_str <-
paste(
    "(binding_mask * (",
    activation_bound,
    ")) + ((1 - binding_mask) * (",
    activation_open,
    "))",
    sep = ""
    ) %>%
    str_replace_all(c(
        "Fb" = "concentration",
        "Kb" = "exp(logKb)",
        "L" = "exp(logL)",
        "E" = "exp(logE)"
        ))

combined_formula <-
paste(
    "response ~ (inhibition_mask * (",
    inhibition_formula_str,
    ")) + ((1 - inhibition_mask) * (",
    activation_formula_str,
    "))",
    sep = ""
    ) %>%
    as.formula()

full_formula <-
    bf(
        combined_formula,
        nl = TRUE,
        logKa + logKb + logD + logE ~ 1 + (1||unique_experiment_id),
        logL ~ 1 + method + (method||unique_experiment_id),
        family = gaussian()
        )

reduced_formula <-
    bf(
        combined_formula,
        nl = TRUE,
        logKa + logKb + logD + logE ~ 1,
        logL ~ 1 + method,
        family = gaussian()
        )

full_priors <-
    c(
        #posterior probability distribution of ec50s from pcf fluorescence
        set_prior("normal(log(1e4), log(5))", nlpar = "logKa", class = "b"),
        set_prior("normal(log(1e4), log(5))", nlpar = "logKb", class = "b"),
        #99% of density within popens of 0.01 and 0.99
        set_prior("normal(log(1), log(3))", nlpar = "logL", class = "b"),
        #test D prior
        set_prior("normal(log(1e-2), log(3))", nlpar = "logD", class = "b", ub = 0),
        #test E prior
        set_prior("normal(log(1e2), log(3))", nlpar = "logE", class = "b", lb = 0),
        #standard cauchy prior for sigmas
        set_prior("cauchy(0, 1)", nlpar = "logKa", class = "sd"),
        set_prior("cauchy(0, 1)", nlpar = "logKb", class = "sd"),
        set_prior("cauchy(0, 1)", nlpar = "logL", class = "sd"),
        set_prior("cauchy(0, 1)", nlpar = "logD", class = "sd"),
        set_prior("cauchy(0, 1)", nlpar = "logE", class = "sd")
        )

brms_iter <- 4000
brms_warmup <- 2000
brms_chains <- 4
brms_thin <- 1
brms_seed <- 2021

brm(
    formula = full_formula,
    prior = full_priors,
    data = concresp,
    chains = brms_chains,
    iter = brms_iter,
    warmup = brms_warmup,
    thin  = brms_thin,
    seed = brms_seed,
    control = list(adapt_delta = 0.95, max_treedepth = 10),
    cores = getOption("mc.cores", 4),
    file = "data/mwc_fits_new_model/activation_combined_fit_bounded_corrected_220521",
    sample_prior = "yes",
    save_all_pars = TRUE
    ) -> test_run

concresp %>%
ungroup() %>%
expand(nesting(method, measure, inhibition_mask, binding_mask), concentration = 10^seq(-8, -2, length.out = 51)) %>%
add_fitted_draws(test_run, re_formula = reduced_formula) %>%
median_qi(.value, .width = .95) -> inferred_underlying

ggplot() +
geom_ribbon(data = inferred_underlying, aes(x = concentration, ymin = .lower, ymax = .upper, colour = interaction(method, measure), fill = interaction(method, measure)), alpha = 0.5) +
geom_quasirandom(data = concresp, aes(x = concentration, y = response, fill = interaction(method, measure)), size = 3, shape = 21, width=0.1) +
scale_colour_brewer(aesthetics = c("colour", "fill"), palette = "Pastel1") +
scale_x_log10() +
facet_wrap(vars(inhibition_mask)) +
coord_cartesian(ylim = c(-0.1, 1.2))

test_run %>%
gather_draws(`b_.*`, regex = TRUE) %>%
separate(.variable, into = c("discard", ".variable", "type"), sep = "_") %>%
mutate(
    .value = case_when(.variable == "logE" ~ log(1/exp(.value)), TRUE ~ .value),
    measure = case_when(.variable %in% c("logD", "logKa") ~ "inhibition", .variable %in% c("logE", "logKb") ~ "activation"),
    .variable = case_when(.variable == "logE" ~ "logD", .variable == "logKb" ~ "logKa", TRUE ~ .variable)
    ) %>%
filter(type == "Intercept") -> draws

ggplot() +
stat_slab(
    data = draws,
    position = position_dodge(0.2),
    aes(y = .variable, x = exp(.value), fill = measure, fill_ramp = stat(cut_cdf_qi(cdf, .width = c(.5, .8, .95), labels = scales::percent_format())))
    ) +
scale_fill_ramp_discrete(range = c(1, 0.2), na.translate = FALSE) +
labs(fill_ramp = "Interval") +
facet_wrap(vars(.variable), scales = "free") +
scale_x_log10()

post <- posterior_samples(test_run, add_chain = T)

mcmc_hex(post)