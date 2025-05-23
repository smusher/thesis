library(tidyverse)
library(brms)
library(tidybayes)
library(bayesplot)
library(ggdist)
library(distributional)
source("scripts/ch4/mwc_modelling_functions.R")

concresp <-
    read_csv("data/combined_drc_data.csv") %>%
    filter(concentration > 0, method == "pcf") %>%
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
    mutate(
        response = case_when(measure == "fluorescence" ~ (1 - (log2(response + 1))) / 0.9, TRUE ~ response),
        binding_mask = case_when(measure == "fluorescence" ~ 1, TRUE ~ 0),
        log_concentration = log10(concentration)
    )

ggplot() +
geom_point(data = concresp, aes(x = concentration, y = response, colour = measure)) +
scale_color_viridis_d() +
facet_wrap(vars(construct)) +
scale_x_log10()

control_data <-
	concresp %>%
	filter(construct %in% c("W311*-GFP+SUR", "W311*,C166S-GFP+SUR"))

mwc_full_formula <-
	mwc_formula_string %>%
	str_replace_all(c(
		"1/2" = "0.5",
		"1/3" = "0.333",
		"2/3" = "0.666",
		"3/2" = "1.5",
		"1/4" = "0.25",
		"Fa" = "concentration",
		"Ka" = "exp(logKa)",
		"L" = "exp(logL)",
		"Da" = "exp(logDa)",
		"Fb" = "0",
		"Kb" = "1",
		"Db" = "1",
		"C" = "1"
		)) %>%
	as.formula()

single_full_formula <-
	single_formula_string %>%
	str_replace_all(c(
		"1/2" = "0.5",
		"1/3" = "0.333",
		"2/3" = "0.666",
		"3/2" = "1.5",
		"1/4" = "0.25",
		"Fa" = "concentration",
		"Ka" = "exp(logKa)",
		"L" = "exp(logL)",
		"Da" = "exp(logDa)",
		"Fb" = "0",
		"Kb" = "1",
		"Db" = "1",
		"C" = "1"
		)) %>%
	as.formula()

mwc_brms_formula <-
	bf(
		mwc_full_formula,
		nl = TRUE,
		logKa + logDa + logL ~ 0 + construct + (construct||unique_experiment_id),
		sigma ~ (1||construct:binding_mask),
		family = gaussian()
		)

mwc_brms_formula_restricted <-
	bf(
		mwc_full_formula,
		nl = TRUE,
		logKa + logDa ~ 0 + construct,
		logL ~ 0 + construct + (construct||unique_experiment_id),
		sigma ~ (1||construct:binding_mask),
		family = gaussian()
		)

single_brms_formula <-
	bf(
		single_full_formula,
		nl = TRUE,
		logKa + logDa + logL ~ 0 + construct + (construct||unique_experiment_id),
		sigma ~ (1||construct:binding_mask),
		family = gaussian()
		)

mwc_brms_priors <-
	c(
		#posterior probability distribution of ec50s from pcf fluorescence
		set_prior("normal(log(1e4), log(5))", nlpar = "logKa", class = "b"),
		#99% of density within popens of 0.01 and 0.99
		set_prior("normal(log(1), log(5))", nlpar = "logL", class = "b"),
		#test D prior
		set_prior("normal(log(0.1), log(5))", nlpar = "logDa", class = "b"),
		#standard cauchy prior for sigmas
		set_prior("cauchy(0, 1)", nlpar = "logKa", class = "sd"),
		set_prior("cauchy(0, 1)", nlpar = "logL", class = "sd"),
		set_prior("cauchy(0, 1)", nlpar = "logDa", class = "sd"),
		set_prior("cauchy(0, 1)", dpar = "sigma", class = "Intercept"),
		set_prior("cauchy(0, 1)", dpar = "sigma", class = "sd")
		)

mwc_brms_priors_restricted <-
	c(
		#posterior probability distribution of ec50s from pcf fluorescence
		set_prior("normal(log(1e4), log(5))", nlpar = "logKa", class = "b"),
		#99% of density within popens of 0.01 and 0.99
		set_prior("normal(log(1), log(5))", nlpar = "logL", class = "b"),
		#test D prior
		set_prior("normal(log(0.1), log(5))", nlpar = "logDa", class = "b"),
		#standard cauchy prior for sigmas
		set_prior("cauchy(0, 1)", nlpar = "logL", class = "sd"),
		set_prior("cauchy(0, 1)", dpar = "sigma", class = "Intercept"),
		set_prior("cauchy(0, 1)", dpar = "sigma", class = "sd")
		)

brms_iter <- 11000
brms_warmup <- 1000
brms_chains <- 4
brms_thin <- 1
brms_seed <- 2021

brm(
	formula = mwc_brms_formula,
	prior = mwc_brms_priors,
	data = control_data,
	chains = brms_chains,
	iter = brms_iter,
	warmup = brms_warmup,
	thin  = brms_thin,
	seed = brms_seed,
	control = list(adapt_delta = 0.95, max_treedepth = 10),
	cores = getOption("mc.cores", 4),
	file = "data/mwc_fits_new_model/control_mwc_fit_270321",
	sample_prior = "yes",
	save_all_pars = TRUE
	) -> control_run

brm(
	formula = single_brms_formula,
	prior = mwc_brms_priors,
	data = control_data,
	chains = brms_chains,
	iter = brms_iter,
	warmup = brms_warmup,
	thin  = brms_thin,
	seed = brms_seed,
	control = list(adapt_delta = 0.95, max_treedepth = 10),
	cores = getOption("mc.cores", 4),
	file = "data/mwc_fits_new_model/control_single_fit_270321",
	sample_prior = "yes",
	save_all_pars = TRUE
	) -> control_run_2

brm(
	formula = mwc_brms_formula,
	prior = mwc_brms_priors,
	data = concresp,
	chains = brms_chains,
	iter = brms_iter,
	warmup = brms_warmup,
	thin  = brms_thin,
	seed = brms_seed,
	control = list(adapt_delta = 0.95, max_treedepth = 10),
	cores = getOption("mc.cores", 4),
	file = "data/mwc_fits_new_model/full_fit_140321_short",
	sample_prior = "yes",
	save_all_pars = TRUE
	) -> test_run

brm(
	formula = mwc_brms_formula_restricted,
	prior = mwc_brms_priors_restricted,
	data = concresp,
	chains = brms_chains,
	iter = brms_iter,
	warmup = brms_warmup,
	thin  = brms_thin,
	seed = brms_seed,
	control = list(adapt_delta = 0.95, max_treedepth = 10),
	cores = getOption("mc.cores", 4),
	file = "data/mwc_fits_new_model/full_fit_110321_short_restricted",
	sample_prior = "yes",
	save_all_pars = TRUE
	) -> test_run_2

brm(
	formula = single_brms_formula,
	prior = mwc_brms_priors,
	data = concresp,
	chains = brms_chains,
	iter = brms_iter,
	warmup = brms_warmup,
	thin  = brms_thin,
	seed = brms_seed,
	control = list(adapt_delta = 0.95, max_treedepth = 10),
	cores = getOption("mc.cores", 4),
	file = "data/mwc_fits_new_model/single_fit_250321_short",
	sample_prior = "yes",
	save_all_pars = TRUE
	) -> test_run_3

concresp %>%
mutate(response = case_when(measure == "fluorescence" ~ 1 - response, TRUE ~ response)) -> concresp_2

concresp %>%
ungroup() %>%
expand(nesting(construct, measure, binding_mask), concentration = 10^seq(-8, -2, length.out = 51)) %>%
add_fitted_draws(test_run, re_formula = NA) %>%
mutate(.value = case_when(measure == "fluorescence" ~ 1 - .value, TRUE ~ .value)) %>%
median_qi(.value, .width = .95) -> inferred_underlying

ggplot() +
geom_ribbon(data = inferred_underlying %>% filter(construct != "W311*-GFP+SUR+MG"), aes(x = concentration, ymin = .lower, ymax = .upper, fill = measure), alpha = 0.5) +
geom_line(data = inferred_underlying %>% filter(construct != "W311*-GFP+SUR+MG"), aes(x = concentration, y= .value, colour = measure)) +
geom_quasirandom(data = concresp_2 %>% filter(construct != "W311*-GFP+SUR+MG"), aes(x = concentration, y = response, fill = measure), size = 3, shape = 21, width=0.1) +
scale_colour_brewer(aesthetics = c("colour", "fill"), palette = "Pastel1") +
scale_x_log10() +
facet_wrap(vars(construct)) +
coord_cartesian(ylim = c(-0.1, 1.2))

concresp %>%
ungroup() %>%
expand(nesting(construct, measure, binding_mask, unique_experiment_id), concentration = 10^seq(-8, -2, length.out = 51)) %>%
add_fitted_draws(test_run, allow_new_levels=TRUE) %>%
mutate(.value = case_when(measure == "fluorescence" ~ 1 - .value, TRUE ~ .value)) %>%
median_qi(.value, .width = .95) -> per_experiment

ggplot() +
geom_ribbon(data = per_experiment %>% filter(construct != "W311*-GFP+SUR+MG"), aes(x = concentration, ymin = .lower, ymax = .upper, colour = measure, fill = measure), alpha = 0.5) +
geom_point(data = concresp_2 %>% filter(construct != "W311*-GFP+SUR+MG"), aes(x = concentration, y = response, fill = measure), size = 3, shape = 21) +
scale_colour_brewer(aesthetics = c("colour", "fill"), palette = "Pastel1") +
scale_x_log10() +
facet_wrap(vars(construct, unique_experiment_id)) +
coord_cartesian(ylim = c(-0.1, 1.2))

test_run %>%
gather_draws(`b_.*`, regex = TRUE) %>%
separate(.variable, into = c(".variable", "construct"), sep = "_construct") %>%
group_by(construct, .variable) %>%
mutate(model = "full") -> draws_1

test_run_2 %>%
gather_draws(`b_.*`, regex = TRUE) %>%
separate(.variable, into = c(".variable", "construct"), sep = "_construct") %>%
group_by(construct, .variable) %>%
mutate(model = "restricted") -> draws_2

test_run_3 %>%
gather_draws(`b_.*`, regex = TRUE) %>%
separate(.variable, into = c(".variable", "construct"), sep = "_construct") %>%
group_by(construct, .variable) %>%
mutate(model = "single") -> draws_3

draws <- bind_rows(draws_1, draws_2, draws_3)

ggplot() +
stat_slab(
    data = draws %>% filter(!is.na(construct)),
    position = position_dodge(0.2),
    aes(y = construct, x = exp(.value), fill = model, fill_ramp = stat(cut_cdf_qi(cdf, .width = c(.5, .8, .95), labels = scales::percent_format())))
    ) +
scale_fill_ramp_discrete(range = c(1, 0.2), na.translate = FALSE) +
labs(fill_ramp = "Interval") +
facet_wrap(vars(.variable), scales = "free") +
scale_x_log10()

test_run %>%
gather_draws(`sd_.*`, regex = TRUE) %>%
separate(.variable, into = c(".variable", "construct"), sep = "_construct") %>%
group_by(construct, .variable) -> draws_2

ggplot() +
stat_slab(
    data = draws_2 %>% filter(!is.na(construct), construct != "__sigma_Intercept"),
    aes(y = construct, x = exp(.value), fill = construct, fill_ramp = stat(cut_cdf_qi(cdf, .width = c(.5, .8, .95), labels = scales::percent_format())))
    ) +
scale_fill_ramp_discrete(range = c(1, 0.2), na.translate = FALSE) +
labs(fill_ramp = "Interval") +
facet_wrap(vars(.variable), scales = "free") +
theme(legend.position = "none") +
scale_x_log10() +
coord_cartesian(xlim = c(1, 30))

prior_samples(test_run) %>%
as_tibble() %>%
pivot_longer(everything(), names_to = ".variable", values_to = ".value") %>%
mutate(construct = "prior") -> priors

ggplot() +
stat_slab(
    data = priors,
    aes(y = construct, x = exp(.value)), colour = "black", fill = NA, linetype = 2
    ) +
stat_slab(
    data = draws_1 %>% filter(construct == "W311MUMGFPPSUR"),
    aes(y = construct, x = exp(.value), fill = construct, fill_ramp = stat(cut_cdf_qi(cdf, .width = c(.5, .8, .95), labels = scales::percent_format())))
    ) +
scale_fill_ramp_discrete(range = c(1, 0.2), na.translate = FALSE) +
labs(fill_ramp = "Interval") +
facet_wrap(vars(.variable), scales = "free") +
theme(legend.position = "none") +
scale_x_log10()

control_data %>%
mutate(response = case_when(measure == "fluorescence" ~ 1 - response, TRUE ~ response)) -> control_data_2

control_data %>%
ungroup() %>%
expand(nesting(construct, measure, binding_mask), concentration = 10^seq(-8, -2, length.out = 51)) %>%
add_fitted_draws(control_run, re_formula = NA) %>%
mutate(.value = case_when(measure == "fluorescence" ~ 1 - .value, TRUE ~ .value)) %>%
median_qi(.value, .width = .95) -> inferred_underlying

control_data %>%
ungroup() %>%
expand(nesting(construct, measure, binding_mask), concentration = 10^seq(-8, -2, length.out = 51)) %>%
add_fitted_draws(control_run_2, re_formula = NA) %>%
mutate(.value = case_when(measure == "fluorescence" ~ 1 - .value, TRUE ~ .value)) %>%
median_qi(.value, .width = .95) -> inferred_underlying_2

ggplot() +
geom_ribbon(data = inferred_underlying_2, aes(x = concentration, ymin = .lower, ymax = .upper, colour = measure, fill = measure), alpha = 0.5) +
geom_quasirandom(data = control_data_2, aes(x = concentration, y = response, fill = measure), size = 3, shape = 21, width=0.1) +
scale_colour_brewer(aesthetics = c("colour", "fill"), palette = "Pastel1") +
scale_x_log10() +
facet_wrap(vars(construct)) +
coord_cartesian(ylim = c(-0.1, 1.2))

control_data %>%
ungroup() %>%
expand(nesting(construct, measure, binding_mask, unique_experiment_id), concentration = 10^seq(-8, -2, length.out = 51)) %>%
add_fitted_draws(control_run, allow_new_levels=TRUE) %>%
mutate(.value = case_when(measure == "fluorescence" ~ 1 - .value, TRUE ~ .value)) %>%
median_qi(.value, .width = .95) -> per_experiment

control_data %>%
ungroup() %>%
expand(nesting(construct, measure, binding_mask, unique_experiment_id), concentration = 10^seq(-8, -2, length.out = 51)) %>%
add_fitted_draws(control_run_2, allow_new_levels=TRUE) %>%
mutate(.value = case_when(measure == "fluorescence" ~ 1 - .value, TRUE ~ .value)) %>%
median_qi(.value, .width = .95) -> per_experiment_2

ggplot() +
geom_ribbon(data = per_experiment_2, aes(x = concentration, ymin = .lower, ymax = .upper, colour = measure, fill = measure), alpha = 0.5) +
geom_point(data = control_data_2, aes(x = concentration, y = response, fill = measure), size = 3, shape = 21) +
scale_colour_brewer(aesthetics = c("colour", "fill"), palette = "Pastel1") +
scale_x_log10() +
facet_wrap(vars(construct, unique_experiment_id)) +
coord_cartesian(ylim = c(-0.1, 1.2))

control_run %>%
gather_draws(`b_.*`, regex = TRUE) %>%
separate(.variable, into = c(".variable", "construct"), sep = "_construct") %>%
group_by(construct, .variable) %>%
mutate(model = "mwc") -> draws_1

control_run_2 %>%
gather_draws(`b_.*`, regex = TRUE) %>%
separate(.variable, into = c(".variable", "construct"), sep = "_construct") %>%
group_by(construct, .variable) %>%
mutate(model = "single") -> draws_2

draws <-
	bind_rows(draws_1, draws_2) %>%
	ungroup() %>%
	select(construct, model, .draw, .variable, .value) %>%
	filter(!is.na(construct)) %>%
mutate(
	.variable = str_replace_all(.variable, "b_log", ""),
	.value = exp(.value)
	) %>%
pivot_wider(names_from = .variable, values_from = .value) %>%
left_join(tibble(expand.grid(construct = c("W311MUMGFPPSUR","W311MUC166SMGFPPSUR") , Fa = 10^seq(-7, -2, length.out = 51))))

draws %>%
filter(model == "mwc") %>%
mutate(
	c0 = 1,
	c1 = 4*Ka*Fa*c0,
	c2 = 3/2*Ka*Fa*c1,
	c3 = 2/3*Ka*Fa*c2,
	c4 = 1/4*Ka*Fa*c3,
	o0 = L * c0,
	o1 = Da*L*c1,
	o2 = Da^2*L*c2,
	o3 = Da^3*L*c3,
	o4 = Da^4*L*c4
	) %>%
rowwise() %>%
mutate(
	all_states = sum(c_across(c0:o4)),
	across(c0:o4, ~ .x / all_states)
	) %>%
pivot_longer(c0:o4, names_to = "state", values_to = "weight") %>%
group_by(construct, model, state, Fa) %>%
median_qi(weight, .width = .95)  -> sw_1

draws %>%
filter(model == "single") %>%
mutate(
	c0 = 1,
	c1 = 4*Ka*Fa*c0,
	c2 = 3/2*Ka*Fa*c1,
	c3 = 2/3*Ka*Fa*c2,
	c4 = 1/4*Ka*Fa*c3,
	o0 = L * c0,
	o1 = Da*L*c1,
	o2 = Da*L*c2,
	o3 = Da*L*c3,
	o4 = Da*L*c4
	) %>%
rowwise() %>%
mutate(
	all_states = sum(c_across(c0:o4)),
	across(c0:o4, ~ .x / all_states)
	) %>%
pivot_longer(c0:o4, names_to = "state", values_to = "weight") %>%
group_by(construct, model, state, Fa) %>%
median_qi(weight, .width = .95) -> sw_2

bind_rows(sw_1, sw_2) %>%
mutate(
	conformation = case_when(str_detect(state, "o") == TRUE ~ "open", TRUE ~ "closed"),
	state = as.numeric(str_replace_all(state, c("o" = "", "c" = "")))
	) -> states_weights

ggplot(states_weights %>% filter(construct == "W311MUC166SMGFPPSUR"), aes(x = Fa, ymin = .lower, ymax = .upper, colour = factor(state), fill = factor(state))) +
geom_ribbon(alpha = 0.5) +
facet_grid(rows = vars(conformation), cols = vars(model), scales = "free") +
scale_fill_brewer(palette = "Blues", aesthetics = c("colour", "fill")) +
scale_x_log10() +
scale_y_log10() +
coord_cartesian(ylim = c(1e-4, 1))

