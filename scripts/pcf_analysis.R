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
    mutate(
        response = case_when(measure == "fluorescence" ~ 1 - response, TRUE ~ response),
        binding_mask = case_when(measure == "fluorescence" ~ 1, TRUE ~ 0)
    ) %>%
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
    	)

control_data <-
	concresp %>%
	filter(construct == "W311*-GFP+SUR")

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

mwc_brms_formula <-
	bf(
		mwc_full_formula,
		nl = TRUE,
		logKa + logDa + logL ~ 0 + construct + (construct||unique_experiment_id),
		sigma ~ (1||construct),
		family = gaussian()
		)

mwc_brms_formula_restricted <-
	bf(
		mwc_full_formula,
		nl = TRUE,
		logKa + logDa ~ 0 + construct,
		logL ~ 0 + construct + (construct||unique_experiment_id),
		sigma ~ (1||construct),
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

brms_iter <- 2000
brms_warmup <- 1000
brms_chains <- 4
brms_thin <- 1
brms_seed <- 2021

brm(
	formula = mwc_brms_formula,
	prior = mwc_brms_priors,
	data = concresp,
	chains = brms_chains,
	iter = brms_iter,
	warmup = brms_warmup,
	thin  = brms_thin,
	seed = brms_seed,
	control = list(adapt_delta = 0.99, max_treedepth = 15),
	cores = getOption("mc.cores", 4),
	file = "data/mwc_fits_new_model/full_fit_030321_short",
	sample_prior = "yes",
	save_all_pars = TRUE
	) %>%
	add_criterion("loo", moment_match=TRUE) -> test_run

brm(
	formula = mwc_brms_formula_restricted,
	prior = mwc_brms_priors_restricted,
	data = concresp,
	chains = brms_chains,
	iter = brms_iter,
	warmup = brms_warmup,
	thin  = brms_thin,
	seed = brms_seed,
	control = list(adapt_delta = 0.99, max_treedepth = 15),
	cores = getOption("mc.cores", 4),
	file = "data/mwc_fits_new_model/full_fit_030321_short_restricted",
	sample_prior = "yes",
	save_all_pars = TRUE
	) %>%
	add_criterion("loo", moment_match=TRUE) -> test_run_2

posterior <- as.array(test_run_2)
dim(posterior)
dimnames(posterior)

mcmc_areas(
  posterior, 
  pars = c(
  	"b_logKa_Intercept", "b_logDa_Intercept", "b_logL_Intercept", 
  	"prior_b_logKa", "prior_b_logDa", "prior_b_logL"),
  prob = 0.8, # 80% intervals
  prob_outer = 0.99, # 99%
  point_est = "median"
)

mcmc_pairs(
	posterior,
	pars = c("b_logKa_Intercept", "b_logDa_Intercept", "b_logL_Intercept"),
    diag_fun="dens",
    off_diag_args = list(size = 0.75, alpha = 0.2)
    )

mcmc_pairs(
	posterior,
	pars = c("b_logKa_Intercept", "b_logDa_Intercept", "b_logL_Intercept"),
    diag_fun="dens",
    off_diag_args = list(size = 0.75, alpha = 0.2),
	transform = list(
		b_logKa_Intercept = exp,
		b_logDa_Intercept = exp,
		b_logL_Intercept = exp
		)
    )

control_data %>%
ungroup() %>%
expand(nesting(measure, binding_mask), concentration = 10^seq(-8, -2, length.out = 51)) %>%
add_fitted_draws(test_run, re_formula = NA) %>%
median_qi(.value, .width = .95) %>%
mutate(model = "full") -> inferred_underlying_1

control_data %>%
ungroup() %>%
expand(nesting(measure, binding_mask), concentration = 10^seq(-8, -2, length.out = 51)) %>%
add_fitted_draws(test_run_2, re_formula = NA) %>%
median_qi(.value, .width = .95) %>%
mutate(model = "restricted") -> inferred_underlying_2

inferred_underlying <-
	bind_rows(inferred_underlying_1, inferred_underlying_2)

ggplot() +
geom_ribbon(data = inferred_underlying, aes(x = concentration, ymin = .lower, ymax = .upper, colour = measure, linetype = model), fill = NA) +
geom_quasirandom(data = control_data, aes(x = concentration, y = response, fill = measure), size = 3, shape = 21, width=0.1) +
scale_colour_brewer(aesthetics = c("colour", "fill"), palette = "Pastel1") +
scale_x_log10()

control_data %>%
ungroup() %>%
expand(nesting(measure, binding_mask, unique_experiment_id), concentration = 10^seq(-8, -2, length.out = 51)) %>%
add_fitted_draws(test_run) %>%
median_qi(.value, .width = .95) %>%
mutate(model = "full") -> per_experiment_1

control_data %>%
ungroup() %>%
expand(nesting(measure, binding_mask, unique_experiment_id), concentration = 10^seq(-8, -2, length.out = 51)) %>%
add_fitted_draws(test_run_2) %>%
median_qi(.value, .width = .95) %>%
mutate(model = "restricted") -> per_experiment_2

per_experiment <-
	bind_rows(per_experiment_1, per_experiment_2)

ggplot() +
geom_ribbon(data = per_experiment, aes(x = concentration, ymin = .lower, ymax = .upper, colour = measure, linetype = model), fill = NA) +
geom_point(data = control_data, aes(x = concentration, y = response, fill = measure), size = 3, shape = 21) +
scale_colour_brewer(aesthetics = c("colour", "fill"), palette = "Pastel1") +
scale_x_log10() +
facet_wrap(vars(unique_experiment_id))

test_run %>%
gather_draws(`b_.*`, `sd_.*`, regex = TRUE) %>%
mutate(model = "full") -> draws_1

test_run_2 %>%
gather_draws(`b_.*`, `sd_.*`, regex = TRUE) %>%
mutate(model = "restricted") -> draws_2

draws <-
	bind_rows(draws_1, draws_2) %>%
	mutate(.value = case_when(
		.variable == "b_logL_Intercept" ~ exp(.value)/(1+exp(.value)),
		TRUE ~ .value)
	)

ggplot() +
stat_slab(
    data = draws,
    aes(y = model, x = .value, fill = model, fill_ramp = stat(cut_cdf_qi(cdf, .width = c(.5, .8, .95), labels = scales::percent_format())))
    ) +
scale_fill_brewer(palette = "Pastel1", aesthetics = c("fill", "colour")) +
scale_fill_ramp_discrete(range = c(1, 0.2), na.translate = FALSE) +
labs(fill_ramp = "Interval") +
facet_wrap(vars(.variable), scales = "free")
