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
    ungroup()

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
		logKa + logDa + logL ~ 1 + (1||unique_experiment_id),
		sigma ~ (1||unique_experiment_id),
		family = gaussian()
		)

mwc_brms_formula_restricted <-
	bf(
		mwc_full_formula,
		nl = TRUE,
		logKa + logDa ~ 1,
		logL ~ 1 + (1||unique_experiment_id),
		sigma ~ (1||unique_experiment_id),
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

brms_iter <- 2000
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
	control = list(adapt_delta = 0.99, max_treedepth = 15),
	cores = getOption("mc.cores", 4),
	file = "data/mwc_fits_new_model/control_fit_230221_short",
	sample_prior = "yes",
	save_pars = save_pars(all = TRUE)
	) -> test_run

posterior <- as.array(test_run)
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
expand(nesting(measure, binding_mask), concentration = 10^seq(-8, -2, length.out = 31), unique_experiment_id = 1000) %>%
add_predicted_draws(test_run, allow_new_levels = TRUE) %>%
median_qi(.prediction, .width = c(.50, .80, .95)) -> predicted_future

control_data %>%
ungroup() %>%
expand(nesting(measure, binding_mask), concentration = 10^seq(-8, -2, length.out = 31)) %>%
add_fitted_draws(test_run, re_formula = NA) %>%
median_qi(.value, .width = c(.50, .80, .95)) -> inferred_underlying

ggplot() +
geom_ribbon(data = inferred_underlying %>% filter(.width == 0.95), aes(x = concentration, ymin = .lower, ymax = .upper, fill = measure, alpha = .width)) +
geom_point(data = control_data, aes(x = concentration, y = response, fill = measure), size = 2, shape = 21) +
scale_colour_brewer(aesthetics = c("colour", "fill"), palette = "Pastel1") +
scale_x_log10()