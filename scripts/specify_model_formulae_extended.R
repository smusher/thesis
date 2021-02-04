library(tidyverse)
library(brms)
library(tidybayes)
library(modelr)
source("scripts/specify_model_formulae.R")

mwc_model <- function(Fa, Ka, Da, Fb, Kb, Db, L, C){

	# closed states, each monomer bound to one ligand only
	s4u <- 1
	
	s3u1a <- 4*Ka*Fa*s4u
	s2u2a <- 3/2*Ka*Fa*s3u1a
	s1u3a <- 2/3*Ka*Fa*s2u2a
	s4a <- 1/4*Ka*Fa*s1u3a
	
	s3u1b <- 4*Kb*Fb*s4u
	s2u2b <- 3/2*Kb*Fb*s3u1b
	s1u3b <- 2/3*Kb*Fb*s2u2b
	s4b <- 1/4*Kb*Fb*s1u3b
	
	s2u1a1b <- 3*Kb*Fb*s3u1a
	s1u1a2b <- 1*Kb*Fb*s2u1a1b
	s1a3b <- 1/3*Kb*Fb*s1u1a2b
	
	s1u2a1b <- 1*Ka*Fa*s2u1a1b
	s2a2b <- 1/2*Kb*Fb*s1u2a1b
	
	s3a1b <- 1/3*Ka*Fa*s1u2a1b
	
	#closed states, monomers bound to multiple ligands
	s3u1ab <- C*Kb*Fb*s3u1a
	
	s2u1a1ab <- 3*Ka*Fa*s3u1ab
	s1u2a1ab <- Ka*Fa*s2u1a1ab
	s3a1ab <- 1/3*Ka*Fa*s1u2a1ab
	
	s2u1b1ab <- 3*Kb*Fb*s3u1ab
	s1u2b1ab <- Kb*Fb*s2u1b1ab
	s3b1ab <- 1/3*Kb*Fb*s1u2b1ab
	
	s1u1a1b1ab <- 2*Ka*Fa*s2u1b1ab
	s2a1b1ab <- 1/2*Ka*Fa*s1u1a1b1ab
	s1a2b1ab <- 1/2*Kb*Fb*s1u1a1b1ab
	
	s2u2ab <- 1/2*C*Kb*Fb*s2u1a1ab
	s1u1a2ab <- 2*Ka*Fa*s2u2ab
	s2a2ab <- 1/2*Ka*Fa*s1u1a2ab
	
	s1a1b2ab <- Kb*Fb*s1u1a2ab
	s1u1b2ab <- 2*Kb*Fb*s2u2ab
	s2b2ab <- 1/2*Kb*Fb*s1u1b2ab
	
	s1u3ab <- 1/3*C*Kb*Fb*s1u1a2ab
	s1a3ab <- Ka*Fa*s1u3ab
	s1b3ab <- Kb*Fb*s1u3ab
	s4ab <- 1/4*C*Kb*Fb*s1a3ab
	
	#open states, each monomer bound to one ligand only
	s4u_o <- L
	
	s3u1a_o <- Da*L*s3u1a
	s2u2a_o <- Da^2*L*s2u2a
	s1u3a_o <- Da^3*L*s1u3a
	s4a_o <- Da^4*L*s4a
	
	s3u1b_o <- Db*L*s3u1b
	s2u2b_o <- Db^2*L*s2u2b
	s1u3b_o <- Db^3*L*s1u3b
	s4b_o <- Db^4*L*s4b
	
	s2u1a1b_o <- Da*Db*L*s2u1a1b
	s1u1a2b_o <- Da*Db^2*L*s1u1a2b
	s1a3b_o <- Da*Db^3*L*s1a3b
	
	s1u2a1b_o <- Da^2*Db*L*s1u2a1b
	s2a2b_o <- Da^2*Db^2*L*s2a2b
	
	s3a1b_o <- Da^3*Db*L*s3a1b
	
	#open states, monomers bound to multiple ligands
	s3u1ab_o <- Da*Db*L*s3u1ab
	
	s2u1a1ab_o <- Da^2*Db*L*s2u1a1ab
	s1u2a1ab_o <- Da^3*Db*L*s1u2a1ab
	s3a1ab_o <- Da^4*Db*L*s3a1ab
	
	s2u1b1ab_o <- Da*Db^2*L*s2u1b1ab
	s1u2b1ab_o <- Da*Db^3*L*s1u2b1ab
	s3b1ab_o <- Da*Db^3*L*s3b1ab
	
	s1u1a1b1ab_o <- Da^2*Db^2*L*s1u1a1b1ab
	s2a1b1ab_o <- Da^3*Db^2*L*s2a1b1ab
	s1a2b1ab_o <- Da^2*Db^3*L*s1a2b1ab
	
	s2u2ab_o <- Da^2*Db^2*L*s2u2ab
	s1u1a2ab_o <- Da^3*Db^2*L*s1u1a2ab
	s2a2ab_o <- Da^4*Db^2*L*s2a2ab
	
	s1a1b2ab_o <- Da^3*Db^3*L*s1a1b2ab
	s1u1b2ab_o <- Da^2*Db^3*L*s1u1b2ab
	s2b2ab_o <- Da^2*Db^4*L*s2b2ab
	
	s1u3ab_o <- Da^3*Db^3*L*s1u3ab
	s1a3ab_o <- Da^4*Db^3*L*s1a3ab
	s1b3ab_o <- Da^3*Db^4*L*s1b3ab
	s4ab_o <- Da^4*Db^4*L*s4ab

	closed_states <-
		c(
		s4u, s3u1a, s2u2a, s1u3a, s4a, s3u1b, s2u2b, s1u3b, s4b,
		s2u1a1b, s1u2a1b, s1u1a2b, s3a1b, s1a3b, s2a2b,
		s3u1ab, s2u1a1ab, s1u2a1ab, s3a1ab, s2u1b1ab, s1u2b1ab,
		s3b1ab, s1u1a1b1ab, s2a1b1ab, s1a2b1ab, s2u2ab, s1u1a2ab,
		s2a2ab, s1a1b2ab, s1u1b2ab, s2b2ab, s1u3ab, s1a3ab, s1b3ab, s4ab
			)

	open_states <-
		c(
		s4u_o, s3u1a_o, s2u2a_o, s1u3a_o, s4a_o, s3u1b_o, s2u2b_o, s1u3b_o, s4b_o,
		s2u1a1b_o, s1u2a1b_o, s1u1a2b_o, s3a1b_o, s1a3b_o, s2a2b_o,
		s3u1ab_o, s2u1a1ab_o, s1u2a1ab_o, s3a1ab_o, s2u1b1ab_o, s1u2b1ab_o,
		s3b1ab_o, s1u1a1b1ab_o, s2a1b1ab_o, s1a2b1ab_o, s2u2ab_o, s1u1a2ab_o,
		s2a2ab_o, s1a1b2ab_o, s1u1b2ab_o, s2b2ab_o, s1u3ab_o, s1a3ab_o, s1b3ab_o, s4ab_o
			)

	a_bound_states <-
		c(
		s3u1a, 2*s2u2a, 3*s1u3a, 4*s4a, s2u1a1b, 2*s1u2a1b, s1u1a2b, 3*s3a1b,
		s1a3b, 2*s2a2b, s3u1ab, 2*s2u1a1ab, 3*s1u2a1ab, 4*s3a1ab, s2u1b1ab, s1u2b1ab,
		s3b1ab, 2*s1u1a1b1ab, 3*s2a1b1ab, 2*s1a2b1ab, 2*s2u2ab, 3*s1u1a2ab,
		4*s2a2ab, 3*s1a1b2ab, 2*s1u1b2ab, 2*s2b2ab, 3*s1u3ab, 4*s1a3ab, 3*s1b3ab, 4*s4ab,
		s3u1a_o, 2*s2u2a_o, 3*s1u3a_o, 4*s4a_o, s2u1a1b_o, 2*s1u2a1b_o, s1u1a2b_o, 3*s3a1b_o,
		s1a3b_o, 2*s2a2b_o, s3u1ab_o, 2*s2u1a1ab_o, 3*s1u2a1ab_o, 4*s3a1ab_o, s2u1b1ab_o, s1u2b1ab_o,
		s3b1ab_o, 2*s1u1a1b1ab_o, 3*s2a1b1ab_o, 2*s1a2b1ab_o, 2*s2u2ab_o, 3*s1u1a2ab_o,
		4*s2a2ab_o, 3*s1a1b2ab_o, 2*s1u1b2ab_o, 2*s2b2ab_o, 3*s1u3ab_o, 4*s1a3ab_o, 3*s1b3ab_o, 4*s4ab_o
			)

	b_bound_states <-
		c(
		s3u1b, 2*s2u2b, 3*s1u3b, 4*s4b,
		s2u1a1b, s1u2a1b, 2*s1u1a2b, s3a1b, 3*s1a3b, 2*s2a2b,
		s3u1ab, s2u1a1ab, s1u2a1ab, s3a1ab, 2*s2u1b1ab, 3*s1u2b1ab,
		4*s3b1ab, 2*s1u1a1b1ab, 2*s2a1b1ab, 3*s1a2b1ab, 2*s2u2ab, 2*s1u1a2ab,
		2*s2a2ab, 3*s1a1b2ab, 3*s1u1b2ab, 4*s2b2ab, 3*s1u3ab, 3*s1a3ab, 4*s1b3ab, 4*s4ab,
		s3u1b_o, 2*s2u2b_o, 3*s1u3b_o, 4*s4b_o,
		s2u1a1b_o, s1u2a1b_o, 2*s1u1a2b_o, s3a1b_o, 3*s1a3b_o, 2*s2a2b_o,
		s3u1ab_o, s2u1a1ab_o, s1u2a1ab_o, s3a1ab_o, 2*s2u1b1ab_o, 3*s1u2b1ab_o,
		4*s3b1ab_o, 2*s1u1a1b1ab_o, 2*s2a1b1ab_o, 3*s1a2b1ab_o, 2*s2u2ab_o, 2*s1u1a2ab_o,
		2*s2a2ab_o, 3*s1a1b2ab_o, 3*s1u1b2ab_o, 4*s2b2ab_o, 3*s1u3ab_o, 3*s1a3ab_o, 4*s1b3ab_o, 4*s4ab_o
		)

	all_gating_states <- c(closed_states, open_states)
	all_binding_states <- 4 * all_gating_states

		tibble("open_fraction" = sum(open_states) / sum(all_gating_states),
		"a_bound_fraction" = sum(a_bound_states) / sum(all_binding_states),
		"b_bound_fraction" = sum(b_bound_states) / sum(all_binding_states)
		) %>% return()

}

atp_seq <- 10^seq(-7, -2, length.out = 31)

tribble(
	~data_gen_process, ~Fa,      ~Ka, ~Da, ~Fb, ~Kb, ~Db, ~L, ~C,
	"control_scheme_1", atp_seq, 1e4, 0.1, 0,    1,   1,   1,  1,
	"control_scheme_2", atp_seq, 1e4, 0.1, 1e-6, 1e5, 25, 0.01, 1,
	"control_scheme_3", atp_seq, 1e4, 0.1, 1e-6, 1e5, 25, 0.01, 0.1
	) %>%
unnest(Fa) %>%
rowwise() %>% mutate(res = purrr::map(Fa, mwc_model, Ka, Da, Fb, Kb, Db, L, C)) %>%
unnest(res) %>%
group_by(data_gen_process) %>%
mutate(open_fraction_norm = open_fraction/max(open_fraction)) %>%
pivot_longer(open_fraction:open_fraction_norm, names_to = "measure", values_to = "response") -> preview

ggplot(preview, aes(x = Fa, y = response, linetype = measure)) +
geom_line()+
scale_x_log10() +
facet_wrap(vars(data_gen_process))

tribble(
	~data_gen_process, ~Fa,                        ~Ka, ~Da, ~Fb,  ~Kb, ~Db, ~L,   ~C,
	"scheme_1_control", c(0,1e-6,1e-5,1e-4,1e-3),  1e4, 0.1, 0,    1,   1,   1,    1,
	"scheme_1_Ka_shift", c(0,1e-6,1e-5,1e-4,1e-3), 1e3, 0.1, 0,    1,   1,   1,    1,
	"scheme_1_L_shift", c(0,1e-6,1e-5,1e-4,1e-3),  1e4, 0.1, 0,    1,   1,   10,   1,
	"scheme_1_D_shift", c(0,1e-6,1e-5,1e-4,1e-3),  1e4, 0.8, 0,    1,   1,   1,    1,
	"scheme_2_control", c(0,1e-6,1e-5,1e-4,1e-3),  1e4, 0.1, 1e-6, 1e5, 25,  0.01, 1,
	"scheme_2_Kb_shift", c(0,1e-6,1e-5,1e-4,1e-3), 1e4, 0.1, 1e-6, 1e6, 25,  0.01, 1,
	"scheme_3a_control", c(0,1e-6,1e-5,1e-4,1e-3), 1e4, 0.1, 1e-6, 1e5, 25,  0.01, 0.1,
	"scheme_3a_Kb_shift",c(0,1e-6,1e-5,1e-4,1e-3), 1e4, 0.1, 1e-6, 1e6, 25,  0.01, 0.1,
	"scheme_3b_control", c(0,1e-6,1e-5,1e-4,1e-3), 1e4, 0.1, 1e-6, 1e5, 25,  0.01, 0.01,
	"scheme_3b_Kb_shift", c(0,1e-6,1e-5,1e-4,1e-3),1e4, 0.1, 1e-6, 1e6, 25,  0.01, 0.01,
	"scheme_3c_control", c(0,1e-6,1e-5,1e-4,1e-3), 1e4, 0.1, 1e-6, 1e5, 25,  0.01, 0.5,
	"scheme_3c_Kb_shift", c(0,1e-6,1e-5,1e-4,1e-3),1e4, 0.1, 1e-6, 1e6, 25,  0.01, 0.5

	) %>%
pivot_longer(Ka:C, names_to = "param", values_to = "value") %>%
rowwise() %>% mutate(draws = purrr::map(log(value), rlnorm, n=10, sd=0.25)) %>%
unnest(draws) %>%
select(-value) %>%
pivot_wider(names_from = param, values_from = draws) %>%
unnest(-Fa) %>%
mutate(experiment = row_number()) %>%
unnest(Fa) %>%
rowwise() %>% mutate(res = purrr::map(Fa, mwc_model, Ka, Da, Fb, Kb, Db, L, C)) %>%
unnest(res) %>%
group_by(data_gen_process, experiment) %>%
mutate(open_fraction_norm = open_fraction/max(open_fraction)) %>%
pivot_longer(open_fraction:open_fraction_norm, names_to = "measure", values_to = "response") -> test

ggplot(test %>% filter(measure %in% c("a_bound_fraction", "open_fraction_norm")), aes(x = Fa, y = response, fill = measure)) +
geom_point(shape = 21)+
scale_x_log10() +
facet_wrap(vars(data_gen_process))

test %>%
filter(measure %in% c("a_bound_fraction", "open_fraction_norm")) %>%
mutate(binding_mask = case_when(measure == "a_bound_fraction" ~ 1, measure == "open_fraction_norm" ~ 0)) %>%
rename(concentration = Fa) -> data_tofit

mwc_formula_transformed <-
    mwc_formula %>%
    str_replace_all("Ka", "(10^logKa)") %>%
    str_replace_all("L", "(10^logL)")

mwc_priors <-
    c(
        prior(uniform(0, 1), nlpar = "D", lb = 0, ub = 1),
        prior(uniform(2, 6), nlpar = "logKa", lb = 2, ub = 6),
        prior(normal(0, 0.7), nlpar = "logL")
    )

mwc_reduced_2 <-
    bf(
        as.formula(mwc_formula_transformed),
        logL ~ 0 + data_gen_process,
        D ~ 0 + data_gen_process,
        logKa ~ 0 + data_gen_process,
        nl = TRUE
    )

brm(
    formula = mwc_reduced_2,
    data = data_tofit,
    prior = mwc_priors,
    family = gaussian(),
    iter = 4000,
    warmup = 2000,
    chains = 4,
    thin = 1,
    seed = 2021,
    control = list(adapt_delta = 0.99, max_treedepth = 15),
    sample_prior = "yes",
    save_all_pars = TRUE,
    cores = getOption("mc.cores", 4),
    file = "data/mwc_fits_new_model/generated_data_seed2021.rds"
    ) -> test_1

data_tofit %>%
data_grid(
    concentration = 10^seq_range(-7:-2, n = 51),
    data_gen_process,
    binding_mask
    ) %>%
add_fitted_draws(test_1, re_formula = NA) %>%
group_by(data_gen_process, binding_mask, concentration) %>%
point_interval(.value, .point = median, .interval = qi, .width = .95) %>%
mutate(
    measure = case_when(binding_mask == 1 ~ "a_bound_fraction", TRUE ~ "open_fraction_norm"),
    .value = case_when(binding_mask == 1 ~ 1 - .value, TRUE ~ .value),
    .lower = case_when(binding_mask == 1 ~ 1 - .lower, TRUE ~ .lower),
    .upper = case_when(binding_mask == 1 ~ 1 - .upper, TRUE ~ .upper)
    ) -> fitted_draws

data_tofit <-
	data_tofit %>%
	mutate(
		response = case_when(binding_mask == 1 ~ 1 - response, TRUE ~ response)
		)

ggplot() +
geom_ribbon(data = fitted_draws, aes(x = concentration, ymin = .lower, ymax = .upper, fill = measure), alpha = 0.5) +
geom_line(data = fitted_draws, aes(x = concentration, y = .value, colour = measure)) +
geom_point(data = data_tofit, aes(x = concentration, y = response, fill = measure), shape = 21, size = 1) +
facet_wrap(vars(data_gen_process)) +
scale_y_continuous(breaks = c(0.0, 0.5, 1.0)) +
scale_x_log10() +
labs(
    x = "[TNP-ATP] (M)",
    y = expression(I/I[max] ~ or ~ F/F[max])
    ) +
scale_colour_brewer(aesthetics = c("colour", "fill"), palette = "Set1")

test_1 %>%
gather_draws(regex=TRUE, `b_.*`) %>%
separate(.variable, into = c("parameter", "data_gen_process"), sep = "_data_gen_process") -> param_fits

ggplot() +
stat_slab(
    data = param_fits,
    aes(
        x = data_gen_process,
        y = .value,
        fill = data_gen_process
        ),
    colour = "black", size = 0.5, normalize = "panels"
    ) +
theme(legend.position = "none") +
facet_wrap(vars(parameter), scales = "free") +
coord_flip() +
scale_colour_brewer(aesthetics = c("colour", "fill"), palette = "Set1")

