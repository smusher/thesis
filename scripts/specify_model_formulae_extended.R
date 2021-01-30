library(tidyverse)
library(patchwork)

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

Fb <- 0
Kb <- 1
Db <- 1
C <- 1

tibble("Fa" = 10^seq(-7, -2, length.out = 31)) %>%
mutate(res = purrr::map(Fa, mwc_model, 1e4, 0.1, Fb, Kb, Db, 1, C)) %>%
unnest(res) %>%
mutate(open_fraction_norm = open_fraction / max(open_fraction)) %>%
pivot_longer(-Fa, names_to = "measure", values_to = "fraction") -> test

ggplot(test, aes(x = Fa, y = fraction, linetype = measure)) +
geom_line() +
scale_x_log10() +
theme_thesis() +
ggtitle("Ka = 1e4, Da = 0.1, L = 1")

ggsave("/home/sam/thesis/figures/chx/simple_model_1.svg")

tibble("Fa" = 10^seq(-7, -2, length.out = 31)) %>%
mutate(res = purrr::map(Fa, mwc_model, 1e4, 0.1, Fb, Kb, Db, 10, C)) %>%
unnest(res) %>%
mutate(open_fraction_norm = open_fraction / max(open_fraction)) %>%
pivot_longer(-Fa, names_to = "measure", values_to = "fraction") -> test

ggplot(test, aes(x = Fa, y = fraction, linetype = measure)) +
geom_line() +
scale_x_log10() +
theme_thesis() +
ggtitle("Ka = 1e4, Da = 0.1, L = 10")

ggsave("/home/sam/thesis/figures/chx/simple_model_2.svg")

tibble("Fa" = 10^seq(-7, -2, length.out = 31)) %>%
mutate(res = purrr::map(Fa, mwc_model, 1e4, 0.1, Fb, Kb, Db, 100, C)) %>%
unnest(res) %>%
mutate(open_fraction_norm = open_fraction / max(open_fraction)) %>%
pivot_longer(-Fa, names_to = "measure", values_to = "fraction") -> test

ggplot(test, aes(x = Fa, y = fraction, linetype = measure)) +
geom_line() +
scale_x_log10() +
theme_thesis() +
ggtitle("Ka = 1e4, Da = 0.1, L = 100")

ggsave("/home/sam/thesis/figures/chx/simple_model_2b.svg")

tibble("Fa" = 10^seq(-7, -2, length.out = 31)) %>%
mutate(res = purrr::map(Fa, mwc_model, 1e4, 0.6, Fb, Kb, Db, 1, C)) %>%
unnest(res) %>%
mutate(open_fraction_norm = open_fraction / max(open_fraction)) %>%
pivot_longer(-Fa, names_to = "measure", values_to = "fraction") -> test

ggplot(test, aes(x = Fa, y = fraction, linetype = measure)) +
geom_line() +
scale_x_log10() +
theme_thesis() +
ggtitle("Ka = 1e4, Da = 0.6, L = 1")

ggsave("/home/sam/thesis/figures/chx/simple_model_3.svg")

tibble("Fa" = 10^seq(-7, -2, length.out = 31)) %>%
mutate(res = purrr::map(Fa, mwc_model, 1e4, 0.95, Fb, Kb, Db, 1, C)) %>%
unnest(res) %>%
mutate(open_fraction_norm = open_fraction / max(open_fraction)) %>%
pivot_longer(-Fa, names_to = "measure", values_to = "fraction") -> test

ggplot(test, aes(x = Fa, y = fraction, linetype = measure)) +
geom_line() +
scale_x_log10() +
theme_thesis() +
ggtitle("Ka = 1e4, Da = 0.95, L = 1")

ggsave("/home/sam/thesis/figures/chx/simple_model_3b.svg")

tibble("Fa" = 10^seq(-7, -2, length.out = 31)) %>%
mutate(res = purrr::map(Fa, mwc_model, Ka=1e4, Da=0.1, Fb=1e-6, Kb=1e5, Db=10, L=0.1, C=1)) %>%
unnest(res) %>%
mutate(open_fraction_norm = open_fraction / max(open_fraction)) %>%
pivot_longer(-Fa, names_to = "measure", values_to = "fraction") -> test

ggplot(test, aes(x = Fa, y = fraction, linetype = measure)) +
geom_line() +
scale_x_log10() +
theme_thesis() +
ggtitle("Ka = 1e4, Da = 0.1, L = 0.1, [b] = 1e-6, Kb = 1e5, Db = 10")

ggsave("/home/sam/thesis/figures/chx/pip_model_1.svg")

tibble("Fa" = 10^seq(-7, -2, length.out = 31)) %>%
mutate(res = purrr::map(Fa, mwc_model, Ka=1e4, Da=0.1, Fb=1e-6, Kb=1e6, Db=10, L=0.1, C=1)) %>%
unnest(res) %>%
mutate(open_fraction_norm = open_fraction / max(open_fraction)) %>%
pivot_longer(-Fa, names_to = "measure", values_to = "fraction") -> test

ggplot(test, aes(x = Fa, y = fraction, linetype = measure)) +
geom_line() +
scale_x_log10() +
theme_thesis() +
ggtitle("Ka = 1e4, Da = 0.1, L = 0.1, [b] = 1e-6, Kb = 1e6, Db = 10")

ggsave("/home/sam/thesis/figures/chx/pip_model_2.svg")

tibble("Fa" = 10^seq(-7, -2, length.out = 31)) %>%
mutate(res = purrr::map(Fa, mwc_model, Ka=1e4, Da=0.1, Fb=1e-6, Kb=1e5, Db=50, L=0.1, C=1)) %>%
unnest(res) %>%
mutate(open_fraction_norm = open_fraction / max(open_fraction)) %>%
pivot_longer(-Fa, names_to = "measure", values_to = "fraction") -> test

ggplot(test, aes(x = Fa, y = fraction, linetype = measure)) +
geom_line() +
scale_x_log10() +
theme_thesis() +
ggtitle("Ka = 1e4, Da = 0.1, L = 0.1, [b] = 1e-6, Kb = 1e5, Db = 50")

ggsave("/home/sam/thesis/figures/chx/pip_model_3.svg")

tibble("Fa" = 10^seq(-7, -2, length.out = 31)) %>%
mutate(res = purrr::map(Fa, mwc_model, Ka=1e4, Da=0.1, Fb=1e-6, Kb=1e5, Db=10, L=0.1, C=0.5)) %>%
unnest(res) %>%
mutate(open_fraction_norm = open_fraction / max(open_fraction)) %>%
pivot_longer(-Fa, names_to = "measure", values_to = "fraction") -> test

ggplot(test, aes(x = Fa, y = fraction, linetype = measure)) +
geom_line() +
scale_x_log10() +
theme_thesis() +
ggtitle("Ka = 1e4, Da = 0.1, L = 0.1, [b] = 1e-6, Kb = 1e5, Db = 10, C = 0.5")

ggsave("/home/sam/thesis/figures/chx/full_model_1.svg")

tibble("Fa" = 10^seq(-7, -2, length.out = 31)) %>%
mutate(res = purrr::map(Fa, mwc_model, Ka=1e4, Da=0.1, Fb=1e-6, Kb=1e5, Db=10, L=0.1, C=0.01)) %>%
unnest(res) %>%
mutate(open_fraction_norm = open_fraction / max(open_fraction)) %>%
pivot_longer(-Fa, names_to = "measure", values_to = "fraction") -> test

ggplot(test, aes(x = Fa, y = fraction, linetype = measure)) +
geom_line() +
scale_x_log10() +
theme_thesis() +
ggtitle("Ka = 1e4, Da = 0.1, L = 0.1, [b] = 1e-6, Kb = 1e5, Db = 10, C = 0.01")

ggsave("/home/sam/thesis/figures/chx/full_model_2.svg")

tibble("Fa" = 10^seq(-7, -2, length.out = 31)) %>%
mutate(res = purrr::map(Fa, mwc_model, Ka=1e4, Da=0.1, Fb=1e-6, Kb=1e6, Db=10, L=0.1, C=0.01)) %>%
unnest(res) %>%
mutate(open_fraction_norm = open_fraction / max(open_fraction)) %>%
pivot_longer(-Fa, names_to = "measure", values_to = "fraction") -> test

ggplot(test, aes(x = Fa, y = fraction, linetype = measure)) +
geom_line() +
scale_x_log10() +
theme_thesis() +
ggtitle("Ka = 1e4, Da = 0.1, L = 0.1, [b] = 1e-6, Kb = 1e6, Db = 10, C = 0.01")

ggsave("/home/sam/thesis/figures/chx/full_model_3.svg")