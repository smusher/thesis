library(tidyverse)
library(patchwork)

#closed states, each monomer bound to one ligand only
4u <- 1

3u1a <- 4*K1*F1*4u
2u2a <- 3/2*K1*F1*3u1a
1u3a <- 2/3*K1*F1*2u2a
4a <- 1/4*K1*F1*1u3a

3u1b <- 4*K2*F2*4u
2u2b <- 3/2*K2*F2*3u1b
1u3b <- 2/3*K2*F2*2u2b
4b <- 1/4*K2*F2*1u3b

2u1a1b <- 3*K2*F2*3u1a
1u1a2b <- 1*K2*F2*2u1a1b
1a3b <- 1/3*K2*F2*1u1a2b

1u2a1b <- 1*K1*F1*2u1a1b
2a2b <- 1/2*K2*F2*1u2a1b

3a1b <- 1/3*K1*F1*1u2a1b

#closed states, monomers bound to multiple ligands
3u1ab <- C*K2*F2*3u1a

2u1a1ab <- 3*K1*F1*3u1ab
1u2a1ab <- K1*F1*2u1a1ab
3a1ab <- 1/3*K1*F1*1u2a1ab

2u1b1ab<- 3*K2*F2*3u1ab
1u2b1ab <- K2*F2*2u1b1ab
3b1ab <- 1/3*K2*F2*1u2b1ab

1u1a1b1ab <- 2*K1*F1*2u1b1ab
2a1b1ab <- 1/2*K1*F1*1u1a1b1ab
1a2b1ab <- 1/2*K2*F2*1u1a1b1ab

2u2ab <- 1/2*C*K2*F2*2u1a1ab
1u1a2ab <- 2*K1*F1*2u2ab
2a2ab <- 1/2*K1*F1*1u1a2ab

1a1b2ab <- K2*F2*1u1a2ab
1u1b2ab <- 2*K2*F2*2u2ab
2b2ab <- 1/2*K2*F2*1u1b2ab

1u3ab <- 1/3*C*K2*F2*1u1a2ab
1a3ab <- K1*F1*1u3ab
1b3ab <- K2*F2*1u3ab
4ab <- 1/4*C*K2*F2*1a3ab

#open states, each monomer bound to one ligand only

4u_o <- L

all_states <-
    list(
        "(1)",
        "(4 * K1^1 * F1^1)",
        "(6 * K1^2 * F1^2)",
        "(4 * K1^3 * F1^3)",
        "(1 * K1^4 * F1^4)",
        "(L)",
        "(L * D1^1 * 4 * K1^1 * F1^1)",
        "(L * D1^2 * 6 * K1^2 * F1^2)",
        "(L * D1^3 * 4 * K1^3 * F1^3)",
        "(L * D1^4 * 1 * K1^4 * F1^4)",
        "(4 * K2^1 * F2^1)",
        "(6 * K2^2 * F2^2)",
        "(4 * K2^3 * F2^3)",
        "(1 * K2^4 * F2^4)",
        "(L * D2^1 * 4 * K2^1 * F2^1)",
        "(L * D2^2 * 6 * K2^2 * F2^2)",
        "(L * D2^3 * 4 * K2^3 * F2^3)",
        "(L * D2^4 * 1 * K2^4 * F2^4)",
        "(4 * K1^1 * F1^1 * K2^1 * F2^1 * C^1)",
        "(6 * K1^2 * F1^2 * K2^2 * F2^2 * C^2)",
        "(4 * K1^3 * F1^3 * K2^3 * F2^3 * C^3)",
        "(1 * K1^4 * F1^4 * K2^4 * F2^4 * C^4)",
        "(L * D1^1 * D2^1 * 4 * K1^1 * F1^1 * K2^1 * F2^1 * C^1)",
        "(L * D1^2 * D2^2 * 6 * K1^2 * F1^2 * K2^2 * F2^2 * C^2)",
        "(L * D1^3 * D2^3 * 4 * K1^3 * F1^3 * K2^3 * F2^3 * C^3)",
        "(L * D1^4 * D2^4 * 1 * K1^4 * F1^4 * K2^4 * F2^4 * C^4)"
    )

open_states <-
    list(
        "(L)",
        "(L * D1^1 * 4 * K1^1 * F1^1)",
        "(L * D1^2 * 6 * K1^2 * F1^2)",
        "(L * D1^3 * 4 * K1^3 * F1^3)",
        "(L * D1^4 * 1 * K1^4 * F1^4)",
        "(L * D2^1 * 4 * K2^1 * F2^1)",
        "(L * D2^2 * 6 * K2^2 * F2^2)",
        "(L * D2^3 * 4 * K2^3 * F2^3)",
        "(L * D2^4 * 1 * K2^4 * F2^4)",
        "(L * D1^1 * D2^1 * 4 * K1^1 * F1^1 * K2^1 * F2^1 * C^1)",
        "(L * D1^2 * D2^2 * 6 * K1^2 * F1^2 * K2^2 * F2^2 * C^2)",
        "(L * D1^3 * D2^3 * 4 * K1^3 * F1^3 * K2^3 * F2^3 * C^3)",
        "(L * D1^4 * D2^4 * 1 * K1^4 * F1^4 * K2^4 * F2^4 * C^4)"
    )

bound_states <-
    list(
        "(1 * 4 * K1^1 * F1^1)",
        "(2 * 6 * K1^2 * F1^2)",
        "(3 * 4 * K1^3 * F1^3)",
        "(4 * 1 * K1^4 * F1^4)",
        "(1 * L * D1^1 * 4 * K1^1 * F1^1)",
        "(2 * L * D1^2 * 6 * K1^2 * F1^2)",
        "(3 * L * D1^3 * 4 * K1^3 * F1^3)",
        "(4 * L * D1^4 * 1 * K1^4 * F1^4)",
        "(1 * 4 * K1^1 * F1^1 * K2^1 * F2^1 * C^1)",
        "(2 * 6 * K1^2 * F1^2 * K2^2 * F2^2 * C^2)",
        "(3 * 4 * K1^3 * F1^3 * K2^3 * F2^3 * C^3)",
        "(4 * 1 * K1^4 * F1^4 * K2^4 * F2^4 * C^4)",
        "(1 * L * D1^1 * D2^1 * 4 * K1^1 * F1^1 * K2^1 * F2^1 * C^1)",
        "(2 * L * D1^2 * D2^2 * 6 * K1^2 * F1^2 * K2^2 * F2^2 * C^2)",
        "(3 * L * D1^3 * D2^3 * 4 * K1^3 * F1^3 * K2^3 * F2^3 * C^3)",
        "(4 * L * D1^4 * D2^4 * 1 * K1^4 * F1^4 * K2^4 * F2^4 * C^4)"
    )

binding_eq <-
    paste(
        "(",
        paste(bound_states, collapse = " + "),
        ") / (4 * (",
        paste(all_states, collapse = " + "),
        "))",
        sep = ""
    )

gating_eq <-
    paste(
        "((",
        paste(open_states, collapse = " + "),
        ") / (",
        paste(all_states, collapse = " + "),
        "))",
        sep = ""
    )

args <- c("F1, F2, L, K1, D1, K2, D2, C")

paste('f <- function(', args, ') { return(' ,binding_eq , ')}', sep='') %>%
parse(text = .) %>%
eval() -> binding_function

paste('f <- function(', args, ') { return(' ,gating_eq , ')}', sep='') %>%
parse(text = .) %>%
eval() -> gating_function

simplest_model <- 
	expand.grid(
		"ATP" = 10^seq(-6, -2, length.out = 31),
		"PIP2" = 0,
		"L" = c(0.1, 1, 10),
		"KA_ATP" = c(10^3, 10^4, 10^5),
		"KA_PIP2" = 1,
		"D_ATP" = c(0.1, 0.5, 0.9),
		"D_PIP2" = 1,
		"C" = 1) %>%
	tibble() %>%
	mutate(
		atp_bound = binding_function(F1=ATP, F2=PIP2, L=L, K1=KA_ATP, D1=D_ATP, K2=KA_PIP2, D2=D_PIP2, C=C),
		pip2_bound = binding_function(F1=PIP2, F2=ATP, L=L, K1=KA_PIP2, D1=D_PIP2, K2=KA_ATP, D2=D_ATP, C=C),
		open = gating_function(F1=ATP, F2=PIP2, L=L, K1=KA_ATP, D1=D_ATP, K2=KA_PIP2, D2=D_PIP2, C=C),
		open_normalised = open / (L/(L+1))
		) %>%
	pivot_longer(atp_bound:open_normalised, names_to = "measure", values_to = "fraction")

ggplot(simplest_model %>% filter(L == 1, D_ATP == 0.5), aes(x = ATP, y = fraction, colour = factor(KA_ATP))) +
geom_line() +
facet_grid(cols = vars(measure), rows = vars(interaction(L, D_ATP)), labeller = label_both) +
scale_x_log10() +
theme_thesis() -> simple_a

ggplot(simplest_model %>% filter(KA_ATP == 10^4, D_ATP == 0.5), aes(x = ATP, y = fraction, colour = factor(L))) +
geom_line() +
facet_grid(cols = vars(measure), rows = vars(interaction(KA_ATP, D_ATP)), labeller = label_both) +
scale_x_log10() +
theme_thesis() -> simple_b

ggplot(simplest_model %>% filter(KA_ATP == 10^4, L == 1), aes(x = ATP, y = fraction, colour = factor(D_ATP))) +
geom_line() +
facet_grid(cols = vars(measure), rows = vars(interaction(L, KA_ATP)), labeller = label_both) +
scale_x_log10() +
theme_thesis() -> simple_c

simple_a/simple_b/simple_c


pip_model <- 
	expand.grid(
		"ATP" = 10^seq(-6, -2, length.out = 31),
		"PIP2" = 10^seq(-6, -2, length.out = 31),
		"L" = c(0.1, 1, 10),
		"KA_ATP" = 10^4,
		"KA_PIP2" = c(10^3, 10^4, 10^5),
		"D_ATP" = 0.5,
		"D_PIP2" = c(2, 5, 10),
		"C" = 1) %>%
	tibble() %>%
	mutate(
		atp_bound = binding_function(F1=ATP, F2=PIP2, L=L, K1=KA_ATP, D1=D_ATP, K2=KA_PIP2, D2=D_PIP2, C=C),
		pip2_bound = binding_function(F1=PIP2, F2=ATP, L=L, K1=KA_PIP2, D1=D_PIP2, K2=KA_ATP, D2=D_ATP, C=C),
		open = gating_function(F1=ATP, F2=PIP2, L=L, K1=KA_ATP, D1=D_ATP, K2=KA_PIP2, D2=D_PIP2, C=C),
		) %>%
	group_by(L, PIP2, KA_ATP, KA_PIP2, D_ATP, D_PIP2, C) %>%
	mutate(open_normalised = open/max(open)) %>%
	pivot_longer(atp_bound:open_normalised, names_to = "measure", values_to = "fraction")

ggplot(pip_model %>% filter(L == 0.1, D_PIP2 == 5, PIP2 == 1e-4), aes(x = ATP, y = fraction, colour = factor(KA_PIP2))) +
geom_line() +
facet_grid(cols = vars(measure), rows = vars(interaction(PIP2, L, D_ATP, D_PIP2)), labeller = label_both) +
scale_x_log10() +
theme_thesis() -> pip_a

ggplot(simplest_model %>% filter(KA_ATP == 10^4, D_ATP == 0.5), aes(x = ATP, y = fraction, colour = factor(L))) +
geom_line() +
facet_grid(cols = vars(measure), rows = vars(interaction(KA_ATP, D_ATP)), labeller = label_both) +
scale_x_log10() +
theme_thesis() -> simple_b

ggplot(simplest_model %>% filter(KA_ATP == 10^4, L == 1), aes(x = ATP, y = fraction, colour = factor(D_ATP))) +
geom_line() +
facet_grid(cols = vars(measure), rows = vars(interaction(L, KA_ATP)), labeller = label_both) +
scale_x_log10() +
theme_thesis() -> simple_c

simple_a/simple_b/simple_c
