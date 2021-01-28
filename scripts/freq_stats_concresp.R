library(tidyverse)
library(broom)
library(cowplot)
library(ggbeeswarm)
library(patchwork)

currents_1 <-
    read_csv("data/electrophys_data.csv") %>%
    select(n, method, measure, construct, nucleotide, concentration, response) %>%
    mutate(experimenter = "Sam")

currents_2 <-
    read_csv("data/natascias_data.csv") %>%
    select(n, method, measure, construct, nucleotide, concentration, response) %>%
    mutate(experimenter = "Natascia")

currents <-
    bind_rows(currents_1, currents_2) %>%
    group_by(experimenter, n, method, measure, construct, nucleotide) %>%
    mutate(unique_experiment_id = cur_group_id())

fluorescence <-
    read_csv("data/unroofed_concresp_data.csv") %>%
    filter(dye < 480) %>%
    select(unique_experiment_id, method, measure, construct, concentration, response) %>%
    mutate(nucleotide = "TNP-ATP", experimenter = "Sam")

pcf_data <-
    read_csv("data/pcf_data.csv") %>%
    filter(dye < 480 | is.na(dye)) %>%
    select(unique_experiment_id, method, measure, construct, concentration, response) %>%
    mutate(nucleotide = "TNP-ATP", experimenter = "Sam")

concresp <-
    bind_rows(currents, fluorescence, pcf_data) %>%
    mutate(log_concentration = log10(concentration)) %>%
    filter(log_concentration > -Inf) %>%
    group_by(unique_experiment_id, experimenter, method, measure, construct, nucleotide) %>%
    mutate(n = n()) %>%
    filter(n > 3)

concresp_summary <-
	concresp %>%
	group_by(method, measure, construct, nucleotide, log_concentration) %>%
	summarise(
		error = qnorm(0.975)*sd(response)/sqrt(length(response)),
		response = mean(response)
		)

hill_equation <- function(log_concentration, ec50, hill, floor){
    floor + ((1 - floor) / (1 + 10^((ec50 - log_concentration) * -hill)))
}

fancy_scientific <- function(l) {
    l <- format(l, scientific = TRUE)
    l <- gsub("^(.*)e", "10^", l)
    parse(text=l)
}

interesting_constructs <- 
	c(
	"WT-GFP+SUR", "W311*-GFP+SUR",
	"E179A-GFP+SUR", "W311*,E179A-GFP+SUR",
	"E179K-GFP+SUR", "W311*,E179K-GFP+SUR",
    "K39A-GFP+SUR", "W311*,K39A-GFP+SUR",
	"K39E-GFP+SUR", "W311*,K39E-GFP+SUR",
	"K39R-GFP+SUR", "W311*,K39R-GFP+SUR"
    )

colour_scheme <-
    c("#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3", "#a6d854", "#ffd92f") %>%
    rep(each = 2)
names(colour_scheme) <- interesting_constructs

free_floor <- readRDS("data/individual_hill_fits_free_floor.rds")
fixed_floor <- readRDS("data/individual_hill_fits_fixed_floor.rds")

current <-
    concresp %>%
    filter(measure == "current", nucleotide == "TNP-ATP", construct %in% interesting_constructs) %>%
    ungroup() %>%
    select(construct, log_concentration, response) %>%
    mutate(construct = factor(construct))

current_global <- nls(
    response ~ floor[construct] + ((1 - floor[construct])  / (1 + 10^((ec50[construct] - log_concentration) * -hill))),
    data = data.frame(current),
    start = list(ec50 = rep(-5, 12),floor = rep(0, 12), hill = 1),
    control = nls.control(warnOnly = TRUE)
    )

fluorescence <-
    concresp %>%
    filter(method == "pcf", measure == "fluorescence", nucleotide == "TNP-ATP", construct %in% interesting_constructs) %>%
    ungroup() %>%
    select(construct, log_concentration, response) %>%
    mutate(construct = factor(construct))

fluorescence_global <- nls(
    response ~ (1 / (1 + 10^((ec50[construct] - log_concentration) * -hill))),
    data = data.frame(fluorescence),
    start = list(ec50 = rep(-5, 6), hill = 1),
    control = nls.control(warnOnly = TRUE)
    )

current_global_tidy <- current_global %>% tidy()
fluorescence_global_tidy <- fluorescence_global %>% tidy()

current_global_ci<-
	current_global %>%
	nlstools::confint2() %>%
	as_tibble() %>%
	mutate(term = current_global_tidy$term) %>%
	rename(ci_lower = `2.5 %`, ci_upper = `97.5 %`)

fluorescence_global_ci<-
	fluorescence_global %>%
	nlstools::confint2() %>%
	as_tibble() %>%
	mutate(term = fluorescence_global_tidy$term) %>%
	rename(ci_lower = `2.5 %`, ci_upper = `97.5 %`)

current_global_params <-
	current_global_tidy %>%
	left_join(current_global_ci, by="term") %>%
	mutate(
		construct = as.integer(str_extract(term, "(?<=ec50).*|(?<=floor).*"))
		) %>%
	mutate(
		construct = factor(construct, labels = levels(current$construct)),
		term = str_replace_all(term, "(?<=ec50).*|(?<=floor).*", "")
		)

current_global_fits <-
	current_global_params %>%
	select(term, estimate, ci_lower, ci_upper, construct) %>%
	pivot_wider(names_from = term, values_from = c(estimate, ci_lower, ci_upper)) %>%
	fill(c(estimate_hill, ci_lower_hill, ci_upper_hill), .direction = "up") %>%
	drop_na() %>%
	full_join(
		expand(current, construct, log_concentration = seq(-7, -2, length.out = 51))
		) %>%
	mutate(
		estimate = hill_equation(log_concentration, estimate_ec50, estimate_hill, estimate_floor),
		lower = hill_equation(log_concentration, ci_lower_ec50, ci_lower_hill, ci_lower_floor),
		upper = hill_equation(log_concentration, ci_upper_ec50, ci_upper_hill, ci_upper_floor)
		)

fluorescence_global_params <-
	fluorescence_global_tidy %>%
	left_join(fluorescence_global_ci, by="term") %>%
	mutate(
		construct = as.integer(str_extract(term, "(?<=ec50).*"))
		) %>%
	mutate(
		construct = factor(construct, labels = levels(fluorescence$construct)),
		term = str_replace_all(term, "(?<=ec50).*", "")
		)

fluorescence_global_fits <-
	fluorescence_global_params %>%
	select(term, estimate, ci_lower, ci_upper, construct) %>%
	pivot_wider(names_from = term, values_from = c(estimate, ci_lower, ci_upper)) %>%
	fill(c(estimate_hill, ci_lower_hill, ci_upper_hill), .direction = "up") %>%
	drop_na() %>%
	full_join(
		expand(fluorescence, construct, log_concentration = seq(-7, -2, length.out = 51))
		) %>%
	mutate(
		estimate = hill_equation(log_concentration, estimate_ec50, estimate_hill, 0),
		lower = hill_equation(log_concentration, ci_lower_ec50, ci_lower_hill, 0),
		upper = hill_equation(log_concentration, ci_upper_ec50, ci_upper_hill, 0)
		)

keepers = vector("logical", nrow(free_floor))
for (i in seq_along(free_floor$fit)){
    keepers[[i]] <- 
        ifelse(
            is.null(free_floor$fit[[i]]$convInfo$isConv),
            FALSE,
            free_floor$fit[[i]]$convInfo$isConv
        )
}

free_fits <-
	free_floor %>%
	add_column(converged = keepers) %>%
    filter(converged == TRUE) %>%
    mutate(fit = map(fit, tidy)) %>%
    unnest(fit) %>%
	filter(construct %in% interesting_constructs, estimate < 2)

fixed_fits <-
	fixed_floor %>%
    mutate(fit = map(fit, tidy)) %>%
    unnest(fit) %>%
	filter(construct %in% interesting_constructs)

free_fits_plot <-
	free_floor %>%
	add_column(converged = keepers) %>%
    filter(converged == TRUE) %>%
    mutate(fit = map(fit, augment, newdata = tibble(log_concentration = seq(-7, -2, length.out = 51)))) %>%
    unnest(fit) %>%
	right_join(free_fits) %>%
	group_by(method, measure, construct, nucleotide, log_concentration) %>%
	summarise(
		error = qnorm(0.975)*sd(.fitted)/sqrt(length(.fitted)),
		response = mean(.fitted)
		)

fixed_fits_plot <-
	fixed_floor %>%
	add_column(converged = keepers) %>%
    filter(converged == TRUE) %>%
    mutate(fit = map(fit, augment, newdata = tibble(log_concentration = seq(-7, -2, length.out = 51)))) %>%
    unnest(fit) %>%
	right_join(free_fits) %>%
	group_by(method, measure, construct, nucleotide, log_concentration) %>%
	summarise(
		error = qnorm(0.975)*sd(.fitted)/sqrt(length(.fitted)),
		response = mean(.fitted)
		)

ggplot() +
geom_ribbon(
	data = free_fits_plot %>% filter(construct == "W311*-GFP+SUR" | construct == "W311*,K39R-GFP+SUR", measure == "current", nucleotide == "TNP-ATP"),
	aes(x = log_concentration, ymin = response - error, ymax = response + error, colour = construct),
	size = 1
	) +
geom_line(
	data = free_fits_plot %>% filter(construct == "W311*-GFP+SUR" | construct == "W311*,K39R-GFP+SUR", measure == "current", nucleotide == "TNP-ATP"),
	aes(x = log_concentration, y = response, colour = construct),
	size = 1
	) +
geom_pointrange(
	data = concresp_summary %>% filter(construct == "W311*-GFP+SUR" | construct == "W311*,K39R-GFP+SUR", measure == "current", nucleotide == "TNP-ATP"),
	aes(x = log_concentration, y = response, ymin = response-error, ymax = response+error, colour = construct)
	) +
theme_cowplot(font_family = "Liberation Sans", font_size = 11) +
scale_colour_manual(values = colour_scheme, aesthetics = c('colour', 'fill')) +
labs(x = "[TNP-ATP] (M)", y = expression(I / I[max]))

ggplot() +
geom_ribbon(
	data = fixed_fits_plot %>% filter(construct == "W311*-GFP+SUR" | construct == "W311*,K39R-GFP+SUR", measure == "fluorescence", method == "pcf"),
	aes(x = log_concentration, ymin = response - error, ymax = response + error, colour = construct),
	size = 1
	) +
geom_line(
	data = fixed_fits_plot %>% filter(construct == "W311*-GFP+SUR" | construct == "W311*,K39R-GFP+SUR", measure == "fluorescence", method == "pcf"),
	aes(x = log_concentration, y = response, colour = construct),
	size = 1
	) +
geom_pointrange(
	data = concresp_summary %>% filter(construct == "W311*-GFP+SUR" | construct == "W311*,K39R-GFP+SUR", measure == "fluorescence", method == "pcf"),
	aes(x = log_concentration, y = response, ymin = response-error, ymax = response+error, colour = construct)
	) +
theme_cowplot(font_family = "Liberation Sans", font_size = 11) +
scale_colour_manual(values = colour_scheme, aesthetics = c('colour', 'fill')) +
labs(x = "[TNP-ATP] (M)", y = expression(F / F[max]))

current_fits <-
	free_fits %>%
	filter(
		term == "ec50",
    	nucleotide != "MG-ATP",
    	measure == "current"
    ) %>%
    mutate(
    	Anap = case_when(stringr::str_detect(construct, "W311*") == TRUE ~ TRUE, TRUE ~ FALSE)
    )

current_fits_summary <-
	current_fits %>%
	group_by(Anap, construct, nucleotide) %>%
	summarise(
		error = qnorm(0.975)*sd(estimate)/sqrt(length(estimate)),
		estimate = mean(estimate)
		)	

fluorescence_fits <-
	fixed_fits %>%
	filter(
		measure == "fluorescence",
		term == "ec50",
    	method == "pcf"
    )

fluorescence_fits_summary <-
	fluorescence_fits %>%
	group_by(construct) %>%
	summarise(
		error = qnorm(0.975)*sd(estimate)/sqrt(length(estimate)),
		estimate = mean(estimate)
		)	

ggplot() +
geom_quasirandom(data = current_fits %>% filter(Anap == FALSE), aes(x = construct, y = 10^estimate, fill = construct), shape = 21, stroke = 1, size = 2.5, alpha = 0.5) +
geom_pointrange(data = current_fits_summary %>% filter(Anap == FALSE), aes(x = construct, ymin = 10^(estimate-error), y = 10^estimate, ymax = 10^(estimate+error), colour = construct)) +
theme_cowplot(font_family = "IBM Plex Sans Condensed", font_size = 15) +
theme(strip.background  = element_rect(fill = "white"), legend.position = "none") +
scale_colour_manual(values = colour_scheme, aesthetics = c('colour', 'fill')) +
scale_y_log10(labels = fancy_scientific, breaks = c(1e-6, 1e-5, 1e-4, 1e-3), limits = c(3e-7, 1e-3)) +
coord_flip() +
facet_grid(rows = vars(nucleotide)) +
labs(x = "Construct", y = "IC50 (M)") -> ic50_1

ggplot() +
geom_quasirandom(data = current_fits %>% filter(Anap == TRUE), aes(x = construct, y = 10^estimate, fill = construct), shape = 21, stroke = 1, size = 2.5, alpha = 0.5) +
geom_pointrange(data = current_fits_summary %>% filter(Anap == TRUE), aes(x = construct, ymin = 10^(estimate-error), y = 10^estimate, ymax = 10^(estimate+error), colour = construct)) +
theme_cowplot(font_family = "IBM Plex Sans Condensed", font_size = 15) +
theme(strip.background  = element_rect(fill = "white"), legend.position = "none") +
scale_colour_manual(values = colour_scheme, aesthetics = c('colour', 'fill')) +
scale_y_log10(labels = fancy_scientific, breaks = c(1e-6, 1e-5, 1e-4, 1e-3), limits = c(3e-7, 1e-3)) +
coord_flip() +
facet_grid(rows = vars(nucleotide)) +
labs(x = "Construct", y = "IC50 (M)") -> ic50_2

ggplot() +
geom_quasirandom(data = fluorescence_fits, aes(x = construct, y = 10^estimate, fill = construct), shape = 21, stroke = 1, size = 2.5, alpha = 0.5) +
geom_pointrange(data = fluorescence_fits_summary, aes(x = construct, ymin = 10^(estimate-error), y = 10^estimate, ymax = 10^(estimate+error), colour = construct)) +
theme_cowplot(font_family = "IBM Plex Sans Condensed", font_size = 15) +
theme(strip.background  = element_rect(fill = "white"), legend.position = "none") +
scale_colour_manual(values = colour_scheme, aesthetics = c('colour', 'fill')) +
scale_y_log10(labels = fancy_scientific, breaks = c(1e-6, 1e-5, 1e-4, 1e-3), limits = c(3e-7, 1e-3)) +
coord_flip() +
labs(x = "Construct", y = "EC50 (M)") -> ec50_1

ic50_1 + ic50_2 + ec50_1 -> ic50_plot

control <-
	current_fits %>%
	ungroup() %>%
	filter(construct %in% c("WT-GFP+SUR", "W311*-GFP+SUR")) %>%
	select(Anap, nucleotide, estimate) %>%
	group_by(Anap, nucleotide) %>%
	summarise(control_estimate = list(estimate))

test <-
	current_fits %>%
	ungroup() %>%
	filter(construct != "WT-GFP+SUR", construct != "W311*-GFP+SUR", unique_experiment_id != 12) %>%
	select(Anap, construct, nucleotide, estimate) %>%
	group_by(Anap, construct, nucleotide) %>%
	summarise(estimate = list(estimate)) %>%
	full_join(control)

t_tests <-
	test %>%
	group_by(Anap, construct, nucleotide) %>%
	mutate(
		p_value = tidy(t.test(unlist(estimate), unlist(control_estimate), data = .))$p.value
		) %>%
	group_by(Anap, nucleotide) %>%
	mutate(
		p_value_adjusted = p.adjust(p_value, method = "bonferroni", n = n())
		)

ic50_1 +
geom_label(
	data = t_tests %>% filter(Anap == FALSE),
	aes(x = construct, y = 1e-3,
		label = format(signif(p_value_adjusted, digits = 3), scientific = FALSE))
	) -> ic50_1a

ic50_2 +
geom_label(
	data = t_tests %>% filter(Anap == TRUE),
	aes(x = construct, y = 1e-6,
		label = format(signif(p_value_adjusted, digits = 3), scientific = FALSE))
	) -> ic50_2a

ic50_1a + ic50_2a -> ic50_plot_a

fluorescence_fits <-
	fluorescence_fits %>%
	mutate(construct = factor(construct, levels = interesting_constructs))

current_fits <-
	current_fits %>%
	mutate(construct = factor(construct, levels = interesting_constructs))

anova_anap_atp <- aov(
		estimate ~ construct,
		data = current_fits %>% filter(nucleotide == "ATP", Anap == TRUE, unique_experiment_id != 12)
	) %>%
multcomp::glht(linfct = multcomp::mcp(construct = "Dunnet")) %>%
tidy(conf.int = TRUE)

anova_atp <- aov(
		estimate ~ construct,
		data = current_fits %>% filter(nucleotide == "ATP", Anap == FALSE)
	) %>%
multcomp::glht(linfct = multcomp::mcp(construct = "Dunnet")) %>%
tidy(conf.int = TRUE) %>%
bind_rows(anova_anap_atp) %>%
mutate(nucleotide = "ATP", measure = "current")

anova_anap_tnpatp <- aov(
		estimate ~ construct,
		data = current_fits %>% filter(nucleotide == "TNP-ATP", Anap == TRUE)
	) %>%
multcomp::glht(linfct = multcomp::mcp(construct = "Dunnet")) %>%
tidy(conf.int = TRUE)

anova_tnpatp <- aov(
		estimate ~ construct,
		data = current_fits %>% filter(nucleotide == "TNP-ATP", Anap == FALSE)
	) %>%
multcomp::glht(linfct = multcomp::mcp(construct = "Dunnet")) %>%
tidy(conf.int = TRUE) %>%
bind_rows(anova_anap_tnpatp) %>%
mutate(nucleotide = "TNP-ATP", measure = "current")

anova_fluor <- aov(
		estimate ~ construct,
		data = fluorescence_fits %>% filter(method == "pcf")
	) %>%
multcomp::glht(linfct = multcomp::mcp(construct = "Dunnet")) %>%
tidy(conf.int = TRUE) %>%
mutate(nucleotide = "TNP-ATP", measure = "fluorescence")

anova_table <-
	bind_rows(anova_atp, anova_tnpatp, anova_fluor)

ggplot(data = anova_atp, aes(x = contrast, ymin = conf.low, y = estimate, ymax = conf.high), shape = 21) +
geom_pointrange() +
geom_label(aes(label = format(signif(adj.p.value, digits = 3), scientific = FALSE)), position = position_nudge(x = 0.25)) +
coord_flip() +
theme_cowplot(font_family = "IBM Plex Sans Condensed", font_size = 15) +
theme(strip.background  = element_rect(fill = "white"), legend.position = "none") +
labs(x = "Construct", y = "Confidence bands (Order of magnitude difference from control)") +
geom_hline(yintercept = 0, linetype = 2, size = 1) +
ggtitle("ATP") +
scale_y_continuous(limits = c(-0.5, 1)) -> aov_1

ggplot(data = anova_tnpatp, aes(x = contrast, ymin = conf.low, y = estimate, ymax = conf.high), shape = 21) +
geom_pointrange() +
geom_label(aes(label = format(signif(adj.p.value, digits = 3), scientific = FALSE)), position = position_nudge(x = 0.25)) +
coord_flip() +
theme_cowplot(font_family = "IBM Plex Sans Condensed", font_size = 15) +
theme(strip.background  = element_rect(fill = "white"), legend.position = "none") +
labs(x = "Construct", y = "Tukey confidence bands (Order of magnitude difference from control)") +
geom_hline(yintercept = 0, linetype = 2, size = 1) +
ggtitle("TNP-ATP") +
scale_y_continuous(limits = c(-0.5, 1)) -> aov_2

aov_1 + aov_2 -> aov_plot