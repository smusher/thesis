library(tidyverse)
library(brms)
library(tidybayes)
library(ggdist)
library(distributional)

surface_expression <-
    bind_rows(
        read_csv("data/anap_construct_surface_expression.csv"),
        read_csv("data/tmd0_surface_expression_one.csv"),
        read_csv("data/tmd0_surface_expression_two.csv")
        ) %>%
    filter(sur_construct %in% c(NA, "SUR")) %>%
    mutate(
        sur1_present = case_when(sur_construct == "SUR" ~ TRUE, is.na(sur_construct) ~ FALSE),
        ha_present = str_detect(kir_construct, "HA"),
        kir_construct = str_replace_all(kir_construct, "-HA", ""),
        log_counts = log10(counts)
        ) %>%
    select(kir_construct, anap_present, sur1_present, ha_present, log_counts) %>%
    drop_na()

ggplot(surface_expression, aes(x = kir_construct, y = log_counts, fill = ha_present)) +
geom_point(shape=21) +
facet_grid(rows=vars(sur1_present), cols=vars(anap_present), labeller = label_both) +
coord_flip()

control <-
    surface_expression %>%
    filter(ha_present == FALSE) %>%
    group_by(kir_construct) %>%
    summarise(control_mean = mean(log_counts))

centered <-
    surface_expression %>%
    left_join(control) %>%
    mutate(
        centered_response = log_counts - control_mean,
        kir_construct = fct_relevel(kir_construct, "WT", "WT-GFP", "F183*", "F183*-GFP", "W311*", "W311*-GFP")
        )

centered_prior_calc <-
    centered %>%
    summarise(mu =  mean(centered_response), sigma = sd(centered_response))

ggplot() +
geom_point(data = centered, aes(x = kir_construct, y = centered_response, fill = ha_present), shape=21) +
stat_dist_slab(data = centered_prior_calc, aes(dist = dist_normal(mu, sigma))) +
facet_grid(rows=vars(sur1_present), cols=vars(anap_present), labeller = label_both) +
coord_flip()

priors <- c(
    prior(normal(0.928, 0.753), class = Intercept),
    prior(cauchy(0, 5), class = sigma),
    prior(cauchy(0, 1), class = sd, group = kir_construct)
    )

brm(
    centered_response ~ 1 + (1 + anap_present * sur1_present * ha_present||kir_construct),
    family = gaussian(),
    data = centered,
    prior = priors,
    cores = getOption("mc.cores", 4),
    sample_prior = "yes",
    save_all_pars = TRUE,
    chains = 4,
    iter = 45000,
    warmup= 5000,
    thin  = 10,
    control = list(adapt_delta = 0.99, max_treedepth = 15),
    file = "/home/sam/thesis/data/other_fits/surface_expression_2.rds"
    ) -> model_1

centered %>%
expand(nesting(kir_construct, anap_present, sur1_present, ha_present)) %>%
add_fitted_draws(model_1) -> plot_1

ggplot() +
stat_slab(data = plot_1, aes(x = kir_construct, y = .value, colour = interaction(anap_present, ha_present)), fill = NA, size = 1) +
geom_point(data = centered, position = position_dodge(width=0.2), aes(x = kir_construct, y = centered_response, fill = interaction(anap_present, ha_present)), shape=21, size = 2) +
facet_grid(cols=vars(sur1_present), labeller = label_both) +
coord_flip() +
scale_colour_brewer(palette = "Pastel1", aesthetics = c("fill", "colour"))

nd_1 <-
    centered %>%
    expand(kir_construct, ha_present, anap_present = TRUE, sur1_present = TRUE)

fitted_draws(model_1, newdata = nd_1) %>%
ungroup() %>%
select(kir_construct, ha_present, .draw, .value) %>%
pivot_wider(names_from = ha_present, names_prefix = "ha_present_", values_from = .value) %>%
mutate(ha_contrast = ha_present_TRUE - ha_present_FALSE) -> ha_contrasts_1

ggplot() +
stat_slab(data = ha_contrasts_1, aes(x = kir_construct, y = 10^ha_contrast, fill = stat(cut_cdf_qi(cdf, .width = c(.5, .8, .95), labels = scales::percent_format())))) +
coord_flip() +
scale_fill_brewer(palette = "Blues", direction = -1, na.translate = FALSE) +
labs(fill = "Interval", x = "Kir6.2 construct", y = "Fold increase in luminescence") +
scale_y_log10()

nd_2 <-
    centered %>%
    expand(kir_construct, ha_present, anap_present = FALSE, sur1_present = TRUE)

fitted_draws(model_1, newdata = nd_2) %>%
ungroup() %>%
select(kir_construct, ha_present, .draw, .value) %>%
pivot_wider(names_from = ha_present, names_prefix = "ha_present_", values_from = .value) %>%
mutate(ha_contrast = ha_present_TRUE - ha_present_FALSE) -> ha_contrasts_2

ggplot() +
stat_slab(data = ha_contrasts_2, aes(x = kir_construct, y = 10^ha_contrast, fill = stat(cut_cdf_qi(cdf, .width = c(.5, .8, .95), labels = scales::percent_format())))) +
coord_flip() +
scale_fill_brewer(palette = "Blues", direction = -1, na.translate = FALSE) +
labs(fill = "Interval", x = "Kir6.2 construct", y = "Fold increase in luminescence") +
scale_y_log10()

nd_3 <-
    centered %>%
    expand(kir_construct, anap_present, ha_present = TRUE, sur1_present = TRUE)

fitted_draws(model_1, newdata = nd_3) %>%
ungroup() %>%
select(kir_construct, anap_present, .draw, .value) %>%
pivot_wider(names_from = anap_present, names_prefix = "anap_present_", values_from = .value) %>%
mutate(anap_contrast = anap_present_TRUE - anap_present_FALSE) -> anap_contrasts_1

ggplot() +
stat_slab(data = anap_contrasts_1, aes(x = kir_construct, y = 10^anap_contrast, fill = stat(cut_cdf_qi(cdf, .width = c(.5, .8, .95), labels = scales::percent_format())))) +
coord_flip() +
scale_fill_brewer(palette = "Blues", direction = -1, na.translate = FALSE) +
labs(fill = "Interval", x = "Kir6.2 construct", y = "Fold increase in luminescence") +
scale_y_log10()
