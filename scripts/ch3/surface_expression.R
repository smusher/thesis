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
        log_counts = log(counts)
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
    summarise(control_mean = mean(log_counts), control_sd = sd(log_counts))

centered <-
    surface_expression %>%
    mutate(
        centered_response = (log_counts - control$control_mean[1]) / control$control_sd[1],
        kir_construct = fct_relevel(kir_construct, "WT", "WT-GFP", "F183*", "F183*-GFP", "W311*", "W311*-GFP")
        )

ggplot(centered, aes(x = kir_construct, y = centered_response, fill = ha_present)) +
geom_point(shape=21) +
facet_grid(rows=vars(sur1_present), cols=vars(anap_present), labeller = label_both) +
coord_flip()

brm(
    centered_response ~ 1 + ha_present + (1 + anap_present * sur1_present||kir_construct),
    family = gaussian(),
    data = centered,
    cores = getOption("mc.cores", 4),
    sample_prior = "yes"
    ) -> model_1

centered %>%
expand(kir_construct, anap_present, sur1_present, ha_present) %>%
add_fitted_draws(model_1, n = 10, allow_new_levels=FALSE, dpar = c("mu", "sigma")) -> plot_1

ggplot() +
geom_point(data = centered, aes(x = kir_construct, y = centered_response, fill = ha_present), shape=21) +
stat_dist_slab(data = plot_1, position = position_nudge(x=0.2), aes(x = kir_construct, dist = dist_normal(mu, sigma), colour = ha_present), fill = NA, size = 0.5) +
facet_grid(cols=vars(sur1_present), rows=vars(anap_present), labeller = label_both) +
coord_flip() +
scale_colour_brewer(palette = "Pastel1", aesthetics = c("fill", "colour"))

hypothesis(model_1, "ha_presentTRUE > Intercept")

hypothesis(model_1, "anap_presentTRUE > 0", scope = "coef", group = "kir_construct")

hypothesis(model_1, "sur1_presentTRUE > 0", scope = "coef", group = "kir_construct")

hypothesis(model_1, "anap_presentTRUE:sur1_presentTRUE > sur1_presentTRUE", scope = "coef", group = "kir_construct")

hypothesis(model_1, "anap_presentTRUE:sur1_presentTRUE > sur1_presentTRUE", scope = "coef", group = "kir_construct")$samples %>%
as_tibble() %>%
rename(
    WT = "H1",
    `WT-GFP` = "H2",
    `F183*` = "H3",
    `F183*-GFP` = "H4",
    `W311*` = "H5",
    `W311*-GFP` = "H6"
    ) %>%
pivot_longer(everything(), names_to = "kir_construct", values_to = "posterior_prob") -> hypothesis_plot

ggplot() +
stat_slab(
    data = hypothesis_plot,
    aes(
        x = kir_construct,
        y = posterior_prob,
        fill = stat(cut_cdf_qi(cdf, .width = c(.5, .8, .95), labels = scales::percent_format()))
        )
    ) +
coord_flip() +
scale_fill_brewer(palette = "Blues", direction = -1, na.translate = FALSE) +
labs(fill = "Interval")