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
    filter(chaperone != "TEA" | is.na(chaperone)) %>%
    mutate(
        sur_construct = case_when(is.na(sur_construct) ~ "None", TRUE ~ sur_construct),
        ha_present = str_detect(kir_construct, "HA"),
        kir_construct = str_replace_all(kir_construct, "-HA", ""),
        log_counts = log10(counts)
        ) %>%
    select(kir_construct, anap_present, sur_construct, ha_present, log_counts) %>%
    drop_na()

ggplot(surface_expression, aes(x = kir_construct, y = log_counts, fill = ha_present)) +
geom_point(shape=21) +
facet_grid(rows=vars(sur_construct), cols=vars(anap_present), labeller = label_both) +
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

write_csv(centered, "data/combined_surface_expression.csv")

centered_prior_calc <-
    centered %>%
    summarise(mu =  mean(centered_response), sigma = sd(centered_response))

ggplot() +
geom_point(data = centered, aes(x = kir_construct, y = centered_response, fill = ha_present), shape=21) +
stat_dist_slab(data = centered_prior_calc, aes(dist = dist_normal(mu, sigma))) +
facet_grid(rows=vars(sur_construct), cols=vars(anap_present), labeller = label_both) +
coord_flip()

priors <- c(
    prior(normal(1, 0.7), class = b),
    prior(cauchy(0, 1), dpar = sigma, class = Intercept),
    prior(cauchy(0, 1), dpar = sigma, class = sd)
    )

form_1 <-
    bf(
        centered_response ~ 0 + anap_present * sur_construct * ha_present * kir_construct,
        sigma ~ (1|ha_present)
        )

brm(
    formula = form_1,
    family = gaussian(),
    data = centered,
    prior = priors,
    cores = getOption("mc.cores", 4),
    sample_prior = "yes",
    save_all_pars = TRUE,
    chains = 4,
    iter = 40000,
    warmup= 20000,
    thin  = 10,
    control = list(adapt_delta = 0.99, max_treedepth = 15),
    file = "/home/sam/thesis/data/other_fits/surface_expression_3.rds"
    ) -> model_1

centered %>%
expand(nesting(kir_construct, anap_present, sur_construct, ha_present)) %>%
add_predicted_draws(model_1) -> plot_1

ggplot() +
stat_slab(data = plot_1 %>% filter(sur_construct != "TMD0_1363"), aes(x = kir_construct, y = .prediction, colour = interaction(anap_present, ha_present, sur_construct)), fill = NA, size = 1) +
geom_point(data = centered %>% filter(sur_construct != "TMD0_1363"), position = position_dodge(width=0.2), aes(x = kir_construct, y = centered_response, fill = interaction(anap_present, ha_present, sur_construct)), shape=21, size = 2) +
coord_flip() +
scale_colour_brewer(palette = "Pastel1", aesthetics = c("fill", "colour"))

nd_1 <-
    centered %>%
    expand(kir_construct, ha_present, anap_present = TRUE, sur_construct = "SUR")

predicted_draws(model_1, newdata = nd_1) %>%
ungroup() %>%
select(kir_construct, ha_present, .draw, .prediction) %>%
pivot_wider(names_from = ha_present, names_prefix = "ha_present_", values_from = .prediction) %>%
mutate(ha_contrast = ha_present_TRUE - ha_present_FALSE) -> ha_contrasts_1

ggplot() +
stat_slab(data = ha_contrasts_1, aes(x = kir_construct, y = 10^ha_contrast, fill = stat(cut_cdf_qi(cdf, .width = c(.5, .8, .95), labels = scales::percent_format())))) +
coord_flip() +
scale_fill_brewer(palette = "Blues", direction = -1, na.translate = FALSE) +
labs(fill = "Interval", x = "Kir6.2 construct", y = "Fold increase in luminescence") +
scale_y_log10()

nd_2 <-
    centered %>%
    expand(kir_construct, ha_present, anap_present = FALSE, sur_construct = "SUR")

predicted_draws(model_1, newdata = nd_2) %>%
ungroup() %>%
select(kir_construct, ha_present, .draw, .prediction) %>%
pivot_wider(names_from = ha_present, names_prefix = "ha_present_", values_from = .prediction) %>%
mutate(ha_contrast = ha_present_TRUE - ha_present_FALSE) -> ha_contrasts_2

ggplot() +
stat_slab(data = ha_contrasts_2, aes(x = kir_construct, y = 10^ha_contrast, fill = stat(cut_cdf_qi(cdf, .width = c(.5, .8, .95), labels = scales::percent_format())))) +
coord_flip() +
scale_fill_brewer(palette = "Blues", direction = -1, na.translate = FALSE) +
labs(fill = "Interval", x = "Kir6.2 construct", y = "Fold increase in luminescence") +
scale_y_log10()

nd_3 <-
    centered %>%
    expand(kir_construct, anap_present, ha_present = TRUE, sur_construct = "SUR")

predicted_draws(model_1, newdata = nd_3) %>%
ungroup() %>%
select(kir_construct, anap_present, .draw, .prediction) %>%
pivot_wider(names_from = anap_present, names_prefix = "anap_present_", values_from = .prediction) %>%
mutate(anap_contrast = anap_present_TRUE - anap_present_FALSE) -> anap_contrasts_1

ggplot() +
stat_slab(data = anap_contrasts_1, aes(x = kir_construct, y = 10^anap_contrast, fill = stat(cut_cdf_qi(cdf, .width = c(.5, .8, .95), labels = scales::percent_format())))) +
coord_flip() +
scale_fill_brewer(palette = "Blues", direction = -1, na.translate = FALSE) +
labs(fill = "Interval", x = "Kir6.2 construct", y = "Fold increase in luminescence") +
scale_y_log10()

nd_4 <-
    centered %>%
    expand(nesting(sur_construct, anap_present = FALSE, kir_construct = "WT-GFP"), ha_present = TRUE)

predicted_draws(model_1, newdata = nd_4) %>%
ungroup() %>%
select(kir_construct, sur_construct, .draw, .prediction) %>%
pivot_wider(names_from = sur_construct, values_from = .prediction) %>%
mutate(across(SUR:TMD0_232, ~.x - None, .names = "contrast_{.col}")) %>%
select(kir_construct, .draw, starts_with("contrast")) %>%
pivot_longer(starts_with("contrast"), names_to = "sur_construct", values_to = "contrast") -> sur_construct_contrasts_1

nd_5 <-
    centered %>%
    expand(nesting(sur_construct, anap_present = TRUE, kir_construct = "W311*-GFP"), ha_present = TRUE)

predicted_draws(model_1, newdata = nd_5) %>%
ungroup() %>%
select(kir_construct, sur_construct, .draw, .prediction) %>%
pivot_wider(names_from = sur_construct, values_from = .prediction) %>%
mutate(across(SUR:TMD0_232, ~.x - None, .names = "contrast_{.col}")) %>%
select(kir_construct, .draw, starts_with("contrast")) %>%
pivot_longer(starts_with("contrast"), names_to = "sur_construct", values_to = "contrast") -> sur_construct_contrasts_2

sur_construct_contrasts <-
    bind_rows(sur_construct_contrasts_1, sur_construct_contrasts_2)

ggplot() +
stat_slab(data =  sur_construct_contrasts %>% filter(sur_construct != "contrast_TMD0_1353"), aes(x = sur_construct, y = 10^contrast, fill = stat(cut_cdf_qi(cdf, .width = c(.5, .8, .95), labels = scales::percent_format())))) +
coord_flip() +
scale_fill_brewer(palette = "Blues", direction = -1, na.translate = FALSE) +
labs(fill = "Interval", x = "SUR construct", y = "Fold increase in luminescence") +
scale_y_log10() +
facet_wrap(vars(kir_construct))
