library(tidyverse)
library(brms)

abfs_1 <- c(
    "/home/sam/data/patch_fluorometry_raw/180524_patching/W311-GFP+SUR/patch1/18524003.abf",
    "/home/sam/data/patch_fluorometry_raw/180524_patching/W311-GFP+SUR/patch4/18524006.abf",
    "/home/sam/data/patch_fluorometry_raw/180524_patching/W311-GFP+SUR/patch6/18524012.abf",
    "/home/sam/data/patch_fluorometry_raw/180524_patching/W311-GFP+SUR/patch8/18524018.abf",
    "/home/sam/data/patch_fluorometry_raw/180529_patching/W311-GFP+SUR/patch2/18529011.abf",
    "/home/sam/data/patch_fluorometry_raw/180529_patching/W311-GFP+SUR/patch4/18529019.abf",
    "/home/sam/data/patch_fluorometry_raw/180529_patching/W311-GFP+SUR/patch5/18529021.abf",
    "/home/sam/data/patch_fluorometry_raw/180529_patching/W311-GFP+SUR/patch7/18529026.abf",
    "/home/sam/data/patch_fluorometry_raw/201102/W311-GFP+SUR1/patch1/20n02004.abf",
    "/home/sam/data/patch_fluorometry_raw/201102/W311-GFP+SUR1/patch2/20n02005.abf",
    "/home/sam/data/patch_fluorometry_raw/201102/W311-GFP+SUR1/patch5/20n02014.abf",
    "/home/sam/data/patch_fluorometry_raw/201102/W311-GFP+SUR1/patch3/20n02005.abf"
)


abfs_2 <- c(
    "/home/sam/data/patch_fluorometry_raw/180724_PATCHING/W311-C166S-GFP/patch1/18724003.abf",
    "/home/sam/data/patch_fluorometry_raw/180731_PATCHING/W311-C166S-GFP/patch1/18731000.abf",
    "/home/sam/data/patch_fluorometry_raw/180731_PATCHING/W311-C166S-GFP/patch3/18731008.abf",
    "/home/sam/data/patch_fluorometry_raw/180731_PATCHING/W311-C166S-GFP/patch4/18731012.abf",
    "/home/sam/data/patch_fluorometry_raw/181212_PATCHING/W311ANAP-C166S-GFP+SUR/patch1/18d12006.abf",
    "/home/sam/data/patch_fluorometry_raw/181212_PATCHING/W311ANAP-C166S-GFP+SUR/patch2/18d13001.abf"
)

abfs_3 <- c(
    "/home/sam/data/patch_fluorometry_raw/201126/WT-GFP+SUR1/20n26005.abf",
    "/home/sam/data/patch_fluorometry_raw/201126/WT-GFP+SUR1/20n26006.abf",
    "/home/sam/data/patch_fluorometry_raw/201126/WT-GFP+SUR1/20n26008.abf",
    "/home/sam/data/patch_fluorometry_raw/201123/WT-GFP+SUR1/20n23001.abf",
    "/home/sam/data/patch_fluorometry_raw/201123/WT-GFP+SUR1/20n23002.abf",
    "/home/sam/data/patch_fluorometry_raw/201123/WT-GFP+SUR1/20n23003.abf",
    "/home/sam/data/patch_fluorometry_raw/201105/WT-GFP+SUR1/20n05028.abf",
    "/home/sam/data/patch_fluorometry_raw/201105/WT-GFP+SUR1/20n05030.abf",
    "/home/sam/data/patch_fluorometry_raw/201105/WT-GFP+SUR1/20n05031.abf",
    "/home/sam/data/patch_fluorometry_raw/201006/WT-GFP+SUR1/20o06007.abf",
    "/home/sam/data/patch_fluorometry_raw/201006/WT-GFP+SUR1/20o06008.abf",
    "/home/sam/data/patch_fluorometry_raw/201006/WT-GFP+SUR1/20o06010.abf"    
)

data_1 <-
    abfs_1 %>%
    map(~ as.data.frame(readABF::readABF(.x))) %>%
    bind_rows(.id = "n") %>%
    as_tibble() %>%
    rename(s = `Time [s]`, pA = `IN 0 [pA]`) %>%
    mutate(
        second_chunks = floor(s)
        ) %>%
    group_by(n, second_chunks) %>%
    mutate(
        pA_average = mean(pA),
        variance = (pA - pA_average)^2,
        ) %>%
    summarise(
        pA = mean(pA_average),
        variance = mean(variance)
        )  %>%
    mutate(construct = "W311*-GFP+SUR")

data_2 <-
    abfs_2 %>%
    map(~ as.data.frame(readABF::readABF(.x))) %>%
    bind_rows(.id = "n") %>%
    as_tibble() %>%
    rename(s = `Time [s]`, pA = `IN 0 [pA]`) %>%
    mutate(
        second_chunks = floor(s)
        ) %>%
    group_by(n, second_chunks) %>%
    mutate(
        pA_average = mean(pA),
        variance = (pA - pA_average)^2,
        ) %>%
    summarise(
        pA = mean(pA_average),
        variance = mean(variance)
        )  %>%
    mutate(construct = "W311*-C166S-GFP+SUR")

data_3 <-
    abfs_3 %>%
    map(~ as.data.frame(readABF::readABF(.x))) %>%
    bind_rows(.id = "n") %>%
    as_tibble() %>%
    rename(s = `Time [s]`, pA = `IN 0 [pA]`) %>%
    mutate(
        second_chunks = floor(s)
        ) %>%
    group_by(n, second_chunks) %>%
    mutate(
        pA_average = mean(pA),
        variance = (pA - pA_average)^2,
        ) %>%
    summarise(
        pA = mean(pA_average),
        variance = mean(variance)
        )  %>%
    mutate(construct = "WT-GFP+SUR")

bind_rows(data_1, data_2, data_3) -> summarised_data

ggplot(summarised_data) +
geom_point(aes(x = second_chunks, y = pA, colour = log(variance))) +
facet_grid(cols = vars(construct), rows = vars(n), scales = "free") +
scale_colour_distiller(palette = "PuOr")

ggplot(summarised_data) +
geom_point(aes(x = interaction(construct, n), y = log(variance), colour = second_chunks)) +
scale_colour_distiller(palette = "PuOr")

filtered_summarised_data <-
    summarised_data %>%
    mutate(rounded_pA = signif(pA, 2)) %>%
    group_by(rounded_pA) %>%
    filter(variance > 0, variance < 20000, second_chunks > 20, variance > quantile(variance, 0.025), variance < quantile(variance, 0.975))

ggplot(filtered_summarised_data) +
geom_point(aes(x = pA, y = variance)) +
facet_wrap(vars(construct, n), scales = "free")

ggplot(filtered_summarised_data) +
geom_point(aes(x = second_chunks, y = rounded_pA, colour = sqrt(variance))) +
facet_grid(rows = vars(construct, n), scales = "free") +
scale_colour_distiller(palette = "PuOr")

barium_data <-
    filtered_summarised_data %>%
    filter(pA < 0, variance < 100) %>%
    group_by(construct, n) %>%
    summarise(
        ba_pA = mean(-pA),
        ba_variance = mean(variance)
        )

filtered_summarised_data %>%
left_join(barium_data) %>%
replace_na(list(ba_pA = 0, ba_variance = 0)) %>%
mutate(corrected_pA = pA - ba_pA, corrected_variance = variance - ba_variance) %>%
filter(corrected_pA >= 0) %>%
group_by(construct, n) %>%
mutate(
    linear = case_when(corrected_pA <= mean(corrected_pA) & corrected_pA <= 2000 ~ TRUE, TRUE ~ FALSE),
    second_chunks = second_chunks - min(second_chunks)
    ) -> cleaned_data

ggplot(cleaned_data) +
geom_point(aes(x = corrected_pA, y = corrected_variance, colour = linear)) +
facet_wrap(vars(construct, n), scales = "free")

fits_n_freei <-
    cleaned_data %>%
    group_by(construct, n) %>%
    nest() %>%
    mutate(
        fit = map(data, ~ nls.multstart::nls_multstart(
            formula = corrected_variance ~ (i * corrected_pA) - ((corrected_pA ^ 2) / nchannels),
            data = data.frame(.),
            iter=500,
            start_lower = list(i = 2, nchannels = 100),
            start_upper = list(i = 6, nchannels = 100000),
            control = nls.control(warnOnly = TRUE)
            )),
        tidy = map(fit, broom::tidy),
        augment = map(fit, broom::augment, newdata = tibble(corrected_pA = seq(0, 10000, length.out = 51)))
        )

fits_n_fixedi <-
    cleaned_data %>%
    group_by(construct, n) %>%
    mutate(max_I = max(corrected_pA)) %>%
    nest() %>%
    mutate(
        fit = map(data, ~ nls.multstart::nls_multstart(
            formula = corrected_variance ~ (4 * corrected_pA) - ((corrected_pA ^ 2) / nchannels),
            data = data.frame(.),
            iter=500,
            start_lower = list(nchannels = 100),
            start_upper = list(nchannels = 10000),
            control = nls.control(warnOnly = TRUE)
            )),
        tidy = map(fit, broom::tidy),
        augment = map(fit, broom::augment, newdata = tibble(corrected_pA = seq(0, 10000, length.out = 51)))
        )

fits_popen <-
    cleaned_data %>%
    filter(linear == TRUE) %>%
    group_by(construct, n) %>%
    nest() %>%
    mutate(
        fit = map(data, ~ lm(
            formula = corrected_variance ~ 0 + corrected_pA,
            data = data.frame(.)
            )),
        tidy = map(fit, broom::tidy),
        augment = map(fit, broom::augment, newdata = tibble(corrected_pA = seq(0, 10000, length.out = 51)))
        )

fits_n_fixedi <-
    cleaned_data %>%
    group_by(construct, n) %>%
    mutate(max_I = max(corrected_pA)) %>%
    nest() %>%
    mutate(
        fit = map(data, ~ nls.multstart::nls_multstart(
            formula = corrected_variance ~ (4 * corrected_pA) - ((corrected_pA ^ 2) / nchannels),
            data = data.frame(.),
            iter=500,
            start_lower = list(nchannels = 100),
            start_upper = list(nchannels = 10000),
            control = nls.control(warnOnly = TRUE)
            )),
        tidy = map(fit, broom::tidy),
        augment = map(fit, broom::augment, newdata = tibble(corrected_pA = seq(0, 10000, length.out = 51)))
        )

bind_rows(
    fits_n_freei %>% unnest(augment) %>% mutate(model = "free_n_free_i"),
    fits_n_fixedi %>% unnest(augment) %>% mutate(model = "free_n_fixed_i"),
    fits_popen %>% unnest(augment) %>% mutate(model = "free_po")
    ) -> plot_fits

ggplot() +
geom_point(data = cleaned_data, aes(x = corrected_pA, y = corrected_variance, colour = linear)) +
geom_line(data = plot_fits, aes(x = corrected_pA, y = .fitted, colour = model)) +
facet_wrap(vars(construct, n), scales = "free") +
coord_cartesian(xlim = c(0, 7500), ylim = c(0, 12000))

bind_rows(
    fits_n_freei %>% unnest(tidy) %>% mutate(model = "free_n_free_i"),
    fits_n_fixedi %>% unnest(tidy) %>% mutate(model = "free_n_fixed_i"),
    fits_popen %>% unnest(tidy) %>% mutate(model = "free_po", term = "popen", estimate = 1 - (estimate / 4))
    ) %>%
select(construct, n, term, estimate, model) -> tidy_fits

cleaned_data %>%
group_by(construct, n) %>%
filter(corrected_pA == max(corrected_pA)) %>%
right_join(tidy_fits %>% pivot_wider(names_from = term, values_from = estimate)) %>%
mutate(
    i = case_when(is.na(i) ~ 4, TRUE ~ i),
    popen = case_when(is.na(popen) ~ corrected_pA / (nchannels * i), TRUE ~ popen)
    ) %>%
pivot_longer(c(i, nchannels, popen), names_to = "term", values_to = "estimate") %>%
group_by(construct, model, n) %>%
filter(is.na(estimate) | (estimate < 1e4 & estimate > 0)) -> tidy_fits_2

tidy_fits_2 %>%
group_by(construct, term, model) %>%
summarise(mu = mean(estimate), sigma = sd(estimate)) %>%
rowwise() %>%
mutate(dist = case_when(
    term == "popen" ~ dist_truncated(dist_normal(mu, sigma), 0, 1),
    TRUE ~ dist_normal(mu, sigma)
    )) -> tidy_fits_3

ggplot() +
geom_point(data = tidy_fits_2, aes(x = construct, y = estimate, colour = model)) +
facet_grid(rows = vars(term), cols = vars(model), scales = "free") +
stat_dist_slab(data = tidy_fits_3, aes(x = construct, dist = dist, colour = model), fill = NA)

brms_data <-
    cleaned_data %>%
    select(n, construct, second_chunks, corrected_pA, corrected_variance) %>%
    unite("construct_n", c(construct, n)) %>%
    mutate(construct_n = factor(construct_n))

noise_formula <-
    bf(
        corrected_variance ~ (4 * corrected_pA) - ((corrected_pA ^ 2) / exp(nchannels)),
        nchannels ~ 0 + Intercept + second_chunks + (1 + second_chunks||construct_n),
        family = gaussian(),
        nl = TRUE
        )

noise_priors <-
    c(
        prior(normal(7, 2), nlpar = nchannels, class = b, coef = Intercept),
        #prior for change in n over time
        prior(normal(-2, 1), nlpar = nchannels, class = b, coef = second_chunks)
    )

brm(
    formula = noise_formula,
    prior = noise_priors,
    data = brms_data %>% filter(construct_n %in% c("WT-GFP+SUR_1", "WT-GFP+SUR_2", "WT-GFP+SUR_3", "WT-GFP+SUR_4", "WT-GFP+SUR_5")),
    chains = 1,
    iter = 2e3,
    warmup = 1e3,
    thin  = 1,
    # seed = brms_seed,
    control = list(adapt_delta = 0.99, max_treedepth = 15),
    # cores = getOption("mc.cores", 4),
    # file = "data/other_fits/mixed_noise_analysis",
    sample_prior = "yes",
    save_all_pars = TRUE
    ) -> fit_both

brms_data %>%
filter(construct_n %in% c("WT-GFP+SUR_1")) %>%
mutate(construct_n = factor(construct_n)) %>%
expand(
    nesting(construct_n, second_chunks),
    corrected_pA = seq(min(corrected_pA), max(corrected_pA), length.out = 21)
    ) %>%
add_fitted_draws(fit_both, n = 5) -> fitted_draws_1

brms_data %>%
filter(construct_n %in% c("WT-GFP+SUR_1")) %>%
mutate(
    construct_n = factor(construct_n)
    ) -> plot_data

ggplot() +
geom_line(data = fitted_draws_1, aes(x = corrected_pA, y = .value, colour = second_chunks, group = interaction(second_chunks, .draw)), size = 0.1) +
geom_point(data = plot_data, aes(x = corrected_pA, y = corrected_variance, colour = second_chunks)) +
coord_cartesian(ylim = c(-1000, 10000)) +
scale_colour_distiller(palette = "PuOr")

abf_1 <-
    readABF::readABF("/home/sam/data/patch_fluorometry_raw/201123/WT-GFP+SUR1/20n23001.abf") %>%
    as.data.frame() %>%
    as_tibble() %>%
    rename(s = `Time [s]`, pA = `IN 0 [pA]`)

abf_1 %>%
mutate(
    time_point = case_when(
    (s > 65 & s < 70) ~ 1,
    (s > 144 & s < 149) ~ 2,
    (s > 215 & s < 220) ~ 3,
    (s > 255 & s < 260) ~ 4,
    (s > 283 & s < 288) ~ 5,
    (s > 440 & s < 444) ~ 6,
    TRUE ~ 0
        )
    ) %>%
filter(time_point !=0) %>%
group_by(time_point) %>%
mutate(
    pA_average = mean(pA),
    variance = (pA - pA_average)^2,
    ) %>%
summarise(
    pA = mean(pA_average),
    variance = mean(variance)
    ) -> wt_1

abf_2 <-
    readABF::readABF("/home/sam/data/patch_fluorometry_raw/201123/WT-GFP+SUR1/20n23002.abf") %>%
    as.data.frame() %>%
    as_tibble() %>%
    rename(s = `Time [s]`, pA = `IN 0 [pA]`)

abf_2 %>%
mutate(
    time_point = case_when(
    (s > 23 & s <28) ~ 1,
    (s > 32 & s < 37) ~ 2,
    (s > 111 & s < 116) ~ 3,
    (s > 166 & s < 171) ~ 4,
    (s > 191 & s < 196) ~ 5,
    (s > 247 & s < 252) ~ 6,
    (s > 276 & s < 281) ~ 7,
    (s > 409 & s < 414) ~ 8,
    TRUE ~ 0
        )
    ) %>%
filter(time_point !=0) %>%
group_by(time_point) %>%
mutate(
    pA_average = mean(pA),
    variance = (pA - pA_average)^2,
    ) %>%
summarise(
    pA = mean(pA_average),
    variance = mean(variance)
    ) -> wt_2

abf_3 <-
    readABF::readABF("/home/sam/data/patch_fluorometry_raw/201123/WT-GFP+SUR1/20n23003.abf") %>%
    as.data.frame() %>%
    as_tibble() %>%
    rename(s = `Time [s]`, pA = `IN 0 [pA]`)

abf_3 %>%
mutate(
    time_point = case_when(
    (s > 45 & s < 50) ~ 1,
    (s > 83 & s < 88) ~ 2,
    (s > 99 & s < 104) ~ 3,
    (s > 130 & s < 135) ~ 4,
    (s > 163 & s < 168) ~ 5,
    (s > 212 & s < 217) ~ 6,
    (s > 363 & s < 368) ~ 7,
    TRUE ~ 0
        )
    ) %>%
filter(time_point !=0) %>%
group_by(time_point) -> abf_3a

abf_3a %>%
mutate(
    pA_average = mean(pA),
    variance = (pA - pA_average)^2,
    ) %>%
summarise(
    pA = mean(pA_average),
    variance = mean(variance)
    ) -> wt_3

wt_1 %>%
filter(time_point == 6) %>%
rename(ba_pA = pA, ba_var = variance) %>%
mutate(ba_pA = -ba_pA) %>%
select(-time_point) %>%
bind_cols(wt_1 %>% filter(time_point !=6)) %>%
mutate(pA = pA - ba_pA, variance = variance - ba_var, n = 1) -> wt_1_cor

wt_2 %>%
filter(time_point == 8) %>%
rename(ba_pA = pA, ba_var = variance) %>%
mutate(ba_pA = -ba_pA) %>%
select(-time_point) %>%
bind_cols(wt_2 %>% filter(time_point !=8)) %>%
mutate(pA = pA - ba_pA, variance = variance - ba_var, n = 2) -> wt_2_cor

wt_3 %>%
filter(time_point == 7) %>%
rename(ba_pA = pA, ba_var = variance) %>%
mutate(ba_pA = -ba_pA) %>%
select(-time_point) %>%
bind_cols(wt_3 %>% filter(time_point !=7)) %>%
mutate(pA = pA - ba_pA, variance = variance - ba_var, n = 3) -> wt_3_cor

wt_cor <-
    bind_rows(wt_1_cor, wt_2_cor, wt_3_cor)

ggplot() +
geom_point(data = wt_cor, aes(x = pA, y = variance, colour = factor(n)))

free <-
    wt_cor %>%
    group_by(n) %>%
    nest() %>%
    mutate(
        fit = map(data, ~ nls.multstart::nls_multstart(
            formula = variance ~ (i * pA) - ((pA ^ 2) / nchannels),
            data = data.frame(.),
            iter=500,
            start_lower = list(i = 2, nchannels = 100),
            start_upper = list(i = 6, nchannels = 100000),
            control = nls.control(warnOnly = TRUE)
            ))) %>%
    mutate(augmented = map(fit, broom::augment, newdata = tibble(pA = seq(0, 7500, length.out = 51))))

fixed <-
    wt_cor %>%
    group_by(n) %>%
    nest() %>%
    mutate(
        fit = map(data, ~ nls.multstart::nls_multstart(
            formula = variance ~ (4 * pA) - ((pA ^ 2) / nchannels),
            data = data.frame(.),
            iter=500,
            start_lower = list(nchannels = 100),
            start_upper = list(nchannels = 100000),
            control = nls.control(warnOnly = TRUE)
            ))) %>%
    mutate(augmented = map(fit, broom::augment, newdata = tibble(pA = seq(0, 7500, length.out = 51))))

ggplot() +
geom_point(data = wt_cor, aes(x = pA, y = variance, colour = factor(n))) +
geom_line(
    data = free %>% unnest(augmented),
    aes(x = pA, y = .fitted, colour = factor(n)), linetype = 1) +
geom_line(
    data = fixed %>% unnest(augmented),
    aes(x = pA, y = .fitted, colour = factor(n)), linetype = 2) +
coord_cartesian(ylim = c(0, 8000)) +
facet_wrap(vars(n))


expand.grid(
    N = 1000,
    i = 4,
    popen = seq(0.1, 0.9, length.out = 6),
    time = seq(0, 1000, length.out = 1001)
    ) %>%
as_tibble() %>%
group_by(popen) %>%
mutate(
    nopen = rbinom(length(time), N, popen),
    pA = nopen*i
    ) -> test_grid

ggplot(test_grid, aes(x = time, y = pA, colour = popen)) +
geom_line() +
facet_wrap(vars(popen))

test_grid %>%
mutate(time_chunks = floor(time/10)) %>%
group_by(N, i, popen, time_chunks) %>%
summarise(nopen = mean(nopen), pA = mean(pA)) -> test_grid_filtered

ggplot(test_grid_filtered, aes(x = time_chunks, y = pA, colour = popen)) +
geom_line() +
facet_wrap(vars(popen))

test_grid %>%
group_by(popen) %>%
mutate(
    pA_average = mean(pA),
    variance = (pA - pA_average)^2,
    ) %>%
summarise(
    pA = mean(pA_average),
    variance = mean(variance)
    ) -> test_grid_summarised

nls(
    formula = variance ~ (i * pA) - ((pA ^ 2) / nchannels),
    data = test_grid_summarised,
    start = list(i=4, nchannels = 1000)
    ) -> fit_1

ggplot() +
geom_point(data = test_grid_summarised, aes(x = pA, y = variance)) +
geom_line(data = fit_1 %>% broom::augment(newdata = tibble(pA = seq(0, 4000, length.out = 51))), aes(x = pA, y = .fitted))

test_grid_filtered %>%
group_by(popen) %>%
mutate(
    pA_average = mean(pA),
    variance = (pA - pA_average)^2,
    ) %>%
summarise(
    pA = mean(pA_average),
    variance = mean(variance)
    ) -> test_grid_filtered_summarised

nls(
    formula = variance ~ (i * pA) - ((pA ^ 2) / nchannels),
    data = test_grid_filtered_summarised,
    start = list(i=4, nchannels = 1000)
    ) -> fit_2

ggplot() +
geom_point(data = test_grid_filtered_summarised, aes(x = pA, y = variance)) +
geom_line(data = fit_2 %>% broom::augment(newdata = tibble(pA = seq(0, 4500, length.out = 51))), aes(x = pA, y = .fitted))

expand.grid(
    N = seq(0, 1000, length.out = 5),
    i = 4,
    popen = seq(0.1, 0.9, length.out = 5),
    time = seq(0, 1000, length.out = 1001)
    ) %>%
as_tibble() %>%
group_by(N) %>%
mutate(
    nopen = rbinom(length(time), N, popen) + rbinom(length(time), 1000-N, popen/10),
    pA = nopen*i
    ) -> test_grid_2

ggplot(test_grid_2, aes(x = time, y = pA, colour = N)) +
geom_line() +
facet_wrap(vars(N, popen))

test_grid_2 %>%
group_by(N, popen) %>%
mutate(
    pA_average = mean(pA),
    variance = (pA - pA_average)^2,
    ) %>%
summarise(
    pA = mean(pA_average),
    variance = mean(variance)
    ) -> test_grid_2_summarised

nls(
    formula = variance ~ (i * pA) - ((pA ^ 2) / nchannels),
    data = test_grid_2_summarised,
    start = list(i=4, nchannels = 1000)
    ) -> fit_3

nls(
    formula = variance ~ (4 * pA) - ((pA ^ 2) / nchannels),
    data = test_grid_2_summarised,
    start = list(nchannels = 1000)
    ) -> fit_4


lm(
    formula = variance ~ 0 + pA,
    data = test_grid_2_summarised %>% filter(pA < 500),
    ) -> fit_5

ggplot() +
geom_point(data = test_grid_2_summarised, aes(x = pA, y = variance)) +
geom_line(data = fit_3 %>% broom::augment(newdata = tibble(pA = seq(0, 4000, length.out = 51))), aes(x = pA, y = .fitted)) +
geom_line(data = fit_4 %>% broom::augment(newdata = tibble(pA = seq(0, 4000, length.out = 51))), aes(x = pA, y = .fitted), linetype = 2) +
geom_line(data = fit_5 %>% broom::augment(newdata = tibble(pA = seq(0, 2000, length.out = 51))), aes(x = pA, y = .fitted), linetype = 3)
