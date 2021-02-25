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
    "/home/sam/data/patch_fluorometry_raw/201105/WT-GFP+SUR1/20n05031.abf"    
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

data_3 <-
    abfs_3 %>%
    map(~ as.data.frame(readABF::readABF(.x))) %>%
    bind_rows(.id = "n") %>%
    as_tibble() %>%
    rename(s = `Time [s]`, pA = `IN 0 [pA]`)

data_3 %>%
group_by(n) %>%
nest() %>%
mutate(fit  = map(data, ~ loess(pA ~ s, data = data.frame(.)))) %>%
uunest(fit, augment)

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
    filter(variance < 1e4, variance > 0, second_chunks > 20)

ggplot(filtered_summarised_data) +
geom_point(aes(x = second_chunks, y = pA, colour = log(variance))) +
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
group_by(construct, n, pA) %>%
filter(corrected_variance < 2*corrected_pA) -> cleaned_data

ggplot(cleaned_data) +
geom_point(aes(x = corrected_pA, y = corrected_variance)) +
facet_wrap(vars(construct, n))

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
    nest() %>%
    mutate(
        fit = map(data, ~ nls.multstart::nls_multstart(
            formula = corrected_variance ~ (4.32 * corrected_pA) - ((corrected_pA ^ 2) / nchannels),
            data = data.frame(.),
            iter=500,
            start_lower = list(nchannels = 100),
            start_upper = list(nchannels = 100000),
            control = nls.control(warnOnly = TRUE)
            )),
        tidy = map(fit, broom::tidy),
        augment = map(fit, broom::augment, newdata = tibble(corrected_pA = seq(0, 10000, length.out = 51)))
        )

fits_popen_freei <-
    cleaned_data %>%
    group_by(construct, n) %>%
    nest() %>%
    mutate(
        fit = map(data, ~ nls.multstart::nls_multstart(
            formula = corrected_variance ~ i * corrected_pA * (1 - popen),
            data = data.frame(.),
            iter=500,
            start_lower = list(i = 2, popen = 0.1),
            start_upper = list(i = 6, popen = 0.9),
            control = nls.control(warnOnly = TRUE)
            )),
        tidy = map(fit, broom::tidy),
        augment = map(fit, broom::augment, newdata = tibble(corrected_pA = seq(0, 10000, length.out = 51)))
        )

fits_popen_fixedi <-
    cleaned_data %>%
    group_by(construct, n) %>%
    nest() %>%
    mutate(
        fit = map(data, ~ nls.multstart::nls_multstart(
            formula = corrected_variance ~ 4.32 * corrected_pA * (1 - popen),
            data = data.frame(.),
            iter=500,
            start_lower = list(popen = 0.1),
            start_upper = list(popen = 0.9),
            control = nls.control(warnOnly = TRUE)
            )),
        tidy = map(fit, broom::tidy),
        augment = map(fit, broom::augment, newdata = tibble(corrected_pA = seq(0, 10000, length.out = 51)))
        )

bind_rows(
    fits_n_freei %>% unnest(augment) %>% mutate(model = "free_n_free_i"),
    fits_n_fixedi %>% unnest(augment) %>% mutate(model = "free_n_fixed_i"),
    fits_popen_freei %>% unnest(augment) %>% mutate(model = "free_po_free_i"),
    fits_popen_fixedi %>% unnest(augment) %>% mutate(model = "free_po_fixed_i")
    ) -> plot_fits

ggplot() +
geom_point(data = cleaned_data, aes(x = corrected_pA, y = corrected_variance)) +
geom_line(data = plot_fits, aes(x = corrected_pA, y = .fitted, colour = model)) +
facet_wrap(vars(construct, n)) +
coord_cartesian(xlim = c(0, 7500), ylim = c(0, 12000))

bind_rows(
    fits_n_freei %>% unnest(tidy) %>% mutate(model = "free_n_free_i"),
    fits_n_fixedi %>% unnest(tidy) %>% mutate(model = "free_n_fixed_i"),
    fits_popen_freei %>% unnest(tidy) %>% mutate(model = "free_po_free_i"),
    fits_popen_fixedi %>% unnest(tidy) %>% mutate(model = "free_po_fixed_i")
    ) -> tidy_fits

cleaned_data %>%
group_by(construct, n) %>%
filter(corrected_pA == max(corrected_pA)) %>%
right_join(tidy_fits %>% pivot_wider(names_from = term, values_from = estimate)) %>%
mutate(
    i = case_when(is.na(i) ~ 4.32, TRUE ~ i),
    popen = case_when(is.na(popen) ~ corrected_pA / (nchannels * i), TRUE ~ popen)
    ) %>%
pivot_longer(c(i, nchannels, popen), names_to = "term", values_to = "estimate") -> tidy_fits_2

ggplot() +
geom_point(data = tidy_fits_2, aes(x = construct, y = estimate, colour = model)) +
facet_grid(rows = vars(term), scales = "free")

noise_formula <-
    bf(
        corrected_variance ~ (i * corrected_pA) - ((corrected_pA ^ 2) / nchannels),
        i ~ 1 + (1||n),
        nchannels ~ 0 + n,
        sigma ~ (1||n),
        nl = TRUE
        )

noise_priors <-
    c(
        #prior for single channel current given -60mV and 72 pS
        set_prior("normal(4.32, 1)", nlpar = "i"),
        #prior for number of channels given above current and observed max current
        set_prior("uniform(0, 100000)", nlpar = "nchannels", lb = "0", ub="100000"),
        set_prior("cauchy(0, 1)", nlpar = "i", class = "sd"),
        set_prior("cauchy(0, 1)", dpar = "sigma", class = "Intercept"),
        set_prior("cauchy(0, 1)", dpar = "sigma", class = "sd")
    )

noise_formula_2 <-
    bf(
        corrected_variance ~ (4.32 * corrected_pA) - ((corrected_pA ^ 2) / nchannels),
        nchannels ~ 0 + n,
        sigma ~ (1||n),
        nl = TRUE
        )

noise_priors_2 <-
    c(
        #prior for number of channels given above current and observed max current
        set_prior("uniform(0, 100000)", nlpar = "nchannels", lb = "0", ub="100000"),
        set_prior("cauchy(0, 1)", dpar = "sigma", class = "Intercept"),
        set_prior("cauchy(0, 1)", dpar = "sigma", class = "sd")
    )

noise_formula_3 <-
    bf(
        corrected_variance ~ 0 + corrected_pA + (corrected_pA||n)
        )

noise_priors_3 <-
    c(
        #prior for number of channels given above current and observed max current
        set_prior("normal(2, 1)", class = "b", lb = "0", ub = "4.32"),
        set_prior("cauchy(0, 1)", class = "sigma"),
        set_prior("cauchy(0, 1)", class = "sd")
    )

brms_iter <- 2000
brms_warmup <- 1000
brms_chains <- 4
brms_thin <- 1
brms_seed <- 2021

brm(
    formula = noise_formula,
    prior = noise_priors,
    data = cleaned_data,
    chains = brms_chains,
    iter = brms_iter,
    warmup = brms_warmup,
    thin  = brms_thin,
    seed = brms_seed,
    control = list(adapt_delta = 0.99, max_treedepth = 15),
    cores = getOption("mc.cores", 4),
    file = "data/other_fits/control_noise_analysis",
    sample_prior = "yes",
    save_all_pars = TRUE
    ) -> fit_i

brm(
    formula = noise_formula_2,
    prior = noise_priors_2,
    data = cleaned_data,
    chains = brms_chains,
    iter = brms_iter,
    warmup = brms_warmup,
    thin  = brms_thin,
    seed = brms_seed,
    control = list(adapt_delta = 0.99, max_treedepth = 15),
    cores = getOption("mc.cores", 4),
    file = "data/other_fits/control_noise_analysis_fixed",
    sample_prior = "yes",
    save_all_pars = TRUE
    ) -> fit_fixed

brm(
    formula = noise_formula_3,
    prior = noise_priors_3,
    data = cleaned_data,
    chains = brms_chains,
    iter = brms_iter,
    warmup = brms_warmup,
    thin  = brms_thin,
    seed = brms_seed,
    control = list(adapt_delta = 0.99, max_treedepth = 15),
    cores = getOption("mc.cores", 4),
    file = "data/other_fits/control_noise_analysis_popen",
    sample_prior = "yes",
    save_all_pars = TRUE
    ) -> fit_popen

cleaned_data %>%
ungroup() %>%
expand(nesting(n), corrected_pA = seq(min(corrected_pA), max(corrected_pA), length.out = 51)) %>%
add_fitted_draws(fit_fixed, re_formula = NA) %>%
median_qi(.value, .width = .95) %>%
mutate(model = "fixed") -> inferred_underlying_1

cleaned_data %>%
ungroup() %>%
expand(nesting(n), corrected_pA = seq(min(corrected_pA), max(corrected_pA), length.out = 51)) %>%
add_fitted_draws(fit_i, re_formula = NA) %>%
median_qi(.value, .width = .95) %>%
mutate(model = "free_i") -> inferred_underlying_2

cleaned_data %>%
ungroup() %>%
expand(nesting(n), corrected_pA = seq(min(corrected_pA), max(corrected_pA), length.out = 51)) %>%
add_fitted_draws(fit_popen) %>%
median_qi(.value, .width = .95) %>%
mutate(model = "popen") -> inferred_underlying_3

inferred_underlying <-
    bind_rows(inferred_underlying_1, inferred_underlying_2, inferred_underlying_3)

ggplot() +
geom_ribbon(data = inferred_underlying, aes(x = corrected_pA, ymin = .lower, ymax = .upper, colour = model), fill = NA) +
geom_point(data = cleaned_data, aes(x = corrected_pA, y = corrected_variance)) +
facet_wrap(vars(n)) +
coord_cartesian(ylim = c(0, 12000))


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
    mutate(augmented = map(fit, augment, newdata = tibble(pA = seq(0, 7500, length.out = 51))))

fixed <-
    wt_cor %>%
    group_by(n) %>%
    nest() %>%
    mutate(
        fit = map(data, ~ nls.multstart::nls_multstart(
            formula = variance ~ (4.32 * pA) - ((pA ^ 2) / nchannels),
            data = data.frame(.),
            iter=500,
            start_lower = list(nchannels = 100),
            start_upper = list(nchannels = 100000),
            control = nls.control(warnOnly = TRUE)
            ))) %>%
    mutate(augmented = map(fit, augment, newdata = tibble(pA = seq(0, 7500, length.out = 51))))

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
