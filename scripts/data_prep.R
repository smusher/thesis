library(tidyverse)

import <- function(directory, pattern) {

    oldwd <- getwd()

    setwd(directory)

    filenames <-
        list.files(pattern = pattern)

    tempdf <-
        filenames %>%
        plyr::llply(
            read_csv, col_types = cols(
                construct = col_character(),
                date = col_date(),
                cumulative_exposure = col_integer(),
                concentration = col_double(),
                dye = col_double(),
                response = col_double(),
                raw_current = col_double(),
                wash_current = col_double(),
                barium_current = col_double(),
                raw_intensity = col_double(),
                corrected_intensity = col_double()
            )
        ) %>%
        bind_rows(.id = "n") %>%
        type_convert()

    setwd(oldwd)

    return(tempdf)
}

w311_current <-
    import("/home/sam/current_analysis/paper/pcf/W311-GFP+SUR/reanalyse", "currents") %>%
    mutate(
        construct = "W311*-GFP+SUR",
        measure = "current",
        )

w311c166s_current <-
    import("/home/sam/current_analysis/paper/pcf/W311-C166S-GFP+SUR/reanalyse", "currents") %>%
    mutate(
        construct = "W311*,C166S-GFP+SUR",
        measure = "current"
        )

w311_surk205e_current <-
    import("/home/sam/current_analysis/paper/pcf/W311-GFP+SUR-K205E", "currents") %>%
    mutate(
        construct = "W311*-GFP+SUR-K205E",
        measure = "current"
        )

w311_surk205a_current <-
    import("/home/sam/current_analysis/paper/pcf/W311-GFP+SUR-K205A", "currents") %>%
    mutate(
        construct = "W311*-GFP+SUR-K205A",
        measure = "current"
        )

w311e179a_current <-
    import("/home/sam/current_analysis/paper/pcf/W311-E179A-GFP+SUR", "currents") %>%
    mutate(
        construct = "W311*,E179A-GFP+SUR",
        measure = "current"
        )

w311e179k_current <-
    import("/home/sam/current_analysis/paper/pcf/W311-E179K-GFP+SUR", "currents") %>%
    mutate(
        construct = "W311*,E179K-GFP+SUR",
        measure = "current"
        )

w311k39r_current <-
    import("/home/sam/current_analysis/paper/pcf/W311-K39R-GFP+SUR", "currents") %>%
    mutate(
        construct = "W311*,K39R-GFP+SUR",
        measure = "current"
        )

w311k39e_current <-
    import("/home/sam/current_analysis/paper/pcf/W311-K39E-GFP+SUR", "currents") %>%
    mutate(
        construct = "W311*,K39E-GFP+SUR",
        measure = "current"
        )

w311k39a_current <-
    import("/home/sam/current_analysis/paper/pcf/W311-K39A-GFP+SUR", "currents") %>%
    mutate(
        construct = "W311*,K39A-GFP+SUR",
        measure = "current"
        )

w311mg_current <-
    import("/home/sam/current_analysis/paper/pcf/W311-GFP+SUR1+MG/mgbath", "currents") %>%
    mutate(
        construct = "W311*-GFP+SUR+MG",
        measure = "current"
        )

w311_fluorescence <-
    import("/home/sam/current_analysis/paper/pcf/W311-GFP+SUR/reanalyse", "intensities") %>%
    mutate(
        construct = "W311*-GFP+SUR",
        measure = "fluorescence",
        response = 1 - response
        )

w311c166s_fluorescence <-
    import("/home/sam/current_analysis/paper/pcf/W311-C166S-GFP+SUR/reanalyse", "intensities") %>%
    mutate(
        construct = "W311*,C166S-GFP+SUR",
        measure = "fluorescence",
        response = 1 - response
        )

w311_surk205e_fluorescence <-
    import("/home/sam/current_analysis/paper/pcf/W311-GFP+SUR-K205E", "intensities") %>%
    mutate(
        construct = "W311*-GFP+SUR-K205E",
        measure = "fluorescence",
        response = 1 - response
        )

w311_surk205a_fluorescence <-
    import("/home/sam/current_analysis/paper/pcf/W311-GFP+SUR-K205A", "intensities") %>%
    mutate(
        construct = "W311*-GFP+SUR-K205A",
        measure = "fluorescence",
        response = 1 - response
        )

w311e179a_fluorescence <-
    import("/home/sam/current_analysis/paper/pcf/W311-E179A-GFP+SUR", "intensities") %>%
    mutate(
        construct = "W311*,E179A-GFP+SUR",
        measure = "fluorescence",
        response = 1 - response
        )

w311e179k_fluorescence <-
    import("/home/sam/current_analysis/paper/pcf/W311-E179K-GFP+SUR", "intensities") %>%
    mutate(
        construct = "W311*,E179K-GFP+SUR",
        measure = "fluorescence",
        response = 1 - response
        )

w311k39r_fluorescence <-
    import("/home/sam/current_analysis/paper/pcf/W311-K39R-GFP+SUR", "intensities") %>%
    mutate(
        construct = "W311*,K39R-GFP+SUR",
        measure = "fluorescence",
        response = 1 - response
        )

w311k39e_fluorescence <-
    import("/home/sam/current_analysis/paper/pcf/W311-K39E-GFP+SUR", "intensities") %>%
    mutate(
        construct = "W311*,K39E-GFP+SUR",
        measure = "fluorescence",
        response = 1 - response
        )

w311k39a_fluorescence <-
    import("/home/sam/current_analysis/paper/pcf/W311-K39A-GFP+SUR", "intensities") %>%
    mutate(
        construct = "W311*,K39A-GFP+SUR",
        measure = "fluorescence",
        response = 1 - response
        )

w311mg_fluorescence <-
    import("/home/sam/current_analysis/paper/pcf/W311-GFP+SUR1+MG/mgbath", "intensities") %>%
    mutate(
        construct = "W311*-GFP+SUR+MG",
        measure = "fluorescence",
        response = 1 - response
        )

combined_pcf <-
    bind_rows(
        w311_current, w311_fluorescence,
        w311c166s_current, w311c166s_fluorescence,
        w311e179a_current, w311e179a_fluorescence,
        w311e179k_current, w311e179k_fluorescence,
        w311k39r_current, w311k39r_fluorescence,
        w311k39a_current, w311k39a_fluorescence,
        w311k39e_current, w311k39e_fluorescence,
        w311_surk205e_current, w311_surk205e_fluorescence,
        w311_surk205a_current, w311_surk205a_fluorescence,
        w311mg_current, w311mg_fluorescence
        ) %>%
    mutate(
        construct = factor(construct)
        ) %>%
    group_by(construct, n) %>%
    mutate(within_method_id = factor(cur_group_id())) %>%
    ungroup() %>%
    filter(dye < 480 | is.na(dye)) %>%
    select(construct, within_method_id, concentration, measure, response) %>%
    mutate(
        method = "pcf",
        nucleotide = "TNP-ATP",
        experimenter = "Sam"
        )

data_1 <-
  import("/home/sam/current_analysis/paper/tnp-atp_doseresponse_unroofed/W311-GFP+SUR", "intensities") %>%
  mutate(construct = "W311*-GFP+SUR")
data_2 <-
  import("/home/sam/current_analysis/paper/tnp-atp_doseresponse_unroofed/W311-GFP", "intensities") %>%
  mutate(construct = "W311*-GFP")
data_3 <-
  import("/home/sam/current_analysis/paper/tnp-atp_doseresponse_unroofed/W311-G334D-GFP+SUR", "intensities") %>%
  mutate(construct = "W311*,G334D-GFP+SUR")
data_5 <-
  import("/home/sam/current_analysis/paper/tnp-atp_doseresponse_unroofed/W311-C166S-GFP+SUR", "intensities") %>%
  mutate(construct = "W311*,C166S-GFP+SUR")
data_6 <-
  import("/home/sam/current_analysis/paper/tnp-atp_doseresponse_unroofed/W311-GFP+TMD0_232", "intensities") %>%
  mutate(construct = "W311*-GFP+TMD0_232")
data_7 <-
  import("/home/sam/current_analysis/paper/tnp-atp_doseresponse_unroofed/W311-GFP+TMD0_195", "intensities") %>%
  mutate(construct = "W311*-GFP+TMD0_195")
data_8 <-
  import("/home/sam/current_analysis/paper/tnp-atp_doseresponse_unroofed/W311-GFP+SUR-K205E", "intensities") %>%
  mutate(construct = "W311*-GFP+SUR-K205E")
data_9 <-
  import("/home/sam/current_analysis/paper/tnp-atp_doseresponse_unroofed/W311+SUR", "intensities") %>%
  mutate(construct = "W311*+SUR")
data_10 <-
  import("/home/sam/current_analysis/paper/tnp-atp_doseresponse_unroofed/W311-GFP+SUR+MG", "intensities") %>%
  mutate(construct = "W311*-GFP+SUR+MG")

combined_unroofed <-
    bind_rows(data_1, data_2, data_3, data_5, data_6, data_7, data_8, data_9, data_10) %>%
    mutate(
      construct = factor(construct),
      response = 1 - response
      ) %>%
    group_by(construct, n) %>%
    mutate(within_method_id = factor(cur_group_id())) %>%
    ungroup() %>%
    filter(dye < 480) %>%
    select(construct, within_method_id, concentration, response) %>%
    mutate(
        method = "unroofed",
        measure = "fluorescence",
        nucleotide = "TNP-ATP",
        experimenter = "Sam"
        )

patch_1 <-
    read_csv("data/electrophys_data.csv") %>%
    select(n, method, measure, construct, nucleotide, concentration, response) %>%
    mutate(experimenter = "Sam")

patch_2 <-
    read_csv("data/natascias_data.csv") %>%
    select(n, method, measure, construct, nucleotide, concentration, response) %>%
    mutate(experimenter = "Natascia")

combined_patch <-
    bind_rows(patch_1, patch_2) %>%
    group_by(experimenter, n, method, measure, construct, nucleotide) %>%
    mutate(within_method_id = factor(cur_group_id())) %>%
    ungroup() %>%
    select(-n)

combined_data <-
    bind_rows(combined_pcf, combined_patch, combined_unroofed) %>%
    group_by(method, within_method_id) %>%
    mutate(unique_experiment_id = factor(cur_group_id())) %>%
    write_csv("data/combined_drc_data.csv")

# ggplot(combined_data %>% unite(experiment, c(method, measure)), aes(x = concentration, y = response, colour = experiment)) +
# geom_point() +
# facet_wrap(vars(unique_experiment_id)) +
# scale_x_log10() +
# theme_void() +
# theme(legend.position = "top") +
# scale_colour_brewer(palette = "Set2")
