library(tidyverse)
library(broom)

data_1 <-
	read_csv("data/confocal_images/per_pixel_intensities.csv") %>%
	select(-X1) %>%
	mutate(
		construct = factor(construct),
		image_stack = as.integer(image_stack)
		) %>%
	filter(
		gfp_intensity != max(gfp_intensity),
		anap_intensity != max(anap_intensity)
		) %>%
	group_by(construct, image_stack) %>%
	mutate(across(
		ends_with("intensity"),
		~(log(.x) - mean((log(.x)))) / sd(log(.x)),
		.names = "{.col}_norm"
		))

ggplot(data = data, aes(x = gfp_intensity_norm, y = anap_intensity_norm)) +
geom_point(size = 0.1, alpha = 0.01) +
geom_smooth(method = lm) +
facet_grid(rows = vars(construct))

linear <- 
	data %>%
	group_by(construct) %>%
	nest() %>%
	mutate(fit = map(data, ~lm(
		anap_intensity_norm ~ gfp_intensity_norm,
		data = data.frame(.)))
	)

linear %>%
mutate(tidied = map(fit, tidy)) %>%
unnest(tidied)

wt_data <-
	read_csv("data/confocal_images/wt-example_profile.csv") %>%
	select(-X1) %>%
	mutate(construct = "WT-GFP+SUR1")
w311_data <-
	read_csv("data/confocal_images/w311-example_profile.csv") %>%
	select(-X1) %>%
	mutate(construct = "W311*-GFP+SUR1")

data_2 <-
	bind_rows(wt_data, w311_data) %>%
	group_by(construct, channel) %>%
	mutate(
		intensity_norm = (intensity - min(intensity)) / (max(intensity) - min(intensity)),
		chunk = ntile(y, 200)
		) %>%
	group_by(construct, channel, chunk, x) %>%
	summarise(
		intensity = mean(intensity),
		intensity_norm = mean(intensity_norm)
		) %>%
	group_by(construct, channel, chunk) %>%
	mutate(intensity_norm = (intensity_norm - min(intensity_norm)) / (max(intensity_norm) - min(intensity_norm)))

chunk_sample <- sample.int(200, 10)

sampled_data <-
	data_2 %>%
	filter(chunk %in% chunk_sample, channel !=0)

ggplot(sampled_data, aes(x = x, y = intensity_norm, colour = factor(channel))) +
geom_line() +
facet_grid(rows = vars(chunk), cols = vars(construct)) +
scale_colour_brewer(palette = "Set1")
