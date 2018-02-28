library(tidyverse)
library(jsonlite)

.cols = c("site", "species", "indel", "nonsynonymous", "synonymous", "activity")
.nsam = 100L
.gametes = "summary.json.gz" %>% fromJSON()

.summary = .gametes %>%
  sample(.nsam) %>%
  {
    tibble(tmpcol = ., gamete = seq_along(.))
  } %>%
  tidyr::unnest() %>%
  tidyr::separate(tmpcol, .cols, ":") %>%
  dplyr::mutate_at(vars(c("site", "species", "indel", "nonsynonymous", "synonymous")), as.integer) %>%
  dplyr::mutate(activity = as.double(activity)) %>%
  dplyr::mutate(dn = nonsynonymous / 200, ds = synonymous / 100) %>%
  dplyr::mutate(dn_ds = dn / ds, distance = (nonsynonymous + synonymous) / 300) %>%
  print()

.summary %>%
  ggplot(aes(dn_ds, group = activity)) +
  geom_histogram(aes(fill = activity), binwidth = 0.1, center = 0) +
  scale_fill_gradientn(colours = rev(head(rainbow(15L), 12L)), breaks = c(0, 0.5, 1)) +
  labs(x = "dN/dS", y = "Copy number") +
  theme_bw()

.summary %>%
  ggplot(aes(distance, dn_ds)) +
  geom_jitter(aes(colour = activity), alpha = 0.3) +
  scale_colour_gradientn(colours = rev(head(rainbow(15L), 12L)), breaks = c(0, 0.5, 1)) +
  geom_hline(aes(yintercept = 1.0), linetype = "dashed") +
  coord_cartesian(xlim = c(0, 1)) +
  labs(x = "Distance from original sequence", y = "dN/dS") +
  theme_bw()

.summary %>%
  dplyr::count(site) %>%
  dplyr::mutate(freq = n / .nsam) %>%
  ggplot(aes(freq)) +
  geom_histogram(binwidth = 0.1, center = 0) +
  coord_cartesian(xlim = c(0, 1)) +
  labs(x = "TE frequency", y = "Site number") +
  theme_bw()
