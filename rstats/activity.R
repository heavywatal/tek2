library(tidyverse)

read_activity = function(indir) {
  file.path(indir, "activity.tsv.gz") %>% read_tsv()
}

ggplot_activity = function(data, popsize=1, hypercolor = NA) {
  .guide = guide_colorbar(title = "Activity\nLevel", barheight = 12, reverse = TRUE)
  .colors = rev(head(rainbow(15L), 12L))
  .values = seq(0, 1, length.out = length(.colors))
  .breaks = c(0, 0.5, 1)
  if (!is.na(hypercolor)) {
    .colors = c(.colors, hypercolor)
    .breaks = c(.breaks, 2)
    .values = c(seq(0, 0.5, length.out = length(.colors)), 1)
  }
  .limits = c(0, max(.breaks))
  dplyr::mutate(data,
    copy_number = copy_number / popsize,
    species = factorize_species(species)
  ) %>%
  ggplot(aes(generation, copy_number)) +
  geom_area(aes(group = interaction(activity, species), fill = activity), position = position_stack(reverse = FALSE)) +
  scale_fill_gradientn(colours = .colors, limits = .limits, values = .values, breaks = .breaks, guide = .guide) +
  labs(x = "Generation", y = "Average Copy Number") +
  wtl::theme_wtl()
}

plot_copynumber_generation = function(data) {
  ggplot(data, aes(generation, copy_number, group = activity)) +
    geom_area(aes(fill = activity)) +
    scale_fill_gradientn(colours = rev(head(rainbow(15L), 12L)), limits = c(0, 1), breaks = c(0, 0.5, 1)) +
    wtl::theme_wtl() +
    theme(legend.position = "bottom")
}
# .counted %>% plot_copynumber_generation()

factorize_species = function(x) {
  factor(x, levels = sort.int(unique(x), decreasing = TRUE))
}


# #######1#########2#########3#########4#########5#########6#########7#########
if (FALSE) {

.tbl_act = .metadata %>%
  # dplyr::filter(!(xi < 6e-4 & upper < 16)) %>%
  # dplyr::filter(lower < 300, upper < 300) %>%
  dplyr::filter(lower == 300L, upper == 300L) %>%
  # sample_n(6L) %>%
  dplyr::mutate(adata = purrr::map(indir, read_activity)) %>%
  tidyr::unnest() %>%
  print()

.p = ggplot_activity(.tbl_act) +
  facet_grid(xi * lower ~ upper * repl)
.p
ggsave("copynumber-activity-species.png", .p, width = 10, height = 10)

.p = .tbl_act %>%
  dplyr::distinct(lower, upper, repl, generation, species) %>%
  dplyr::count(lower, upper, repl, species) %>%
  dplyr::mutate(n = n * 50L) %>%
  ggplot(aes(species, n)) +
  geom_col() +
  facet_grid(lower ~ upper * repl) +
  labs(x = "species ID", y = "sojourn time") +
  wtl::theme_wtl()
.p
ggsave("species-sojourn-time.png", .p, width = 7, height = 7)

.nested = .tbl_act %>%
  tidyr::nest(-xi) %>%
  dplyr::mutate(plt = purrr::map2(data, paste0("xi = ", xi), ~{
    plot_copynumber_generation(.x) +
      facet_grid(alpha + desc(nu) ~ lambda, labeller = label_both) +
      labs(title = .y)
  }))

.p = cowplot::plot_grid(plotlist = .nested$plt, nrow = 1)
.p
.nested$plt[[1]] + coord_cartesian(ylim = c(0, 370))

ggsave("fig-s2.png", .p, width = 15, height = 12)

} # if (FALSE)
