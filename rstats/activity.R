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
