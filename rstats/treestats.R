library(tidyverse)
library(ape)
# devtools::install_github("https://github.com/heavywatal/aptree")
# wtl::refresh("aptree")
library(aptree)

ggplot_tetree = function(data, colorbar=TRUE, hypercolor = NA) {
  .max_copy_number = max(data[["copy_number"]], na.rm = TRUE)
  .size_breaks = c(32, 16, 8, 4, 2, 1) %>% {.[. < .max_copy_number]}
  .colors = rev(head(rainbow(15L), 12L))
  .values = seq(0, 1, length.out = length(.colors))
  .breaks = c(0, 0.5, 1)
  if (!is.na(hypercolor)) {
    .colors = c(.colors, hypercolor)
    .breaks = c(.breaks, 2)
    .values = c(seq(0, 0.5, length.out = length(.colors)), 1)
  }
  .limits = c(0, max(.breaks))
  .ylim = data$yend %>% {c(min(.), max(.) * 1.1)}
  .size_guide = guide_legend(order = 10, override.aes = list(stroke = 0, alpha = 0.4))
  .oaes = list(alpha = 0.4, fill = "#000000", stroke = 1, size = 5)
  .pch_guide = guide_legend(order = 20, title = "Frequency", override.aes = .oaes)
  .col_guide = if (colorbar) {
    guide_colourbar(order = 30, title = "Acitivity\nLevel", barheight = 12, reverse = TRUE)
  } else {FALSE}
  filter_active = function(x) dplyr::filter(x, activity > 0)
  filter_major = function(x) dplyr::filter(x, is_major)
  ggplot(data %>% dplyr::mutate(is_major = dplyr::coalesce(is_major, FALSE))) +
  geom_point(data = filter_active, aes(xend, yend, colour = activity, size = copy_number), pch = 16, alpha = 0.6, stroke = 1) +
  geom_point(data = filter_major, aes(xend, yend, size = copy_number), pch = 1, colour = "#000000", alpha = 0.5, stroke = 1) +
  geom_point(data = filter_major, aes(xend, yend, colour = activity, size = copy_number), pch = 1, alpha = 0.6, stroke = 1) +
  geom_point(aes(xend, yend, shape = is_major), alpha = 0) +
  geom_segment(aes(x, y, xend = xend, yend = yend), size = 0.28) +
  scale_shape_manual(values = c(`TRUE` = 21, `FALSE` = 16), guide = .pch_guide, labels = c("< 0.5", "≥ 0.5")) +
  scale_colour_gradientn(colours = .colors, limits = .limits, breaks = .breaks, values = .values, guide = .col_guide) +
  scale_size(limit = c(1, .max_copy_number), breaks = .size_breaks, range = c(3, 12), name = "Copy\nNumber", guide = .size_guide) +
  labs(x = NULL, y = NULL) +
  coord_fixed(ylim = .ylim)
}

label_both_tree = function(labels, multi_line = TRUE) {
  value <- label_value(labels, multi_line = multi_line)
  variable <- labels %>% dplyr::rename(t = generation) %>% names()
  out <- vector("list", length(value))
  for (i in seq_along(out)) {
      out[[i]] <- paste(variable[[i]], value[[i]], sep = " = ")
  }
  out
}

add_phylo = function(.tbl, root=NULL) {
  .tbl %>%
    dplyr::mutate(distmat = purrr::map(seqs, Biostrings::stringDist, method = "hamming")) %>%
    dplyr::mutate(phylo = purrr::map(distmat, ~{
      if (length(.x) > 2) {
        .phy = ape::fastme.ols(.x)
        if (!is.null(root) && (root %in% .phy$tip.label)) {
          .phy = ape::root(.phy, root, resolve.root = TRUE)
        }
        .phy
      } else {NA}
    }))
}

nest_metadata = function(.tidy_metadata, dss) {
  .useqs = unique_dss(dss)
  .tidy_metadata %>%
    tidyr::nest(-individual) %>%
    dplyr::mutate(seqs = purrr::map(data, ~{.useqs[.x$label]})) %>%
    add_phylo(root = origin_name)
}

eval_treeshape = function(.tblphy) {
  dplyr::transmute(.tblphy,
    generation,
    shape_pda = purrr::map_dbl(phylo, ~{
      if (identical(.x, NA)) return(NA)
      aptree::as.treeshape(.x, "pda") %>%
      aptree::shape.statistic(norm = "pda")
    }),
    bimodality = purrr::map_dbl(distmat, wtl::bimodality)
  )
}

make_range_stats = function(generation=0L) {
  tibble::tibble(
    generation,
    stat = rep(c("bimodality", "shape_pda"), each = 2L),
    value = c(0, 1, 2, -10))
}

label_value_tr = function(labels, multi_line = TRUE) {
  lapply(labels, function(x) {c(bimodality = "BI", shape_pda = "Shape stat. (PDA)")[x]})
}

ggplot_evolution = function(.tblstats, only_bi = FALSE) {
  .x = tidyr::gather(.tblstats, stat, value, -generation)
  .blank = make_range_stats(min(.x$generation))
  .hline = tibble::tibble(stat = c("bimodality", "shape_pda"), value = c(5 / 9, NA))
  if (only_bi) {
    .x = .x %>% dplyr::filter(stat == "bimodality")
    .blank = .blank %>% dplyr::filter(stat == "bimodality")
    .hline = .hline %>% dplyr::filter(stat == "bimodality")
  }
  ggplot(.x, aes(generation, value)) +
  geom_blank(data = .blank) +
  geom_line(size = 0.8) +
  geom_hline(data = .hline, aes(yintercept = value), colour = "red", linetype = "dashed") +
  facet_grid(stat ~ ., scale = "free_y", switch = "y", label = label_value_tr) +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    strip.placement = "outside",
    strip.background = element_blank(),
    strip.text = element_text(size = 12, family = "Helvetica"),
  )
}
