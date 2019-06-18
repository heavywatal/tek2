if (FALSE) {

# @ilp15% ~/git/tek2/run.py -p2 -j6 te2fig4 -r96

fig4cand = .metadata %>% dplyr::transmute(n, repl, lapply(indir, read_activity)) %>% print()
.p = fig4cand %>% tidyr::unnest() %>%
  ggplot_activity(1000, "#aa0000") +
  coord_cartesian(xlim = c(3000, 6000)) +
  facet_wrap(~n + repl)
.p
ggsave("fig4_candidates.png", .p, width = 16, height = 16)

}
# #######1#########2#########3#########4#########5#########6#########7#########

source("~/git/tek2/rstats/read.R")
.indirs = fs::dir_ls(regexp = "^n\\d+", type = "directory")
.metadata = read_metadata(.indirs) %>% print()

.df4 = .metadata %>% dplyr::filter(n == 1000L, str_detect(indir, "_063$")) %>% print()

df4act = .df4$indir %>% read_activity() %>% print()
df4act %>% ggplot_activity(1000, "#aa0000")

.lbound = 3700L; .ubound = 5500L
.gens4 = c(4000L, 4300L, 4600L, 4900L, 5200L)

# .timepoints_df = df4act %>%
#   dplyr::filter(generation %in% .gens4) %>%
#   dplyr::group_by(generation) %>%
#   dplyr::summarise(copy_number = sum(copy_number) / 1000) %>%
#   dplyr::mutate(label = tolower(as.roman(dplyr::row_number()))) %>%
#   print()

fig4act = df4act %>%
  ggplot_activity(1000, "#660000") +
  geom_vline(xintercept = .gens4, linetype = "dashed", colour = "#666666") +
  # geom_point(data = .timepoints_df, position = position_nudge(0, 1), shape = 25, fill = "#000000") +
  # geom_text(data = .timepoints_df, aes(label = label), position = position_nudge(-40, 2)) +
  coord_cartesian(xlim = c(.lbound, .ubound), ylim = c(0, 60), expand = FALSE) +
  scale_x_continuous(breaks = .gens4) +
  scale_y_continuous(breaks = c(0, 20, 40)) +
  wtl::theme_wtl(base_size = 13) +
  theme(panel.grid.major = element_blank())
fig4act
ggsave("fig4_activity.pdf", fig4act, width = 4, height = 5, scale = 1)

df4stat = .df4$indir %>%
  read_fastas(interval = 100L, from = .lbound, .ubound) %>%
  add_phylo() %>%
  eval_treeshape() %>%
  print()

.timepoints_df = df4stat %>%
  dplyr::filter(generation %in% .gens4) %>%
  dplyr::transmute(generation, value = bimodality, stat = "bimodality") %>%
  print()

fig4stat = df4stat %>%
  ggplot_evolution(only_bi = TRUE) +
  geom_vline(aes(xintercept = generation), .timepoints_df, linetype = "dashed", colour = "#666666") +
  geom_point(data = .timepoints_df, colour = "red", size = 3) +
  coord_cartesian(xlim = c(.lbound, .ubound), expand = FALSE) +
  theme(
    panel.grid.major.x = element_blank(),
    axis.title = element_blank(), axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )
fig4stat

fig4left = cowplot::plot_grid(
  fig4stat, fig4act + theme(legend.position = "none"),
  ncol = 1L, rel_heights = c(1, 1.2), align = "v", axis = "lr"
) + coord_fixed(ratio = 1.8)
fig4left
ggsave("fig4_left.pdf", fig4left, width = 3, height = 7, scale = 1)

.df4inds = fs::path(.df4$indir, sprintf("generation_%05d.fa.gz", .gens4)) %>%
  setNames(str_extract(., "(?<=generation_)\\d+")) %>%
  purrr::map_dfr(read_individuals, .id = "generation") %>%
  dplyr::mutate(generation = as.integer(generation)) %>%
  dplyr::filter(individual != "total") %>%
  print()

.df4unrooted = .df4inds %>%
  purrr::pmap_df(function(phylo, data, generation, individual, ...) {
    wtl::ape_layout_unrooted(phylo) %>%
      dplyr::left_join(data, by = "label") %>%
      dplyr::mutate(generation = generation, individual = individual)
  }) %>%
  print()

.p = .df4unrooted %>%
  ggplot_tetree(hypercolor = "#660000") +
  facet_grid(generation ~ individual) +
  theme_classic()
.p
ggsave("fig4_tree_candidates.pdf", .p, width = 9.9, height = 7, scale = 2, device = cairo_pdf)

.inds = c("0", "4", "4", "4", "4")
df4delegates = tibble::tibble(generation = .gens4, individual = .inds) %>%
  dplyr::left_join(.df4unrooted, by = c("generation", "individual")) %>%
  print()

df4delegates %>% dplyr::select(matches("^x|^y")) %>% dplyr::summarise_all(list(~min(.)))
.annot_base = df4delegates %>% dplyr::distinct(generation) %>% tail(1L) %>% print()
.dsegm = .annot_base %>% dplyr::mutate(x = 0, xend = 1.8, y = -2, yend = -2) %>% print()
.dtext = .annot_base %>% dplyr::mutate(x = 0.9, y = -2.3, label = "0.005") %>% print()

fig4trees = ggplot_tetree(df4delegates, colorbar = TRUE, hypercolor = "#660000") +
  geom_segment(data = .dsegm, aes(x, y, xend = xend, yend = yend), size = 1) +
  geom_text(data = .dtext, aes(x, y, label = label)) +
  facet_wrap(~ generation, labeller = label_both_tree, nrow = 2L) +
  theme_classic(base_size = 13) +
  theme(
    axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(),
    panel.spacing.y = grid::unit(4, "lines"),
    strip.background = element_blank(),
    strip.text.x = element_text(size = 13, hjust = 0.2, vjust = 0)
  )
fig4trees
ggsave("fig4_trees.pdf", fig4trees, width = 9, height = 7, scale = 1, family = "Helvetica", device = cairo_pdf)

fig4 = cowplot::plot_grid(fig4left, fig4trees, nrow = 1L, rel_widths = c(1, 3))
fig4
ggsave("fig4.pdf", fig4, width = 12, height = 7, family = "Helvetica", device = cairo_pdf)
