if (FALSE) {

# @ilp15% ~/git/teaposon/run.py -p2 -j6 te2fig4 -r96

fig4cand = .metadata %>% dplyr::transmute(n, repl, lapply(indir, read_activity)) %>% print()
.p = fig4cand %>% tidyr::unnest() %>%
  ggplot_activity(1000, "#aa0000") +
  coord_cartesian(xlim = c(3000, 6000)) +
  facet_wrap(~n + repl)
.p
ggsave("fig4_candidates.png", .p, width = 16, height = 16)

}
# #######1#########2#########3#########4#########5#########6#########7#########

source("~/git/teaposon/rstats/read.R")
.indirs = fs::dir_ls(regexp = "^n\\d+", type = "directory")
.metadata = read_metadata(.indirs) %>% print()

.df4 = .metadata %>% dplyr::filter(n == 1000L, str_detect(indir, "_063$")) %>% print()

df4act = .df4 %>%
  dplyr::transmute(data = lapply(indir, read_activity)) %>%
  tidyr::unnest() %>%
  print()
df4act %>% ggplot_activity(1000, "#aa0000")

.gens4 = c(4000L, 4300L, 4600L, 4900L, 5200L)
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
  ggplot_tetree(hypercolor = "#aa0000") +
  facet_grid(generation ~ individual) +
  theme_classic()
.p
ggsave("fig4_tree_candidates.pdf", .p, width = 9.9, height = 7, scale = 2, device = cairo_pdf)

.timepoints_df = df4act %>%
  dplyr::filter(generation %in% .gens4) %>%
  dplyr::group_by(generation) %>%
  dplyr::summarise(copy_number = sum(copy_number) / 1000) %>%
  dplyr::mutate(label = tolower(as.roman(dplyr::row_number()))) %>%
  print()

fig4act = df4act %>%
  ggplot_activity(1000, "#aa0000") +
  geom_vline(xintercept = .gens4, linetype = "dashed") +
  # geom_point(data = .timepoints_df, position = position_nudge(0, 3), shape = 25, fill = "#000000") +
  # geom_text(data = .timepoints_df, aes(label = label), position = position_nudge(-40, 2)) +
  coord_cartesian(xlim = c(3700, 5500), ylim = c(0, 60), expand = FALSE) +
  scale_x_continuous(breaks = .gens4) +
  scale_y_continuous(breaks = c(0, 20, 40)) +
  theme(panel.grid.major = element_blank())
fig4act
ggsave("fig4_activity.pdf", fig4act, width = 4, height = 5, scale = 1)

.inds = c("0", "4", "4", "4", "4")
df4delegates = tibble::tibble(generation = .gens4, individual = .inds) %>%
  dplyr::left_join(.df4unrooted, by = c("generation", "individual")) %>%
  print()

df4delegates %>% dplyr::select(matches("^x|^y")) %>% dplyr::summarise_all(list(~min(.)))
.annot_base = df4delegates %>% dplyr::distinct(generation) %>% tail(1L) %>% print()
.dsegm = .annot_base %>% dplyr::mutate(x = 0, xend = 1.8, y = -2, yend = -2) %>% print()
.dtext = .annot_base %>% dplyr::mutate(x = 0.9, y = -2.3, label = "0.005") %>% print()

fig4trees = ggplot_tetree(df4delegates, colorbar = FALSE, hypercolor = "#aa0000") +
  geom_segment(data = .dsegm, aes(x, y, xend = xend, yend = yend), size = 1) +
  geom_text(data = .dtext, aes(x, y, label = label)) +
  facet_wrap(~ generation, labeller = label_both_tree, nrow = 1L) +
  theme_classic() +
  theme(
    axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(),
    strip.background = element_blank(),
    strip.text.x = element_text(size = 12, hjust = 0.2, vjust = 0)
  )
fig4trees
ggsave("fig4_trees.pdf", fig4trees, width = 12, height = 4, scale = 1, family = "Helvetica", device = cairo_pdf)

fig4 = cowplot::plot_grid(fig4act, fig4trees, nrow = 1L, rel_widths = c(1, 3))
fig4
ggsave("fig4.pdf", fig4, width = 15, height = 4, family = "Helvetica", device = cairo_pdf)
