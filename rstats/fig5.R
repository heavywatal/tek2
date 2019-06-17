if (FALSE) {

fig5acti = .metadata %>%
  dplyr::mutate(
    title = sprintf("n=%d xi=%.0e c=%d lower=%d upper=%d (%d)", n, xi, coexist, lower, upper, repl),
    adata = purrr::map(indir, read_activity),
    aplot = purrr::map2(adata, n, ggplot_activity),
    aplot = purrr::map2(aplot, title, ~.x + labs(title = .y) + theme(legend.position = "none"))
  ) %>%
  print()
# fig5acti$aplot[[1]]
ggsave("fig5_activity.pdf", fig5acti$aplot, width = 8, height = 4)

.fig5actidy = .metadata %>%
  dplyr::mutate(
    adata = purrr::map(indir, read_activity),
    indir = NULL
  ) %>%
  tidyr::unnest() %>%
  dplyr::mutate(copy_number = copy_number / n) %>%
  print()

.p5facet = .fig5actidy %>%
  dplyr::filter(xi == 0.001) %>%
  dplyr::filter(generation %% 1000L == 0L) %>%
  ggplot_activity() +
  facet_grid(coexist + n + lower ~ upper + repl, scale = "free_y")
ggsave("fig5facet_possible.pdf", .p5facet, width = 30, height = 20)
# ggsave("fig5facet.pdf", .p5facet, width=30, height=20)

.p5facet_narrow = .fig5actidy %>%
  dplyr::filter(xi == 0.001, n == 500L, upper < 30L) %>%
  dplyr::filter(generation %% 1000L == 0L) %>%
  ggplot_activity() +
  facet_grid(coexist + lower ~ upper + repl, scale = "free_y")
ggsave("fig5facet_narrow.pdf", .p5facet_narrow, width = 20, height = 12)

}
# #######1#########2#########3#########4#########5#########6#########7#########

source("~/git/teaposon/rstats/read.R")
.indirs = fs::dir_ls(regexp = "^n\\d+", type = "directory")
.metadata = read_metadata(.indirs) %>% print()

fig5df = .metadata %>%
  dplyr::filter(coexist == 2L, xi == 0.001, n == 1000L, lower == 6L, upper == 24L, repl == 3L) %>%
  dplyr::mutate(
    adata = purrr::map(indir, read_activity),
    tdata = purrr::map(indir, ~{
      message(.x)
      read_fastas(.x, interval = 100L) %>%
        add_phylo() %>%
        eval_treeshape()
    })
  ) %>%
  print()

saveRDS(fig5df, "fig5.rds")

.gens5 = c(17000L, 22000L, 27000L, 30000L, 40000L)
.infiles5 = fs::path(fig5df$indir, sprintf("generation_%05d.fa.gz", .gens5))
.inds_df5 = .infiles5 %>%
  setNames(str_extract(., "(?<=generation_)\\d+")) %>%
  purrr::map_dfr(read_individuals, .id = "generation") %>%
  dplyr::mutate(generation = as.integer(generation)) %>%
  dplyr::filter(individual != "total") %>%
  print()

.timepoints_df = fig5df$tdata[[1]] %>%
  dplyr::filter(generation %in% .gens5) %>%
  dplyr::transmute(generation, value = bimodality, stat = "bimodality") %>%
  print()

fig5plts = fig5df %>%
  dplyr::mutate(
    tplot = purrr::map(tdata, ~{
      .xmin = min(.x$generation)
      ggplot_evolution(.x %>% dplyr::filter(generation > 500L), only_bi = TRUE) +
        coord_cartesian(xlim = c(.xmin, max(.x$generation))) +
        geom_point(data = .timepoints_df, colour = "red", size = 3) +
        theme(
          axis.title = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank()
        )
    }),
    aplot = purrr::map2(adata, n, ~{
      .p = ggplot_activity(.x, .y)
      if (max(.x$copy_number) / .y < 100) {
        .p = .p + theme(axis.title.y = element_text(margin = margin(r = 11.5)))
      }
      .p
    }),
    top = purrr::map2(tplot, aplot, ~{
      cowplot::plot_grid(.x, .y + theme(legend.position = "none"),
        ncol = 1, rel_heights = c(1, 1.2), align = "v", axis = "lr")
    })
  ) %>%
  print()
fig5plts$top[[1]]

.df_unrooted5 = .inds_df5 %>%
  purrr::pmap_df(function(phylo, data, generation, individual, ...) {
    wtl::ape_layout_unrooted(phylo) %>%
      dplyr::left_join(data, by = "label") %>%
      dplyr::mutate(generation = generation, individual = individual)
  }) %>%
  print()

.p = .df_unrooted5 %>%
  ggplot_tetree() +
  facet_grid(generation ~ individual) +
  theme_classic()
.p
ggsave("fig5_tree_candidates.pdf", .p, width = 9.9, height = 7, scale = 2, device = cairo_pdf)

.inds = c("0", "0", "5", "2", "7")
.delegates = .timepoints_df %>%
  dplyr::transmute(generation, BI = sprintf("%.3f", value), individual = .inds) %>%
  print()
.delegates_df = .delegates %>%
  dplyr::left_join(.df_unrooted5 %>% tidyr::nest(-generation, -individual), by = c("generation", "individual")) %>%
  tidyr::unnest() %>%
  print()

.annot_base = .delegates_df %>% dplyr::distinct(generation, BI) %>% head(1L) %>% print()
.dsegm = .annot_base %>% dplyr::mutate(x = 5, xend = 8, y = -10, yend = -10) %>% print()
.dtext = .annot_base %>% dplyr::mutate(x = 6.5, y = -12, label = "0.01") %>% print()

.fig5_bottom = ggplot_tetree(.delegates_df, colorbar = FALSE) +
  geom_segment(data = .dsegm, aes(x, y, xend = xend, yend = yend), size = 1) +
  geom_text(data = .dtext, aes(x, y, label = label)) +
  facet_wrap(~ generation + BI, labeller = label_both_tree, nrow = 1L) +
  theme_classic() +
  theme(
    axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(),
    strip.background = element_blank(),
    strip.text.x = element_text(size = 12, hjust = 0, margin = margin(6, 12, 3, 12)),
    legend.box.margin = margin(b = 30)
  )
.fig5_bottom
ggsave("fig5_bottom.pdf", .fig5_bottom, width = 10, height = 4, family = "Helvetica", device = cairo_pdf)

.guide_acti = guide_colorbar(title = "Activity\nLevel", barheight = 11, reverse = TRUE)
.cbar = cowplot::get_legend(fig5plts$aplot[[1L]] + guides(fill = .guide_acti))

fig5 = cowplot::plot_grid(
  fig5plts$top[[1L]],
  cowplot::plot_grid(.fig5_bottom, .cbar, nrow = 1L, rel_widths = c(1, 0.1)),
  ncol = 1L, rel_heights = c(1.3, 1), labels = c("A", "B"), scale = 0.95
)
fig5
ggsave("fig5.pdf", fig5, width = 10, height = 8, family = "Helvetica", device = cairo_pdf)
