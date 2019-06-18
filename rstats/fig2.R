if (FALSE) {

# @ILP23% ~/git/tek2/run.py -p2 -j6 te2fig2 -r96

fig2cand = .metadata %>%
  dplyr::filter(xi > 6e-4) %>%
  dplyr::mutate(
    adata = purrr::map(indir, read_activity),
    tdata = purrr::map(indir, ~{
      message(.x)
      read_fastas(.x, interval = 500L) %>%
        add_phylo() %>%
        eval_treeshape()
    })
  ) %>%
  print()

fig2cand_plts = fig2cand %>%
  dplyr::mutate(
    title = sprintf("n=%d xi=%.0e (%d)", n, xi, repl),
    aplot = purrr::map2(adata, n, ggplot_activity),
    tplot = purrr::map(tdata, ~{
      ggplot_evolution(.x, only_bi = TRUE) +
        theme(
          axis.title = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank()
        )
    }),
  ) %>%
  purrr::pmap(function(aplot, tplot, title, ...) {
    .title = cowplot::ggdraw() + cowplot::draw_label(title, x = 0.1, hjust = 0)
    aplot = aplot + theme(legend.position = "none")
    cowplot::plot_grid(.title, tplot, aplot, ncol = 1, rel_heights = c(0.2, 1, 1), align = "v", axis = "lr")
  })
ggsave("fig2_candidates.pdf", fig2cand_plts, width = 8, height = 4)

# rsync -auv ilp23.local:~/working/te2-fig2/ ~/working/tek/te2-fig2/

fig2windows = tibble::tribble(
  ~repl, ~lbound, ~ubound,
    43L,  25000L,  45000L,
    46L,  10000L,  30000L,
    47L,  15000L,  30000L,
    65L,  10000L,  30000L,
    73L,  10000L,  27000L
) %>% print()

fig2candidates = .metadata %>%
  dplyr::right_join(fig2windows, by = "repl") %>%
  dplyr::mutate(
    adata = purrr::pmap(., function(indir, lbound, ubound, ...) {
      read_activity(indir) %>% dplyr::filter(lbound <= generation, generation <= ubound)
    }),
    tdata = purrr::pmap(., function(indir, lbound, ubound, ...) {
      message(indir)
      read_fastas(indir, interval = 100L, from = lbound, to = ubound) %>%
        add_phylo() %>%
        eval_treeshape()
    })
  ) %>%
  print()

fig2candplt = fig2candidates %>% dplyr::transmute(
  aplot = purrr::map2(adata, n, ggplot_activity),
  tplot = purrr::map(tdata, ~{
    ggplot_evolution(.x, only_bi = TRUE) +
      theme(
        axis.title = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()
      )
  }),
  plt = purrr::pmap(list(aplot, tplot, repl), function(.x, .y, .z) {
    .title = cowplot::ggdraw() + cowplot::draw_label(sprintf("repl=%d", .z), x = 0.1, hjust = 0)
    .x = .x + theme(legend.position = "none")
    cowplot::plot_grid(.title, .y, .x, ncol = 1, rel_heights = c(0.2, 1, 1), align = "v", axis = "lr")
  })
)
# ggsave("fig2left.pdf", fig2candplt$plt, width=4, height=9, family="Helvetica")
.pg = cowplot::plot_grid(plotlist = fig2candplt$plt, ncol = 3L)
ggsave("fig2_finalists.pdf", .pg, width = 7, height = 9.9, family = "Helvetica", scale = 1.5)

}
# #######1#########2#########3#########4#########5#########6#########7#########

source("~/git/tek2/rstats/read.R")
.indirs = fs::dir_ls(regexp = "^n\\d+", type = "directory")
.metadata = read_metadata(.indirs) %>% print()

 .repl = 43L; .lbound = 25000L; .ubound = 45000L

df2 = .metadata %>%
  dplyr::filter(n == 1000L, repl == .repl) %>%
  dplyr::mutate(
    adata = purrr::map(indir, ~{read_activity(.) %>% dplyr::filter(.lbound <= generation, generation <= .ubound)}),
    tdata = purrr::map(indir, ~{
      read_fastas(.x, interval = 100L, from = .lbound, to = .ubound) %>%
        add_phylo() %>%
        eval_treeshape()
    })
  ) %>%
  print()

df2$tdata[[1]] %>% dplyr::arrange(bimodality)
df2$tdata[[1]] %>% dplyr::arrange(abs(bimodality - 5 / 9))

.gens = list()
.gens[["43"]] = c(25500L, 29700L, 33900L, 38400L, 41200L, 42600L)
.timepoints_df = df2$tdata[[1]] %>%
  dplyr::filter(generation %in% .gens[[as.character(.repl)]]) %>%
  dplyr::transmute(generation, value = bimodality, stat = "bimodality") %>%
  print()

.fig2_left = df2 %>%
  dplyr::transmute(
    aplot = purrr::map2(adata, n, ~{
      ggplot_activity(.x, .y) +
        geom_vline(aes(xintercept = generation), .timepoints_df, linetype = "dashed", colour = "#666666") +
        theme(panel.grid.major.x = element_blank())
    }),
    tplot = purrr::map(tdata, ~{
      ggplot_evolution(.x, only_bi = TRUE) +
        geom_vline(aes(xintercept = generation), .timepoints_df, linetype = "dashed", colour = "#666666") +
        geom_point(data = .timepoints_df, colour = "red", size = 3) +
        theme(
          panel.grid.major.x = element_blank(),
          axis.title = element_blank(), axis.text.x = element_blank(),
          axis.ticks.x = element_blank()
        )
    })
  ) %>%
  purrr::pmap(function(aplot, tplot, ...) {
    .scale_x = scale_x_continuous(breaks = c(20000, 30000, 40000))
    aplot = aplot + theme(legend.position = "none") + .scale_x
    cowplot::plot_grid(tplot + .scale_x, aplot, ncol = 1, rel_heights = c(1, 1), align = "v", axis = "lr")
  }) %>%
  purrr::pluck(1L) +
  coord_fixed(ratio = 1.0)
.fig2_left
ggsave("fig2_left.pdf", .fig2_left, width = 4, height = 7, family = "Helvetica")

.infiles = fs::path(.metadata$indir[.repl], sprintf("generation_%05d.fa.gz", .gens[[as.character(.repl)]]))
.all_inds_df = .infiles %>%
  setNames(str_extract(., "(?<=generation_)\\d+")) %>%
  purrr::map_dfr(read_individuals, .id = "generation") %>%
  dplyr::mutate(generation = as.integer(generation)) %>%
  print()

.total_df = .all_inds_df %>% dplyr::filter(individual == "total") %>% print()
.inds_df = .all_inds_df %>% dplyr::filter(individual != "total") %>% print()

.df_unrooted = .inds_df %>%
# .df_unrooted_total = .total_df %>%
  purrr::pmap_df(function(phylo, data, generation, individual, ...) {
    wtl::ape_layout_unrooted(phylo) %>%
      dplyr::left_join(data, by = "label") %>%
      dplyr::mutate(generation = generation, individual = individual)
  }) %>%
  print()

.p = .df_unrooted %>%
# .p = .df_unrooted_total %>%
  ggplot_tetree() +
  facet_grid(generation ~ individual) +
  # facet_wrap(~ generation)+
  theme_classic()
.p
ggsave(sprintf("fig2_tree_candidates_%d.pdf", .repl), .p, width = 9.9, height = 7, scale = 2)
# ggsave("fig2_right_unrooted_total.pdf", .p, width = 9.9, height=7, scale=2)

.inds = as.character(c(2, 7, 7, 1, 8, 1))
.delegates = .timepoints_df %>%
  dplyr::transmute(generation, BI = sprintf("%.3f", value), individual = .inds) %>%
  print()
.delegates_df = .delegates %>%
  dplyr::left_join(.df_unrooted %>% tidyr::nest(-generation, -individual), by = c("generation", "individual")) %>%
  tidyr::unnest() %>%
  print()

.annot_base = .delegates_df %>% dplyr::distinct(generation, BI) %>% tail(1L) %>% print()
.dsegm = .annot_base %>% dplyr::mutate(x = 5, xend = 8, y = -7.5, yend = -7.5) %>% print()
.dtext = .annot_base %>% dplyr::mutate(x = 6.5, y = -8.7, label = "0.01") %>% print()

.fig2_right = ggplot_tetree(.delegates_df) +
  geom_segment(data = .dsegm, aes(x, y, xend = xend, yend = yend), size = 1) +
  geom_text(data = .dtext, aes(x, y, label = label)) +
  facet_wrap(~ generation + BI, labeller = label_both_tree, nrow = 2L) +
  theme_classic() +
  theme(
    axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(),
    strip.background = element_blank(),
    strip.text.x = element_text(size = 12, hjust = 0, margin = margin(6, 12, 3, 12))
  )
# .fig2_right
ggsave("fig2_right.pdf", .fig2_right, width = 9, height = 7, family = "Helvetica", device = cairo_pdf)

fig2 = cowplot::plot_grid(.fig2_left, .fig2_right, rel_widths = c(4, 9))
fig2
ggsave("fig2.pdf", fig2, width = 12, height = 7, family = "Helvetica", device = cairo_pdf)
