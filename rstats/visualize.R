library(tidyverse)
library(wtl)
loadNamespace("cowplot")

tek_params = c("n", "alpha", "beta", "lambda", "xi", "nu", "lower", "upper", "coexist")
extract_params = function(filename, params=tek_params) {
  patterns = sprintf("_%s([^_]+)_", params)
  str_match(paste0("_", filename), patterns)[, 2] %>%
    parse_double() %>%
    set_names(params) %>%
    as.list() %>%
    as_tibble()
}
# .infiles[1] %>% extract_params()

# .indirs = wtl::command_args()$args
# .indirs = "."
.indirs = fs::dir_ls(regexp = "^n\\d+", type = "directory")
.indirs = fs::dir_ls(regexp = "^r\\d+", type = "directory")

.metadata = .indirs %>%
  str_subset("_\\d+$") %>%
  set_names() %>%
  purrr::map_dfr(extract_params, .id = "indir") %>%
  dplyr::group_by(!!!rlang::syms(tek_params)) %>%
  dplyr::mutate(repl = seq_len(n())) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(n, xi, coexist, lower, upper, repl) %>%
  print()

.metadata %>% dplyr::count(!!!rlang::syms(tek_params))

# .metadata$indir[[1]] %>% read_activity() %>% ggplot_activity()

source("~/git/teaposon/rstats/activity.R")
source("~/git/teaposon/rstats/biostrings.R")
source("~/git/teaposon/rstats/treestats.R")

# #######1#########2#########3#########4#########5#########6#########7#########

fig1acti = .metadata %>%
  # dplyr::filter(xi > 6e-4) %>%
  dplyr::mutate(
    title = sprintf("n=%d xi=%.0e (%d)", n, xi, repl),
    adata = purrr::map(indir, read_activity),
    aplot = purrr::map2(adata, n, ggplot_activity),
    aplot = purrr::map2(aplot, title, ~.x + labs(title = .y) + theme(legend.position = "none"))
  ) %>%
  print()
# fig1acti$aplot[[1]]
ggsave("fig1_activity.pdf", fig1acti$aplot, width = 8, height = 4)

# @ILP23% ~/git/teaposon/run.py -p2 -j6 te2fig2 -r96

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

# @ilp15% ~/git/teaposon/run.py -p2 -j6 te2fig4 -r96

fig4cand = .metadata %>% dplyr::transmute(n, repl, lapply(indir, read_activity)) %>% print()
.p = fig4cand %>% tidyr::unnest() %>%
  ggplot_activity(1000, "#aa0000") +
  coord_cartesian(xlim = c(3000, 6000)) +
  facet_wrap(~n + repl)
.p
ggsave("fig4_candidates.png", .p, width = 16, height = 16)


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


.metadata %>% dplyr::count(!!!rlang::syms(tek_params))

# @ILP22% ~/git/teaposon/run.py -p2 -j6 -r6 te2fig6
# rsync -auv ilp22.local:~/working/te2-fig6/ ~/working/tek/te2-fig6/

df6act = .metadata %>%
  dplyr::mutate(
    label = sprintf("n=%d c=%d l=%d u=%d repl=%d", n, coexist, lower, upper, repl),
    adata = purrr::map(indir, read_activity),
    indir = NULL
  ) %>%
  tidyr::unnest() %>%
  dplyr::mutate(copy_number = copy_number / n) %>%
  print()

fig6candidates = df6act %>%
  dplyr::filter((generation %% 2000) == 0) %>%
  ggplot_activity() +
  facet_wrap(~label, ncol = 6L) +
  theme_minimal()
fig6candidates
ggsave("fig6candidates.pdf", fig6candidates, width = 7, height = 9.9, scale = 2)

fig6candidates_c8 = df6act %>%
  dplyr::filter(coexist == 8L) %>%
  dplyr::filter((generation %% 500) == 0) %>%
  ggplot_activity() +
  facet_wrap(~label, ncol = 6L) +
  theme_minimal()
fig6candidates_c8
ggsave("fig6candidates-c8.pdf", fig6candidates_c8, width = 9.9, height = 7, scale = 2)

.fig6repl = tibble::tribble(
     ~n, ~coexist, ~lower, ~upper, ~repl,
  # 1000L,       8L,     6L,    18L,    1L,
  1000L,       8L,     6L,    18L,    3L,
  1000L,       8L,     6L,    30L,    1L,
  1000L,       8L,     9L,    18L,    5L,
  1000L,       8L,     9L,    30L,    1L,
) %>% print()

df6 = .metadata %>%
  dplyr::right_join(.fig6repl, by = names(.fig6repl)) %>%
  dplyr::mutate(
    label = sprintf("n=%d l=%d u=%d repl=%d", n, lower, upper, repl),
    adata = purrr::map(indir, read_activity),
    indir = NULL
  ) %>%
  tidyr::unnest() %>%
  dplyr::mutate(copy_number = copy_number / n) %>%
  print()

fig6 = df6 %>%
  dplyr::rename(lowerbound = lower, upperbound = upper) %>%
  ggplot_activity() +
  facet_grid(upperbound ~ lowerbound, labeller = label_both) +
  theme(strip.background = element_blank())
fig6
ggsave("fig6.pdf", fig6, width = 10, height = 4)

# ~/git/teaposon/run.py -p3 -j2 -r10 te2fig6

# #######1#########2#########3#########4#########5#########6#########7#########

fig1df = .metadata %>%
  dplyr::filter(n == 1000L) %>%
  dplyr::mutate(
    adata = wtl::mcmap(indir, read_activity),
    tdata = wtl::mcmap(indir, ~{
      message(.x)
      read_fastas(.x, interval = 100L) %>%
        add_phylo() %>%
        eval_treeshape()
    })
  ) %>%
  print()

saveRDS(fig1df, "fig1.rds")

# rsync -auv ilp15.local:~/working/te2-fig1/ ~/working/tek/te2-fig1/
# fig1df = readRDS("fig1.rds")

fig1plts = fig1df %>%
  dplyr::arrange(desc(xi)) %>%
  dplyr::mutate(
    tplot = purrr::map(tdata, ~{
      .xmin = min(.x$generation)
      ggplot_evolution(.x %>% dplyr::filter(generation > 500L), only_bi = TRUE) +
        coord_cartesian(xlim = c(.xmin, max(.x$generation))) +
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
    plt = purrr::map2(tplot, aplot, ~{
      cowplot::plot_grid(.x, .y + theme(legend.position = "none", axis.title.x = element_blank()),
        ncol = 1, rel_heights = c(1, 1), align = "v", axis = "lr")
    })
  ) %>%
  print()
fig1plts$plt[[1]]

fig1 = cowplot::plot_grid(
  cowplot::plot_grid(
    cowplot::plot_grid(plotlist = fig1plts$plt, ncol = 1, labels = c("A", "B", "C"), scale = 0.95),
    cowplot::get_legend(fig1plts$aplot[[1]]),
    nrow = 1L, rel_widths = c(1, 0.1)
  ),
  cowplot::ggdraw() + cowplot::draw_label("Generation"),
  ncol = 1L, rel_heights = c(3, 0.1)
)
ggsave("fig1.pdf", fig1, width = 10, height = 10)

# #######1#########2#########3#########4#########5#########6#########7#########

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
# .gens = c(10600L, 15200L, 19100L, 20300L)
# .gens[["2"]] = c(10600L, 12300L, 14600L, 16800L, 18700L, 20500L)
# .gens[["48"]] = c(25000L, 26900L, 30100L, 33400L, 34100L, 34700L)
# .gens[["69"]] = c(42600L, 43900L, 45200L, 46300L, 47100L, 47800L)
# .gens[["77"]] = c(44900L, 45700L, 46500L, 47500L, 48700L, 49600L)
.gens[["43"]] = c(25500L, 29700L, 33900L, 38400L, 41200L, 42600L)
.timepoints_df = df2$tdata[[1]] %>%
  dplyr::filter(generation %in% .gens[[as.character(.repl)]]) %>%
  dplyr::transmute(generation, value = bimodality, stat = "bimodality") %>%
  print()

.fig2_left = df2 %>%
  dplyr::transmute(
    aplot = purrr::map2(adata, n, ggplot_activity),
    tplot = purrr::map(tdata, ~{
      ggplot_evolution(.x, only_bi = TRUE) +
        geom_point(data = .timepoints_df, colour = "red", size = 3) +
        theme(
          axis.title = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank()
        )
    })
  ) %>%
  purrr::pmap(function(aplot, tplot, ...) {
    .scale_x = scale_x_continuous(breaks = c(20000, 30000, 40000))
    aplot = aplot + theme(legend.position = "none") + .scale_x
    cowplot::plot_grid(tplot + .scale_x, aplot, ncol = 1, rel_heights = c(1, 1), align = "v", axis = "lr")
  }) %>%
  purrr::pluck(1L)
ggsave("fig2_left.pdf", .fig2_left, width = 4, height = 9, family = "Helvetica")


# #######1#########2#########3#########4

read_individuals = function(infile) {
  message(infile)
  .seqs = read_tek_fasta(infile, metadata = TRUE)
  .mcols_all = tidy_metadata(.seqs, add_root = FALSE)
  .nested_mcols = nest_metadata(.mcols_all, .seqs)
}

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

fig2 = cowplot::plot_grid(.fig2_left, .fig2_right, rel_widths = c(3.2, 9))
ggsave("fig2.pdf", fig2, width = 12, height = 7, family = "Helvetica", device = cairo_pdf)

# #######1#########2#########3#########4

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

# #######1#########2#########3#########4

fig5df$adata[[1]] %>%
  dplyr::count(species, generation, wt = copy_number) %>%
  dplyr::group_by(species) %>%
  dplyr::group_modify(~head(.))


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


# #######1#########2#########3#########4#########5#########6#########7#########
# Check continuity of each sequence

.delegates_df %>%
# .df_unrooted %>%
  dplyr::filter(!is.na(label)) %>%
  dplyr::count(generation, label, activity, wt = copy_number) %>%
  ggplot(aes(generation, n, group = label)) +
  geom_line(aes(colour = label), alpha = 0.5, size = 3) +
  # scale_colour_gradientn(colours = rev(head(rainbow(15L), 12L)), limits = c(0, 1), breaks = c(0, 0.5, 1)) +
  theme_bw() + theme(legend.position = "none")
