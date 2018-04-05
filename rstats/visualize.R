library(tidyverse)
library(wtl)
loadNamespace("cowplot")

extract_params = function(filename, params=c("alpha", "beta", "lambda", "xi", "nu", "lower", "upper")) {
  patterns = sprintf("_%s([^_]+)_", params)
  str_match(paste0("_", filename), patterns)[, 2] %>%
    parse_double() %>%
    set_names(params) %>%
    as.list() %>%
    as_tibble()
}
# .infiles[1] %>% extract_params()

popsize = 500
# popsize = 1000
.indirs = list.dirs(full.names = FALSE, recursive = FALSE)
# .indirs = wtl::command_args()$args
# .indirs = "."

.metadata = .indirs %>%
  str_subset("_\\d+$") %>%
  set_names() %>%
  purrr::map_dfr(extract_params, .id = "indir") %>%
  dplyr::group_by(xi, lower, upper) %>%
  dplyr::mutate(repl = seq_len(n())) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(xi, lower, upper, repl) %>%
  print()

.metadata %>% distinct(xi, lower, upper)

# .metadata$indir[[1]] %>% read_activity() %>% ggplot_activity()

source('~/git/tek-evolution/rstats/activity.R')
source('~/git/tek-evolution/rstats/treestats.R')

.out = .metadata %>%
  # dplyr::filter(!(xi < 6e-4 & upper < 16)) %>%
  # head(3L) %>%
  dplyr::mutate(
    title = sprintf('xi=%.0e lower=%d upper=%d (%d)', xi, lower, upper, repl),
    adata = purrr::map(indir, read_activity),
    aplot = purrr::map(adata, ggplot_activity, popsize=popsize),
    tplot = purrr::map2(indir, title, ~{
      read_ggplot_treestats(.x, interval=500L, title=.y)+theme(
        axis.title = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()
      )
    })
  ) %>%
  print()

.plt = .out %>% {purrr::map2(.$aplot, .$tplot, ~{
  cowplot::plot_grid(.y, .x, ncol=1, rel_heights=c(1, 1), align='v', axis='lr')
})}

# cowplot::plot_grid(plotlist=.plt)
.gtable = gridExtra::marrangeGrob(.plt, nrow=1, ncol=3, top=NULL)
ggsave('copynumber-treestats.pdf', .gtable, width=9.9, height=7)

.repl = 5L; .lbound = 25000L; .ubound = 31000L
.repl = 2L; .lbound = 10000L; .ubound = 21000L
.out = .metadata %>%
  dplyr::filter(repl == .repl) %>%
  dplyr::mutate(
    adata = purrr::map(indir, ~{read_activity(.) %>% dplyr::filter(.lbound <= generation, generation <= .ubound)}),
    aplot = purrr::map(adata, ggplot_activity),
    tdata = purrr::map(indir, ~{
      read_fastas(indir, interval=100L) %>%
        dplyr::filter(.lbound <= generation, generation <= .ubound) %>%
        add_phylo() %>%
        eval_treeshape()
    }),
    tplot = purrr::map(tdata, ~{
      ggplot_evolution(.x)+
        geom_point(data=.points_data, colour='red', size=3)+
        theme(
          axis.title = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank()
        )
    }),
    plt = purrr::map2(aplot, tplot, ~{
      cowplot::plot_grid(.y, .x, ncol=1, rel_heights=c(1, 1), align='v', axis='lr')
    })
  ) %>%
  print()
.fig2_left = .out$plt[[1]]
.fig2_left
ggsave('fig2_left.pdf', .fig2_left, width=3, height=9)

.out$tdata[[1]] %>% dplyr::filter(generation < 12000L) %>% dplyr::arrange(bimodality)
.out$tdata[[1]] %>% dplyr::filter(generation > 18000L) %>% dplyr::arrange(bimodality)
.out$tdata[[1]] %>% dplyr::arrange(abs(bimodality - 5/9))

.gens = c(10600L, 15200L, 19100L, 20300L)
.gens = c(10600L, 13800L, 18700L, 20500L)
.points_data = .out$tdata[[1]] %>%
  dplyr::filter(generation %in% .gens) %>%
  dplyr::transmute(generation, value = bimodality, stat = "bimodality")


# #######1#########2#########3#########4#########5#########6#########7#########

source('~/git/tek-evolution/rstats/treeplot.R')

.focus = .metadata %>% dplyr::filter(xi == max(xi), lower == 9, upper == 24) %>% print()
.focus = .metadata %>% dplyr::filter(xi == max(xi), lower == 300L, upper == 300L) %>% print()

.infiles = fs::path(.focus$indir[2], sprintf('generation_%05d.fa.gz', seq(10000L, 20000L, by=2000L)))
.infiles = fs::path(.focus$indir[5], sprintf('generation_%05d.fa.gz', seq(25000L, 30000L, by=1000L)))

.infiles = fs::path(.focus$indir[2], sprintf('generation_%05d.fa.gz', .gens))

read_individuals = function(infile) {
  message(infile)
  .seqs = read_tek_fasta(infile, metadata=TRUE)
  .mcols_all = tidy_metadata(.seqs)
  .nested_mcols = nest_metadata(.mcols_all, .seqs)
}

.all_inds_df = .infiles %>%
  setNames(str_extract(.,'(?<=generation_)\\d+')) %>%
  purrr::map_dfr(read_individuals, .id='generation') %>%
  dplyr::mutate(generation = as.integer(generation)) %>%
  print()

fortify_phylo_tbl = function(.tbl, layout = "rectangular") {
  purrr::pmap_df(.tbl, function(phylo, data, generation, individual, ...) {
    ggtree::fortify(phylo, layout = layout) %>%
      dplyr::left_join(data, by="label") %>%
      dplyr::mutate(generation = generation, individual = individual)
  })
}

.inds_df = .all_inds_df %>%
  dplyr::filter(individual != 'total') %>%
  print()

.max_copy_number = .inds_df$data %>% purrr::map_int(~max(.x$copy_number)) %>% max() %>% print()

.df_rect = .inds_df %>% fortify_phylo_tbl() %>% print()
# .df_eqangle = .inds_df %>% fortify_phylo_tbl(layout="equal_angle") %>% print()
# .df_daylight = .inds_df %>% fortify_phylo_tbl(layout="daylight") %>% print()

.p = ggplot(.df_eqangle, aes(x, y)) +
  # geom_tree(layout="rectangular") +
  geom_tree(layout="equal_angle") +
  # geom_tree(layout="daylight") +
  geom_tippoint(aes(colour=activity, size=copy_number), alpha=0.6) +
  scale_colour_gradientn(colours = rev(head(rainbow(15L), 12L)), limits = c(0, 1), breaks = c(0, 0.5, 1)) +
  scale_size(limit=c(1, .max_copy_number), range=c(3, 12)) +
  geom_tippoint(aes(subset=is_major, size=copy_number), pch=1, colour='#000000')+
  geom_tippoint(aes(subset=is_fixed, size=copy_number), pch=4, colour='#000000')+
  facet_grid(generation ~ individual) +
  # facet_wrap(~ generation) +
  ggtree::theme_tree()
.p
ggsave("fig2_right_candidates.pdf", .p, width = 9.9, height=7, scale=2)
ggsave("fig2_right_candidates_eqangle.pdf", .p, width = 9.9, height=7, scale=2)

ggsave("fig2_right_total.pdf", .p, width = 2, height=7, scale=2)
ggsave("fig2_right_total_eqangle.pdf", .p, width = 2, height=7, scale=2)
ggsave("fig2_right_total_daylight.pdf", .p, width = 2, height=7, scale=2)

ggsave("tek-5-rectangular.pdf", .p, width = 9.9, height=7, scale=2)
ggsave("tek-5-equal_angle.pdf", .p, width = 9.9, height=7, scale=2)

# #######1#########2#########3#########4#########5#########6#########7#########

Rprof()
plot_individuals(.fagz, layout='unrooted')
Rprof(NULL)
summaryRprof()

fs::dir_ls(.focus$indir[2], regex='\\d+\\.fa\\.gz$')
.infiles = fs::path(.focus$indir[2], sprintf('generation_%05d.fa.gz', seq(0, 40000, by=4000)[-1]))
# .plts = purrr::map(.infiles, plot_individuals)
.plts = wtl::map_par(.infiles, plot_individuals)
.plts[[3]]
ggsave('individual_trees.pdf', .plts, width=10, height=10)

.fagz = fs::path(.focus$indir[2], "generation_16000.fa.gz")
.seqs = .fagz %>% read_tek_fasta(metadata=TRUE)
# .mcols = tidy_mcols(.seqs) %>% print()
# count_holders(.mcols)
.mcols_all = tidy_metadata(.seqs) %>% print()
.nested_mcols = .mcols_all %>% nest_metadata(.seqs) %>% print()
.max_copy_number = dplyr::filter(.nested_mcols, individual == 'total')$data[[1]]$copy_number %>% max()

.plts = .nested_mcols %>%
  # head(1) %>%
  purrr::pmap(ggtree_tek, max_copy_number = .max_copy_number)
.cow = cowplot::plot_grid(plotlist=.plts)
.cow
ggsave('individual_trees_rect.png', .cow, width=10, height=10)


.nested_mcols$data[[2]]

.phylo = .nested_mcols$phylo[[1]]
.p = ggtree_hide_root(.phylo) + geom_tiplab()
.p
ggplot(.phylo)$data

# #######1#########2#########3#########4
# draw sampled nodes on master tree

.useqs = unique_dss(.seqs)
.phylo = .useqs %>% Biostrings::stringDist(method='hamming') %>% ape::fastme.ols()
.rooted = ape::root(.phylo, '0x0', resolve.root=TRUE)
.gt = ggtree(.rooted, layout='rectangular', alpha=0.5)

make_facettable(.gt, .mcols_all, .max_copy_number)+
facet_wrap(~individual)
