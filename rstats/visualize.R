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
    aplot = purrr::map(adata, ggplot_activity),
    tplot = purrr::map2(indir, title, ~{
      read_ggplot_treestats(.x, interval=500L, title=.y)+theme(
        axis.title = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()
      )
    })
  ) %>%
  print()

.out2 = .out %>%
  dplyr::mutate(plt = purrr::map2(aplot, tplot, ~{
    cowplot::plot_grid(.y, .x, ncol=1, rel_heights=c(1, 1), align='v', axis='lr')
  }))

# cowplot::plot_grid(plotlist=.out2$plt)
.gtable = gridExtra::marrangeGrob(.out2$plt, nrow=1, ncol=3, top=NULL)
ggsave('copynumber-treestats.pdf', .gtable, width=9.9, height=7)

# #######1#########2#########3#########4#########5#########6#########7#########

source('~/git/tek-evolution/rstats/treeplot.R')

.focus = .metadata %>% dplyr::filter(xi == max(xi), lower == 9, upper == 24) %>% print()
.focus = .metadata %>% dplyr::filter(xi == max(xi), lower == 300L, upper == 300L) %>% print()

plot_individuals(.fagz)

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
