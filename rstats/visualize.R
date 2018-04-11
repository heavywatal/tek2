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
source('~/git/tek-evolution/rstats/biostrings.R')
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

# .repl = 5L; .lbound = 25000L; .ubound = 31000L
.repl = 2L; .lbound = 10000L; .ubound = 21000L

.fig2df = .metadata %>%
  dplyr::filter(repl == .repl) %>%
  dplyr::mutate(
    adata = purrr::map(indir, ~{read_activity(.) %>% dplyr::filter(.lbound <= generation, generation <= .ubound)}),
    tdata = purrr::map(indir, ~{
      read_fastas(indir, interval=100L) %>%
        dplyr::filter(.lbound <= generation, generation <= .ubound) %>%
        add_phylo() %>%
        eval_treeshape()
    })
  ) %>%
  print()

.fig2df$tdata[[1]] %>% dplyr::filter(generation < 12000L) %>% dplyr::arrange(bimodality)
.fig2df$tdata[[1]] %>% dplyr::filter(generation > 18000L) %>% dplyr::arrange(bimodality)
.fig2df$tdata[[1]] %>% dplyr::arrange(abs(bimodality - 5/9))

# .gens = c(10600L, 15200L, 19100L, 20300L)
.gens = c(10600L, 12300L, 14600L, 16800L, 18700L, 20500L)
.timepoints_df = .fig2df$tdata[[1]] %>%
  dplyr::filter(generation %in% .gens) %>%
  dplyr::transmute(generation, value = bimodality, stat = "bimodality") %>%
  print()

.out = .fig2df %>% dplyr::transmute(
  aplot = purrr::map(adata, ggplot_activity, popsize=popsize),
  tplot = purrr::map(tdata, ~{
    ggplot_evolution(.x, only_bi = TRUE)+
      geom_point(data=.timepoints_df, colour='red', size=3)+
      theme(
        axis.title = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()
      )
  }),
  plt = purrr::map2(aplot, tplot, ~{
    .legend = cowplot::get_legend(.x + guides(fill = guide_colorbar(barheight = 12)))
    .x = .x + theme(legend.position = "none")
    cowplot::plot_grid(
      cowplot::plot_grid(.y, .x, ncol=1, rel_heights=c(1, 1), align='v', axis='lr'),
      .legend,
      rel_widths = c(1, 0.3)
    )
  })
)
.fig2_left = .out$plt[[1]]
.fig2_left
ggsave('fig2_left.pdf', .fig2_left, width=4, height=9, family="Helvetica")


# #######1#########2#########3#########4#########5#########6#########7#########

.focus = .metadata %>% dplyr::filter(xi == max(xi), lower == 9, upper == 24) %>% print()
.focus = .metadata %>% dplyr::filter(xi == max(xi), lower == 300L, upper == 300L) %>% print()

.infiles = fs::path(.focus$indir[2], sprintf('generation_%05d.fa.gz', seq(10000L, 20000L, by=2000L)))
.infiles = fs::path(.focus$indir[5], sprintf('generation_%05d.fa.gz', seq(25000L, 30000L, by=1000L)))

.infiles = fs::path(.metadata$indir[2], sprintf('generation_%05d.fa.gz', .gens))

read_individuals = function(infile) {
  message(infile)
  .seqs = read_tek_fasta(infile, metadata=TRUE)
  .mcols_all = tidy_metadata(.seqs, add_root=FALSE)
  .nested_mcols = nest_metadata(.mcols_all, .seqs)
}

.all_inds_df = .infiles %>%
  setNames(str_extract(.,'(?<=generation_)\\d+')) %>%
  purrr::map_dfr(read_individuals, .id='generation') %>%
  dplyr::mutate(generation = as.integer(generation)) %>%
  print()

.phylo = .all_inds_df$phylo[[33L]]
.phylo = .all_inds_df$phylo[[23L]]
.total_df = .all_inds_df %>% dplyr::filter(individual == 'total') %>% print()
.inds_df = .all_inds_df %>% dplyr::filter(individual != 'total') %>% print()

.max_copy_number = .inds_df$data %>% purrr::map_int(~max(.x$copy_number)) %>% max() %>% print()

.df_unrooted = .inds_df %>%
# .df_unrooted_total = .total_df %>%
  purrr::pmap_df(function(phylo, data, generation, individual, ...) {
    wtl::ape_layout_unrooted(phylo) %>%
      dplyr::left_join(data, by="label") %>%
      dplyr::mutate(generation = generation, individual = individual)
  }) %>%
  print()

ggplot_tetree = function(data) {
  .ylim = data$yend %>% {c(min(.), max(.) * 1.1)}
  filter_active = function(x) dplyr::filter(x, activity > 0)
  filter_major = function(x) dplyr::filter(x, is_major)
  ggplot(data)+
  geom_point(data=filter_active, aes(xend, yend, colour=activity, size=copy_number), pch=16, alpha=0.6) +
  geom_point(data=filter_major, aes(xend, yend, size=copy_number), pch=1, colour='#000000', alpha=0.3, stroke=1) +
  geom_point(data=filter_major, aes(xend, yend, colour=activity, size=copy_number), pch=1, alpha=0.6, stroke=1) +
  geom_segment(aes(x, y, xend=xend, yend=yend), size=0.25) +
  scale_colour_gradientn(colours = rev(head(rainbow(15L), 12L)), limits = c(0, 1), breaks = c(0, 0.5, 1), guide=FALSE) +
  scale_size(limit=c(1, .max_copy_number), range=c(3, 12), name = "Copy\nNumber", breaks=c(9, 6, 3, 1)) +
  labs(x=NULL, y=NULL) +
  coord_fixed(ylim = .ylim)
}

.p = .df_unrooted %>%
# .p = .df_unrooted_total %>%
  ggplot_tetree()+
  facet_grid(generation ~ individual)+
  # facet_wrap(~ generation)+
  theme_classic()
.p
ggsave("fig2_right_unrooted_candidates.pdf", .p, width = 9.9, height=7, scale=2)
# ggsave("fig2_right_unrooted_total.pdf", .p, width = 9.9, height=7, scale=2)

.gens
.inds = c("1", "4", "2", "7", "0", "6")
.delegates = .timepoints_df %>%
  dplyr::transmute(generation, BI = sprintf("%.3f", value), individual = .inds) %>%
  print()
.delegates_df = .delegates %>%
  dplyr::left_join(.df_unrooted %>% tidyr::nest(-generation, -individual), by = c("generation", "individual")) %>%
  tidyr::unnest() %>%
  print()

.fig2_right = ggplot_tetree(.delegates_df)+
  facet_wrap(~ generation + BI, labeller=label_both)+
  theme_classic()+
  theme(
    axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(),
    strip.background = element_blank(),
    strip.text.x = element_text(size=12, hjust=0, margin = margin(6, 12, 3, 12))
  )
.fig2_right
ggsave('fig2_right.pdf', .fig2_right, width=7, height=4, family="Helvetica")

.fig2 = cowplot::plot_grid(.fig2_left, .fig2_right, rel_widths=c(2, 6))
ggsave('fig2.pdf', .fig2, width=14, height=7, family="Helvetica")


# #######1#########2#########3#########4#########5#########6#########7#########
# Check continuity of each sequence

.delegates_df %>%
# .df_unrooted %>%
  dplyr::filter(!is.na(label)) %>%
  dplyr::count(generation, label, activity, wt=copy_number) %>%
  ggplot(aes(generation, n, group=label))+
  geom_line(aes(colour=label), alpha=0.5, size=3)+
  # scale_colour_gradientn(colours = rev(head(rainbow(15L), 12L)), limits = c(0, 1), breaks = c(0, 0.5, 1)) +
  theme_bw()+theme(legend.position='none')
