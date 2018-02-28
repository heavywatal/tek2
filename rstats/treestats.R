library(tidyverse)
library(Biostrings)
library(ape)
library(aptree)  # library(apTreeshape)
# wtl::refresh('aptree')

parse_fasta_header = function(x) {
  str_match_all(x, '(\\w+)=(\\S+)') %>%
  purrr::map_dfr(~{
    tibble::tibble(key = .x[,2], value = .x[,3]) %>%
    tidyr::spread(key, value, convert=FALSE)
  }) %>%
  dplyr::mutate_at(vars(activity, dn, ds), as.double) %>%
  dplyr::mutate_at(vars(copy_number, indel, species), as.integer)
}

read_tek_fasta = function(file, metadata=FALSE, nrec=-1L, skip=0L) {
  .dss = Biostrings::readDNAStringSet(file, nrec=nrec, skip=skip)
  .names = names(.dss)
  if (metadata) {
    mcols(.dss) = parse_fasta_header(.names)
    names(.dss) = mcols(.dss)$te
  } else {
    names(.dss) = str_extract(.names, "(?<=te=)\\S+")
  }
  .dss
}

read_fastas = function(dir, interval = 1000L) {
  .fastas = fs::dir_ls(dir, regexp='generation_\\d+\\.fa\\.gz$')
  tibble::tibble(
    path = .fastas,
    infile = fs::path_file(path),
    generation = as.integer(readr::parse_number(infile))) %>%
  dplyr::filter((generation %% interval) == 0L) %>%
  dplyr::transmute(
    generation,
    seqs = purrr::map(path, read_tek_fasta)
  )
}
# .tbl = read_fastas('lower10_upper30_20180130T172206_00') %>% print()

add_phylo = function(.tbl, root=NULL) {
  .tbl %>%
    dplyr::mutate(distmat = purrr::map(seqs, Biostrings::stringDist, method='hamming')) %>%
    dplyr::mutate(phylo = purrr::map(distmat, ~{
      if (length(.x) > 2) {
        .phy = ape::fastme.ols(.x)
        if (!is.null(root)) {
          .phy = ape::root(.phy, root, resolve.root=TRUE)
        }
        .phy
      } else {NA}
    }))
}
# .tblphy = .tbl %>% add_phylo() %>% print()
# .phy = .tblphy$phylo[[5L]]
# library(ggtree)
# ggtree(.phy)

eval_treeshape = function(.tblphy) {
  dplyr::transmute(.tblphy,
    generation,
    shape_pda = purrr::map_dbl(phylo, ~{
      if (identical(.x, NA)) return(NA)
      aptree::as.treeshape(.x, 'pda') %>%
      aptree::shape.statistic(norm = 'pda')
    }),
    bimodality = purrr::map_dbl(distmat, wtl::bimodality)
  )
}
# .tblstats = .tblphy %>% eval_treeshape() %>% print()

range_stats = tibble::tibble(
  generation = 0L,
  stat = rep(c('bimodality', 'shape_pda'), each=2L),
  value = c(0, 1, 2, -10)) %>%
  print()

ggplot_evolution = function(.tblstats) {
  .x = tidyr::gather(.tblstats, stat, value, -generation)
  ggplot(.x, aes(generation, value))+
  geom_blank(data=range_stats)+
  geom_line()+
  facet_grid(stat ~ ., scale='free_y', switch='y')+
  theme_bw()+
  theme(
    panel.grid.minor = element_blank(),
    strip.placement = "outside",
    strip.background = element_blank()
  )
}
# .tblstats %>% ggplot_evolution()

read_ggplot_treestats = function(indir, interval=2000L, title='') {
  read_fastas(indir, interval=interval) %>%
    add_phylo() %>%
    eval_treeshape() %>%
    ggplot_evolution()+
    labs(title=title)
}


# #######1#########2#########3#########4#########5#########6#########7#########
if (FALSE) {

main = function(indir, interval=2000L) {
  message(indir)
  label = fs::path_file(indir)
  read_ggplot_treestats(indir, interval=interval, title=label)
}

root = '~/working/te2-0207'
indirs = fs::dir_ls(root, regexp="\\d+$", type="directory") %>%
    str_subset('xi10e|xi5e') %>%
    print()
# indirs[1] %>% main()
.plt = purrr::map(indirs, main, interval = 1000L)
ggsave('shapestats-0207.pdf', .plt, height=5, width=5)

eval_treeshape_all = function(.tblphy) {
  .tblphy %>%
    dplyr::mutate(sstat = purrr::map(phylo, ~{
      aptree::as.treeshape(.x, 'pda') %>%
      aptree::calc_stat_all()
    })) %>%
    dplyr::select(generation, sstat) %>%
    tidyr::unnest()
}
# .tblstats = .tblphy %>% eval_treeshape_all() %>% print()
# .tblstats %>% ggplot_evolution()

library(doParallel)

root = '~/working/tek-spec-0130'
indirs = fs::dir_ls(root, regexp="\\d+$", type="directory") %>%
    str_subset('xi10e|xi5e') %>%
    print()
wtl::map_par(indirs, main, interval = 500L)

} # if (FALSE)
