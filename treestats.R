library(tidyverse)
library(Biostrings)
library(ggtree)
library(ape)
library(aptree)  # library(apTreeshape)
wtl::refresh('aptree')

read_fastas = function(dir, interval = 1000L) {
  .fastas = fs::dir_ls(dir, regexp='generation_\\d+\\.fa\\.gz$')
  tibble::tibble(
    path = .fastas,
    infile = fs::path_file(path),
    generation = as.integer(readr::parse_number(infile))) %>%
  dplyr::filter((generation %% interval) == 0L) %>%
  dplyr::transmute(
    generation,
    seqs = purrr::map(path, ~{
      Biostrings::readBStringSet(.x) %>%
      setNames(str_extract(names(.), "te=\\S+"))
    })
  )
}
# .tbl = read_fastas('lower10_upper30_20180130T172206_00') %>% print()

add_phylo = function(.tbl) {
  .tbl %>%
    dplyr::mutate(distmat = purrr::map(seqs, Biostrings::stringDist, method='hamming')) %>%
    dplyr::mutate(phylo = purrr::map(distmat, ~{
      if (length(.x) > 2) {ape::fastme.ols(.x)} else {NA}
    }))
}
# .tblphy = .tbl %>% add_phylo() %>% print()
# .phy = .tblphy$phylo[[5L]]
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
  facet_grid(stat ~ ., scale='free_y')+
  theme_bw()
}
# .tblstats %>% ggplot_evolution()

main = function(.x, interval=2000L) {
  message(.x)
  label = fs::path_file(.x)
  outfile = paste0("shapestats-", label, ".png")
  .p = read_fastas(.x, interval) %>%
    add_phylo() %>%
    eval_treeshape() %>%
    ggplot_evolution()
  ggsave(outfile, .p, width=7, height=7)
}

root = '~/working/te2-0207'
indirs = fs::dir_ls(root, regexp="\\d+$", type="directory") %>%
    str_subset('xi10e|xi5e') %>%
    print()
# indirs[1] %>% main()
purrr::walk(indirs, main, interval = 500L)

# #######1#########2#########3#########4#########5#########6#########7#########

eval_treeshape = function(.tblphy) {
  .tblphy %>%
    dplyr::mutate(sstat = purrr::map(phylo, ~{
      aptree::as.treeshape(.x, 'pda') %>%
      aptree::calc_stat_all()
    })) %>%
    dplyr::select(generation, sstat) %>%
    tidyr::unnest()
}
# .tblstats = .tblphy %>% eval_treeshape() %>% print()

ggplot_evolution = function(.tblstats) {
  ggplot(.tblstats, aes(generation, stat))+
  geom_line()+
  facet_grid(norm ~ index, scale='free_y')+
  theme_bw()
}
# .tblstats %>% ggplot_evolution()

library(doParallel)

root = '~/working/tek-spec-0130'
indirs = fs::dir_ls(root, regexp="\\d+$", type="directory") %>%
    str_subset('xi10e|xi5e') %>%
    print()
wtl::map_par(indirs, main, interval = 500L)
