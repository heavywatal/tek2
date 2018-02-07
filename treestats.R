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
    dplyr::mutate(phylo = purrr::map(distmat, ape::fastme.ols))
}
# .tblphy = .tbl %>% add_phylo() %>% print()
# .phy = .tblphy$phylo[[5L]]
# ggtree(.phy)

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

# library(doParallel)
main = function(root, interval = 2000L) {
  indirs = fs::dir_ls(root, regexp="\\d+$", type="directory")
  # purrr::walk(indirs, ~{
  # wtl::map_par(indirs, ~{
    message(.x)
    label = fs::path_file(.x)
    outfile = paste0("shapestats-", label, ".png")
    .p = read_fastas(.x, interval) %>%
      add_phylo() %>%
      eval_treeshape() %>%
      ggplot_evolution()
    ggsave(outfile, .p, width=7, height=7)
  })
}

root = '~/working/tek-spec-0130'
main(root, 500L)
