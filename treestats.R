library(tidyverse)
library(Biostrings)
library(ape)
library(apTreeshape)

.fastas = fs::dir_ls(glob='generation_*.fa.gz')
# Biostrings::readBStringSet(.fastas[[3L]])

.tbl = tibble::tibble(
    infile = .fastas,
    seqset = purrr::map(infile, ~{
      Biostrings::readBStringSet(.x) %>%
      setNames(str_extract(names(.), "te=\\S+"))
    })
  ) %>% print()
# .tbl$seqset[[3L]] %>% Biostrings::stringDist(method='hamming')

.tblphy = .tbl %>%
  # tail(-1L) %>%  # to few nodes
  dplyr::mutate(distmat = purrr::map(seqset, Biostrings::stringDist, method='hamming')) %>%
  print() %>%
  dplyr::mutate(phy = purrr::map(distmat, ape::nj)) %>%
  print()

# .phy = .tblphy$phy[[3L]]
# plot(.phy, type = "u")
# .ts = apTreeshape::as.treeshape(.phy, model = 'pda') %>% print()
# apTreeshape::shape.statistic(.ts, norm = 'pda')

crossing_model_norm = function(.phy) {
  .shapes = tibble::tibble(
    model = c('biased', 'pda', 'yule'),
    tree_shape = purrr::map(model, ~{apTreeshape::as.treeshape(.phy, .x)})
  )
  tidyr::crossing(model = .shapes$model, norm = c('null', 'pda', 'yule')) %>%
    dplyr::left_join(.shapes, by = 'model') %>%
    dplyr::mutate(shape_stat = purrr::map2_dbl(tree_shape, norm, ~{
      if (.y == 'null') .y = NULL
      apTreeshape::shape.statistic(.x, norm=.y)
    }))
}
# crossing_model_norm(.tblphy$phy[[1]])

.tblstats = .tblphy %>%
  dplyr::mutate(sstat = purrr::map(phy, crossing_model_norm)) %>%
  dplyr::transmute(generation = as.integer(readr::parse_number(infile)), sstat) %>%
  tidyr::unnest() %>%
  print()

.tblstats %>%
  dplyr::select(-tree_shape) %>%
  ggplot(aes(generation, shape_stat, colour=model))+
  geom_line()+
  facet_grid(norm ~ ., scale='free_y')+
  theme_bw()
