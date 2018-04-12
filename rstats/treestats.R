library(tidyverse)
library(ape)
library(aptree)  # library(apTreeshape)
# wtl::refresh('aptree')

add_phylo = function(.tbl, root=NULL) {
  .tbl %>%
    dplyr::mutate(distmat = purrr::map(seqs, Biostrings::stringDist, method='hamming')) %>%
    dplyr::mutate(phylo = purrr::map(distmat, ~{
      if (length(.x) > 2) {
        .phy = ape::fastme.ols(.x)
        if (!is.null(root) && (root %in% .phy$tip.label)) {
          .phy = ape::root(.phy, root, resolve.root=TRUE)
        }
        .phy
      } else {NA}
    }))
}

nest_metadata = function(.tidy_metadata, dss) {
  .useqs = unique_dss(dss)
  .tidy_metadata %>%
    tidyr::nest(-individual) %>%
    dplyr::mutate(seqs = purrr::map(data, ~{.useqs[.x$label]})) %>%
    add_phylo(root = origin_name)
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

make_range_stats = function(generation=0L) {
  tibble::tibble(
    generation,
    stat = rep(c('bimodality', 'shape_pda'), each=2L),
    value = c(0, 1, 2, -10))
}

label_value_tr = function(labels, multi_line = TRUE) {
  lapply(labels, function(x) {c(bimodality = 'BI', shape_pda = 'Shape stat. (PDA)')[x]})
}

ggplot_evolution = function(.tblstats, only_bi = FALSE) {
  .x = tidyr::gather(.tblstats, stat, value, -generation)
  .blank = make_range_stats(min(.x$generation))
  .hline = tibble::tibble(stat = c('bimodality', 'shape_pda'), value = c(5/9, NA))
  if (only_bi) {
    .x = .x %>% dplyr::filter(stat == "bimodality")
    .blank = .blank %>% dplyr::filter(stat == "bimodality")
    .hline = .hline %>% dplyr::filter(stat == "bimodality")
  }
  ggplot(.x, aes(generation, value))+
  geom_blank(data = .blank)+
  geom_line()+
  geom_hline(data = .hline, aes(yintercept=value), colour='red', linetype='dashed')+
  facet_grid(stat ~ ., scale='free_y', switch='y', label = label_value_tr)+
  theme_bw()+
  theme(
    panel.grid.minor = element_blank(),
    strip.placement = "outside",
    strip.background = element_blank(),
    strip.text = element_text(size = 12, family = 'Helvetica'),
  )
}
# .tblstats %>% ggplot_evolution()

# #######1#########2#########3#########4#########5#########6#########7#########
if (FALSE) {

main = function(indir, interval=2000L) {
  message(indir)
  label = fs::path_file(indir)
  read_fastas(indir, interval=interval) %>%
    add_phylo() %>%
    eval_treeshape() %>%
    ggplot_evolution()+
    labs(title=label)
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
