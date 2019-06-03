# library(ggtree)
wtl::refresh("ggtree")

mrca_node = function(phylo) {
  root_node = which(phylo$tip.label == origin_name)
  edges = tibble::as_tibble(phylo$edge)
  root_parent = dplyr::filter(edges, V2 == root_node)$V1
  dplyr::filter(edges, V1 == root_parent, V1 != V2, V2 != root_node)$V2
}
# mrca_node(.phylo)

ggtree_hide_root = function(phylo, layout) {
  mrca = mrca_node(phylo)
  ggtree(ggtree::groupClade(phylo, .node=mrca), aes(alpha=group), layout=layout)+
  scale_alpha_manual(values=c("0" = 0.0, "1" = 0.7))
}

list_geom_tippoint = function(max_copy_number) {
  list(
    geom_tippoint(aes(colour=activity, size=copy_number), alpha=0.6),
    scale_colour_gradientn(colours = rev(head(rainbow(15L), 12L)), limits = c(0, 1), breaks = c(0, 0.5, 1)),
    scale_size(limit=c(1, max_copy_number), range=c(3, 12))
  )
}

ggtree_tek = function(phylo, data, individual="", layout="rectangular", max_copy_number=12, ...) {
  copy_number_origin = dplyr::filter(data, label == origin_name)$copy_number %>% sum()
  ggtr = if (copy_number_origin > 0L) {
    ggtree(phylo, layout=layout)
  } else {
    ggtree_hide_root(phylo, layout=layout)
  }
  ggtr %<+% data +
  list_geom_tippoint(max_copy_number)+
  geom_tippoint(aes(subset=is_major, size=copy_number), pch=1, colour="#000000")+
  geom_tippoint(aes(subset=is_fixed, size=copy_number), pch=4, colour="#000000")+
  labs(title = paste("sample", individual))
}

plot_individuals = function(infile, layout="rectangular") {
  message(infile)
  .seqs = read_tek_fasta(infile, metadata=TRUE)
  .mcols_all = tidy_metadata(.seqs)
  .nested_mcols = nest_metadata(.mcols_all, .seqs)
  .max_copy_number = dplyr::filter(.nested_mcols, individual == "total")$data[[1]]$copy_number %>% max()
  .plts = .nested_mcols %>% purrr::pmap(ggtree_tek, max_copy_number = .max_copy_number, layout = layout)
  cowplot::plot_grid(plotlist=.plts)
}

make_facettable = function(gt, dd, max_copy_number) {
  gt$data = dd %>%
    tidyr::nest(-individual) %>%
    {purrr::map2_dfr(.$data, .$individual, ~{
      (gt %<+% .x)$data %>% dplyr::mutate(individual = .y)
    })}
  gt + list_geom_tippoint(max_copy_number)
}


# #######1#########2#########3#########4#########5#########6#########7#########
if (FALSE) {

fortify_phylo_tbl = function(.tbl, layout = "rectangular") {
  purrr::pmap_df(.tbl, function(phylo, data, generation, individual, ...) {
    ggtree::fortify(phylo, layout = layout) %>%
      dplyr::left_join(data, by="label") %>%
      dplyr::mutate(generation = generation, individual = individual)
  })
}

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
  geom_tippoint(aes(subset=is_major, size=copy_number), pch=1, colour="#000000")+
  geom_tippoint(aes(subset=is_fixed, size=copy_number), pch=4, colour="#000000")+
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
plot_individuals(.fagz, layout="daylight")
Rprof(NULL)
summaryRprof()

fs::dir_ls(.focus$indir[2], regex="\\d+\\.fa\\.gz$")
.infiles = fs::path(.focus$indir[2], sprintf("generation_%05d.fa.gz", seq(0, 40000, by=4000)[-1]))
# .plts = purrr::map(.infiles, plot_individuals)
.plts = wtl::map_par(.infiles, plot_individuals)
.plts[[3]]
ggsave("individual_trees.pdf", .plts, width=10, height=10)

.fagz = fs::path(.focus$indir[2], "generation_16000.fa.gz")
.seqs = .fagz %>% read_tek_fasta(metadata=TRUE)
# .mcols = tidy_mcols(.seqs) %>% print()
# count_holders(.mcols)
.mcols_all = tidy_metadata(.seqs) %>% print()
.nested_mcols = .mcols_all %>% nest_metadata(.seqs) %>% print()
.max_copy_number = dplyr::filter(.nested_mcols, individual == "total")$data[[1]]$copy_number %>% max()

.plts = .nested_mcols %>%
  # head(1) %>%
  purrr::pmap(ggtree_tek, max_copy_number = .max_copy_number)
.cow = cowplot::plot_grid(plotlist=.plts)
.cow
ggsave("individual_trees_rect.png", .cow, width=10, height=10)

.nested_mcols$data[[2]]

.phylo = .nested_mcols$phylo[[1]]
.p = ggtree_hide_root(.phylo) + geom_tiplab()
.p
ggplot(.phylo)$data


# #######1#########2#########3#########4#########5#########6#########7#########
# draw sampled nodes on master tree

.useqs = unique_dss(.seqs)
.phylo = .useqs %>% Biostrings::stringDist(method="hamming") %>% ape::fastme.ols()
.rooted = ape::root(.phylo, "0x0", resolve.root=TRUE)
.gt = ggtree(.rooted, layout="rectangular", alpha=0.5)

make_facettable(.gt, .mcols_all, .max_copy_number)+
facet_wrap(~individual)

# #######1#########2#########3#########4#########5#########6#########7#########
# frequency spectrum

.mcols_holders = .mcols %>% count_holders()
.mcols_holders$copy_number %>% wtl::gghist()

.freqdf = tibble(
  num_holders = .mcols_holders$copy_number,
  num_copies = .mcols_total$copy_number) %>% print()

ggplot(.freqdf, aes(num_holders, num_copies))+
  geom_jitter(width=0.1, height=0.1, alpha=0.3)+
  geom_abline(intercept=0, slope=1)+
  theme_bw()

}
