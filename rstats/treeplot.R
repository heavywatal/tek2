# library(ggtree)
wtl::refresh('ggtree')

origin_seq = DNAStringSet(c('0x0' = str_dup("A", 300L)))
origin_name = names(origin_seq)

#' @param x BStringSet with metadata
tidy_mcols = function(x) {
  mcols(x) %>% as.data.frame() %>% as_tibble() %>%
    dplyr::select(label=te, everything())
}

summarise_mcols = function(.mcols, .id = 'total') {
  .mcols %>%
    group_by(label, activity, dn, ds, indel, species) %>%
    summarise(copy_number = sum(copy_number)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(individual = .id)
}

count_holders = function(.mcols, .id = 'holders') {
  .mcols %>%
    dplyr::count(label, activity, dn, ds, indel, species) %>%
    dplyr::rename(copy_number = n) %>%
    dplyr::arrange(desc(copy_number), label) %>%
    dplyr::mutate(individual = .id)
}

freq_in_samples = function(.mcols) {
  .mcols %>%
    count_holders() %>%
    dplyr::transmute(
      label,
      is_major = (copy_number > 5),
      is_fixed = (copy_number == max(copy_number))
    )
}

mrca_node = function(phylo) {
  root_node = which(phylo$tip.label == origin_name)
  edges = tibble::as_tibble(phylo$edge)
  root_parent = dplyr::filter(edges, V2 == root_node)$V1
  dplyr::filter(edges, V1 == root_parent, V1 != V2, V2 != root_node)$V2
}
# mrca_node(.phylo)

ggtree_hide_root = function(phylo, layout='rectangular') {
  mrca = mrca_node(phylo)
  ggtree(ggtree::groupClade(phylo, .node=mrca), aes(alpha=group), layout=layout) %>%
  ggtree::viewClade(mrca)+
  scale_alpha_manual(values=c('0' = 0.0, '1' = 1.0))
}

list_geom_tippoint = function(max_copy_number) {
  list(
    geom_tippoint(aes(colour=activity, size=copy_number), alpha=0.6),
    scale_colour_gradientn(colours = rev(head(rainbow(15L), 12L)), breaks = c(0, 0.5, 1)),
    scale_size(limit=c(1, max_copy_number), range=c(2, 12))
  )
}

ggtree_tek = function(phylo, data, individual='', max_copy_number=12, ...) {
  ggtree_hide_root(phylo, layout='rectangular') %<+% data +
  list_geom_tippoint(max_copy_number)+
  geom_tippoint(aes(subset=is_major, size=copy_number), pch=1, colour='#000000')+
  geom_tippoint(aes(subset=is_fixed, size=copy_number * 0.4), pch=4, colour='#000000')+
  labs(title = paste('sample', individual))
}

tidy_metadata = function(dss) {
  .mcols = tidy_mcols(dss)
  # validate uniqueness of TE address
  .mcols %>%
    dplyr::select(-individual, -copy_number) %>%
    dplyr::distinct() %>%
    {stopifnot(!any(.$label %>% duplicated()))}
  .inds = unique(.mcols$individual)
  .mcols = .mcols %>% add_row(label=origin_name, activity=1.0, copy_number=0L, dn=0.0, ds=0.0, indel=0L, individual=.inds, species=0L)
  .mcols_total = summarise_mcols(.mcols)
  .freq_cols = freq_in_samples(.mcols)
  .mcols %>%
    dplyr::bind_rows(.mcols_total) %>%
    # dplyr::bind_rows(.mcols_holders) %>%
    dplyr::left_join(.freq_cols, by='label')
}

unique_dss = function(dss) {
  dss = c(dss, origin_seq) %>% {.[!duplicated(names(.))]}
  mcols(dss) = NULL
  dss
}

nest_metadata = function(.tidy_metadata, dss) {
  .useqs = unique_dss(dss)
  .tidy_metadata %>%
    tidyr::nest(-individual) %>%
    dplyr::mutate(seqs = purrr::map(data, ~{.useqs[.x$label]})) %>%
    add_phylo(root = origin_name)
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

# #######1#########2#########3#########4
# layout

plot.phylo(.phylo, type='unrooted', show.tip.label=FALSE)

.daylight = ggtree(.phylo, layout='daylight')
.daylight

}
