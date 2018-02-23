library(tidyverse)
library(jsonlite)
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

plot_copynumber_generation = function(data) {
  ggplot(data, aes(generation, copy_number, group = activity)) +
    geom_area(aes(fill = activity)) +
    scale_fill_gradientn(colours = rev(head(rainbow(15L), 12L)), breaks = c(0, 0.5, 1)) +
    wtl::theme_wtl() +
    theme(legend.position = "bottom")
}
# .counted %>% plot_copynumber_generation()

# #######1#########2#########3#########4#########5#########6#########7#########

popsize = 500
.indirs = list.dirs(full.names = FALSE, recursive = FALSE)
.indirs = wtl::command_args()$args
.indirs = "."

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

read_activity = function(indir) {
  file.path(indir, "activity.tsv.gz") %>%
  read_tsv() %>%
  dplyr::mutate(copy_number = copy_number / popsize)
}

ggplot_activity = function(data) {
  dplyr::mutate(data, species = factorize_species(species)) %>%
  ggplot(aes(generation, copy_number)) +
  geom_area(aes(group = interaction(activity, species), fill = activity), position = position_stack(reverse = FALSE)) +
  scale_fill_gradientn(colours = rev(head(rainbow(15L), 12L)), breaks = c(0, 0.5, 1)) +
  wtl::theme_wtl() +
  theme(legend.position = "none")
}

factorize_species = function(x) {
  factor(x, levels=sort.int(unique(x), decreasing=TRUE))
}
# .metadata$indir[[1]] %>% read_activity() %>% ggplot_activity()

source('~/git/tek-evolution/treestats.R')

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

# library(ggtree)
wtl::refresh('ggtree')

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
    dplyr::mutate(individual = .id)
}

sort_by_individual = function(.mcols, .useqs) {
  .mcols %>%
    tidyr::nest(-individual) %>%
    dplyr::mutate(seqs = purrr::map(data, ~{.useqs[.x$label]}))
}

.focus = .metadata %>% dplyr::filter(xi == max(xi), lower == 9, upper == 24) %>% print()

.fagz = fs::path(.focus$indir[1], "generation_30000.fa.gz")
.seqs = .fagz %>% read_tek_fasta(metadata=TRUE)
.mcols = tidy_mcols(.seqs) %>% print()
.mcols_total = summarise_mcols(.mcols) %>% print()
.useqs = .seqs %>% {.[!duplicated(names(.))]} %>% {.[.mcols_total$label]} %>% print()
mcols(.useq) = NULL
.phylo = .useqs %>% Biostrings::stringDist(method='hamming') %>% ape::fastme.ols()
.gt = ggtree(.phylo, layout='equal_angle', alpha=0.5)
.max_copy_number = max(.mcols_total$copy_number)

.mcols_holders = .mcols %>% count_holders() %>% print()
.mcols_all = .mcols %>%
  dplyr::bind_rows(.mcols_all, .mcols_holders) %>%
  print()

# validate uniqueness of TE address
.mcols %>%
  dplyr::select(-individual, -copy_number) %>%
  dplyr::distinct() %>%
  {stopifnot(!any(.$label %>% duplicated()))}

.mcols_holders$copy_number %>% wtl::gghist()

.freqdf = tibble(
  num_holders = .mcols_holders$copy_number,
  num_copies = .mcols_total$copy_number) %>% print()

ggplot(.freqdf, aes(num_holders, num_copies))+
  geom_jitter(width=0.1, height=0.1, alpha=0.3)+
  geom_abline(intercept=0, slope=1)+
  theme_bw()

.tippoint = list(
  geom_tippoint(aes(colour=activity, size=copy_number), alpha=0.6),
  scale_colour_gradientn(colours = rev(head(rainbow(15L), 12L)), breaks = c(0, 0.5, 1)),
  scale_size(limit=c(1, max_copy_number), range=c(2, 12))
)

make_facettable = function(gt, dd) {
  gt$data = dd %>%
    tidyr::nest(-individual) %>%
    {purrr::map2_dfr(.$data, .$individual, ~{
      (gt %<+% .x)$data %>% dplyr::mutate(individual = .y)
    })}
  gt
}

make_facettable(.gt, .mcols_all)+
  .tippoint+
  facet_wrap(~individual)


ggtree_tek = function(data, individual, gt, max_copy_number, ...) {
  gt %<+% data +
  .tippoint+
  labs(title = paste('sample', individual))
}

.nested_mcols = .mcols_all %>%
  tidyr::nest(-individual) %>%
  print()

.plts = .nested_mcols %>% purrr::pmap(ggtree_tek, gt=.gt, max_copy_number=.max_copy_number)
.cow = cowplot::plot_grid(plotlist=.plts)
.cow
ggsave('individual_trees.png', .cow, width=10, height=10)

plot.phylo(.phylo, type='unrooted', show.tip.label=FALSE)

.daylight = ggtree(.phylo, layout='daylight')
.daylight

# #######1#########2#########3#########4#########5#########6#########7#########

.tbl_act = .metadata %>%
  # dplyr::filter(!(xi < 6e-4 & upper < 16)) %>%
  dplyr::filter(lower < 300, upper < 300) %>%
  # sample_n(6L) %>%
  dplyr::mutate(adata = purrr::map(indir, read_activity)) %>%
  tidyr::unnest() %>%
  print()

.p = ggplot_activity(.tbl_act)+
  facet_grid(xi * lower ~ upper * repl)
.p
ggsave("copynumber-activity-species.png", .p, width = 10, height = 10)

.p = .counted %>%
  dplyr::distinct(lower, upper, repl, generation, species) %>%
  dplyr::count(lower, upper, repl, species) %>%
  dplyr::mutate(n = n * 50L) %>%
  ggplot(aes(species, n)) +
  geom_col() +
  facet_grid(lower ~ upper * repl) +
  labs(x = "species ID", y = "sojourn time") +
  wtl::theme_wtl()
.p
ggsave("species-sojourn-time.png", .p, width = 7, height = 7)

.nested = .counted %>%
  tidyr::nest(-xi) %>%
  dplyr::mutate(plt = purrr::map2(data, paste0("xi = ", xi), ~{
    plot_copynumber_generation(.x) +
      facet_grid(alpha + desc(nu) ~ lambda, labeller = label_both) +
      labs(title = .y)
  }))

.p = cowplot::plot_grid(plotlist = .nested$plt, nrow = 1)
.p
.nested$plt[[1]] + coord_cartesian(ylim = c(0, 370))

ggsave("fig-s2.png", .p, width = 15, height = 12)

# #######1#########2#########3#########4#########5#########6#########7#########

.tbl_fitness = read_tsv("fitness.tsv.gz") %>% print()
.p = .tbl_fitness %>%
  dplyr::filter(generation >= 400L) %>%
  {
    .max = . %>% group_by(generation) %>% summarise(fitness = max(fitness))
    ggplot(., aes(fitness, group = generation, colour = generation)) +
      geom_line(aes(y = ..count..), stat = "bin", binwidth = 0.02, size = 1, alpha = 0.5) +
      geom_point(data = .max, y = -10, size = 2, alpha = 0.5) +
      wtl::scale_colour_gradientb("Spectral", 4L) +
      coord_cartesian(xlim = c(0, 1), ylim = c(0, 500)) +
      labs(title = getwd() %>% basename()) +
      theme_bw() + theme(legend.position = "bottom")
  }
.p
ggsave("fitness_history.pdf", .p, width = 5, height = 5)

# #######1#########2#########3#########4#########5#########6#########7#########

w = c(0.04, 0.16)
n = 500L
r = 5000L

.fast = map_int(seq_len(r), ~{
  v = sample(w, 2L * n, replace = TRUE)
  v[v > max(w) * runif(2L * n)] %>% head(n) %>% table() %>% head(1)
}) %>% {
  tibble(p = . / n, method = "fast")
}
.honest = map_int(seq_len(r), ~{
  v = sample(w, 12L * n, replace = TRUE)
  v[v > runif(12L * n)] %>% head(n) %>% table() %>% head(1)
}) %>% {
  tibble(p = . / n, method = "honest")
}

.p = bind_rows(.fast, .honest) %>% {
  ggplot(., aes(p, group = method, fill = method)) +
    geom_density(alpha = 0.5) +
    theme_bw()
}
.p
ggsave("selection_methods.pdf", .p, width = 4, height = 4)

# #######1#########2#########3#########4#########5#########6#########7#########

.cols = c("site", "species", "indel", "nonsynonymous", "synonymous", "activity")
.nsam = 100L
.gametes = "summary.json.gz" %>% fromJSON()

.summary = .gametes %>%
  sample(.nsam) %>%
  {
    tibble(tmpcol = ., gamete = seq_along(.))
  } %>%
  tidyr::unnest() %>%
  tidyr::separate(tmpcol, .cols, ":") %>%
  dplyr::mutate_at(vars(c("site", "species", "indel", "nonsynonymous", "synonymous")), as.integer) %>%
  dplyr::mutate(activity = as.double(activity)) %>%
  dplyr::mutate(dn = nonsynonymous / 200, ds = synonymous / 100) %>%
  dplyr::mutate(dn_ds = dn / ds, distance = (nonsynonymous + synonymous) / 300) %>%
  print()

.summary %>%
  ggplot(aes(dn_ds, group = activity)) +
  geom_histogram(aes(fill = activity), binwidth = 0.1, center = 0) +
  scale_fill_gradientn(colours = rev(head(rainbow(15L), 12L)), breaks = c(0, 0.5, 1)) +
  labs(x = "dN/dS", y = "Copy number") +
  theme_bw()

.summary %>%
  ggplot(aes(distance, dn_ds)) +
  geom_jitter(aes(colour = activity), alpha = 0.3) +
  scale_colour_gradientn(colours = rev(head(rainbow(15L), 12L)), breaks = c(0, 0.5, 1)) +
  geom_hline(aes(yintercept = 1.0), linetype = "dashed") +
  coord_cartesian(xlim = c(0, 1)) +
  labs(x = "Distance from original sequence", y = "dN/dS") +
  theme_bw()

.summary %>%
  dplyr::count(site) %>%
  dplyr::mutate(freq = n / .nsam) %>%
  ggplot(aes(freq)) +
  geom_histogram(binwidth = 0.1, center = 0) +
  coord_cartesian(xlim = c(0, 1)) +
  labs(x = "TE frequency", y = "Site number") +
  theme_bw()
