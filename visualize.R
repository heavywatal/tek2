library(tidyverse)
library(jsonlite)
library(wtl)
loadNamespace('cowplot')

# #######1#########2#########3#########4#########5#########6#########7#########

extract_params = function(filename, params=c('alpha', 'beta', 'lambda', 'xi', 'nu')) {
    patterns = sprintf('_%s([^_]+)_', params)
    str_match(paste0('_', filename), patterns)[,2] %>%
    parse_double() %>%
    set_names(params) %>%
    as.list() %>%
    as_tibble()
}
# .infiles[1] %>% extract_params()

read_history_json = function(path) {
    message(path)
    jsonlite::read_json(path, simplifyVector=TRUE) %>%
    {tibble(generation= as.integer(names(.)), data =.)} %>%
    dplyr::arrange(generation)
}
# .infiles[1] %>% read_history_json()

count_activity = function(x) {
    popsize = length(x) / 2L
    .names = c('site', 'indel', 'nonsynonymous', 'synonymous', 'activity')
    purrr::compact(x) %>%
    purrr::flatten_chr() %>%
    {tibble::tibble(tmpcol=.)} %>%
    tidyr::separate(tmpcol, .names, sep=':', convert=FALSE) %>%
    dplyr::group_by(activity= as.double(activity)) %>%
    dplyr::summarise(copy_number= n() / popsize)
}
# .histories$data[[4]] %>% count_activity()

count_per_individual = function(x) {
    lengths(x) %>%
    matrix(ncol=2L) %>%
    rowSums() %>%
    {tibble::tibble(copy_number=as.integer(.))} %>%
    dplyr::count(copy_number)
}
# .histories$data[[4]] %>% count_per_individual()

plot_copynumber_generation = function(data) {
    ggplot(data, aes(generation, copy_number, group=activity))+
    geom_area(aes(fill=activity))+
    scale_fill_gradientn(colours=rev(head(rainbow(15L), 12L)), breaks=c(0, 0.5, 1))+
    wtl::theme_wtl()+
    theme(legend.position='bottom')
}
# .counted %>% plot_copynumber_generation()

# #######1#########2#########3#########4#########5#########6#########7#########

.indirs = list.dirs(full.names=FALSE, recursive=FALSE)
.indirs = wtl::command_args()$args
.indirs = '.'

.infiles = file.path(.indirs, 'history.json.gz') %>%
    purrr::keep(file.exists) %>%
    print()

.histories = .infiles %>% set_names() %>%
    map_df(extract_params, .id='infile') %>%
    dplyr::mutate(data= purrr::map(infile, read_history_json)) %>%
    tidyr::unnest() %>%
    print()

.counted = .histories %>%
    dplyr::mutate(data= purrr::map(data, count_activity)) %>%
    tidyr::unnest() %>%
    print()

.nested = .counted %>%
    tidyr::nest(-xi) %>%
    dplyr::mutate(plt= purrr::map2(data, paste0('xi = ', xi), ~{
        plot_copynumber_generation(.x)+
        facet_grid(alpha + desc(nu) ~ lambda, labeller=label_both)+
        labs(title=.y)
    }))

.p = cowplot::plot_grid(plotlist=.nested$plt, nrow=1)
.p
.nested$plt[[1]] + coord_cartesian(ylim=c(0, 370))

ggsave('fig-s2.png', .p, width=15, height=12)


.per_individual = .histories %>%
    dplyr::mutate(data= purrr::map(data, count_per_individual)) %>%
    unnest() %>%
    print()

.p = .per_individual %>%
    dplyr::filter(0L <= generation, generation <= 2000L) %>%
    ggplot(aes(copy_number, n, group=generation))+
    geom_line(aes(colour=generation), size=1, alpha=0.8)+
    wtl::scale_colour_gradientb('Spectral', 4L)+
    theme_bw()+theme(legend.position='bottom')
.p
ggsave('copies_per_individual.pdf', .p, width=5, height=5)

# #######1#########2#########3#########4#########5#########6#########7#########

.tbl_fitness = read_tsv('fitness.tsv.gz') %>% print()
.p = .tbl_fitness %>%
    dplyr::filter(generation >= 400L) %>%
    {
      .max = . %>% group_by(generation) %>% summarise(fitness=max(fitness))
      ggplot(., aes(fitness, group=generation, colour=generation))+
      geom_line(aes(y=..count..), stat='bin', binwidth=0.02, size=1, alpha=0.5)+
      geom_point(data=.max, y=-10, size=2, alpha=0.5)+
      wtl::scale_colour_gradientb('Spectral', 4L)+
      coord_cartesian(xlim=c(0, 1), ylim=c(0, 500))+
      labs(title=getwd() %>% basename())+
      theme_bw()+theme(legend.position='bottom')
    }
.p
ggsave('fitness_history.pdf', .p, width=5, height=5)

# #######1#########2#########3#########4#########5#########6#########7#########

w = c(0.04, 0.16)
n = 500L
r = 5000L

.fast = map_int(seq_len(r), ~{
    v = sample(w,  2L * n, replace=TRUE)
    v[v > max(w) * runif(2L * n)] %>% head(n) %>% table() %>% head(1)
}) %>% {tibble(p= . / n, method='fast')}
.honest = map_int(seq_len(r), ~{
    v = sample(w, 12L * n, replace=TRUE)
    v[v > runif(12L * n)] %>% head(n) %>% table() %>% head(1)
}) %>% {tibble(p= . / n, method='honest')}

.p = bind_rows(.fast, .honest) %>% {
    ggplot(., aes(p, group=method, fill=method))+
    geom_density(alpha=0.5)+
    theme_bw()
}
.p
ggsave('selection_methods.pdf', .p, width=4, height=4)

# #######1#########2#########3#########4#########5#########6#########7#########

.cols = c('position', 'indel', 'nonsynonymous', 'synonymous', 'activity')
.nsam = 100L
.gametes = 'summary.json.gz' %>% fromJSON()

.summary = .gametes %>% sample(.nsam) %>% flatten_chr() %>%
    {tibble(x=.)} %>%
    tidyr::separate(x, .cols, ':') %>%
    dplyr::mutate_all(as.numeric) %>%
    dplyr::mutate(dn= nonsynonymous / 200, ds= synonymous / 100) %>%
    dplyr::mutate(dn_ds= dn / ds, distance= (nonsynonymous + synonymous) / 300) %>%
    print()

.summary %>%
    ggplot(aes(dn_ds, group=activity))+
    geom_histogram(aes(fill=activity), binwidth=0.1, center=0)+
    scale_fill_gradientn(colours=rev(head(rainbow(15L), 12L)), breaks=c(0, 0.5, 1))+
    labs(x= 'dN/dS', y= 'Copy number')+
    theme_bw()

.summary %>%
    ggplot(aes(distance, dn_ds))+
    geom_jitter(aes(colour= activity), alpha=0.3)+
    scale_colour_gradientn(colours=rev(head(rainbow(15L), 12L)), breaks=c(0, 0.5, 1))+
    geom_hline(aes(yintercept= 1.0), linetype='dashed')+
    coord_cartesian(xlim=c(0, 1))+
    labs(x= 'Distance from original sequence', y= 'dN/dS')+
    theme_bw()

.summary %>%
    dplyr::count(position) %>%
    dplyr::mutate(freq = n / .nsam) %>%
    ggplot(aes(freq))+
    geom_histogram(binwidth=0.1, center=0)+
    coord_cartesian(xlim=c(0, 1))+
    labs(x= 'TE frequency', y= 'Site number')+
    theme_bw()

# #######1#########2#########3#########4#########5#########6#########7#########

library(Biostrings)
library(ape)

.fa = Biostrings::readBStringSet('sequence.fa.gz') %>%
    {setNames(., str_extract(names(.), '^\\S+'))} %>%
    print()

.mat = .fa %>% sample(16L) %>% stringDist('hamming')
.cl = hclust(.mat, 'single'); plot(.cl)
.cl = hclust(.mat, 'complete'); plot(.cl)
.cl = hclust(.mat, 'average'); plot(.cl)
.cl = hclust(.mat, 'mcquitty'); plot(.cl)

.tree = nj(.mat)
plot(.tree, type='u')
