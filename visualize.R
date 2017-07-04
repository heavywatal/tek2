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

.indirs = '.'
.indirs = list.dirs(full.names=FALSE, recursive=FALSE)
.indirs = wtl::command_args()$args

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

ggsave('fig-s2.png', .p, width=15, height=12)


.per_individual = .histories %>%
    dplyr::mutate(data= purrr::map(data, count_per_individual)) %>%
    unnest() %>%
    print()

.per_individual %>%
    ggplot(aes(copy_number, n, group=generation))+
    geom_line(aes(colour=generation))+
    theme_bw()
