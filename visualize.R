library(tidyverse)
library(wtl)
loadNamespace('cowplot')

.indirs = wtl::command_args()$args
if (length(.indirs) < 1L) {
    .indirs = list.dirs(full.names=FALSE, recursive=FALSE)
}

.infiles = file.path(.indirs, 'activities.tsv') %>%
    purrr::keep(file.exists)

extract_params = function(filename, params=c('alpha', 'beta', 'lambda', 'xi', 'nu')) {
    patterns = sprintf('_%s([^_]+)_', params)
    str_match(paste0('_', filename), patterns)[,2] %>%
    parse_double() %>%
    set_names(params) %>%
    as.list() %>%
    as_tibble()
}
extract_params(.infiles[1])

.data = .infiles %>% set_names() %>%
    map_df(extract_params, .id='infile') %>%
    dplyr::mutate(data= purrr::map(infile, read_tsv)) %>%
    tidyr::unnest() %>% print()

plot_copynumber_time = function(.data, .title='') {
    ggplot(.data, aes(time, copies, group=activity))+
    geom_area(aes(fill=activity))+
    facet_grid(alpha + desc(nu) ~ lambda, labeller=label_both)+
    # wtl::scale_fill_gradientb('Spectral', reverse=TRUE)+
    scale_fill_gradientn(colours=rev(head(rainbow(15L), 12L)), breaks=c(0, 0.5, 1))+
    labs(title=.title)+
    wtl::theme_wtl()+
    theme(legend.position='bottom')
}

.nested = .data %>%
    dplyr::mutate(copies= copies / 500) %>%
    tidyr::nest(-xi) %>%
    dplyr::mutate(plt= purrr::map2(data, paste0('xi = ', xi), plot_copynumber_time))

.p = cowplot::plot_grid(plotlist=.nested$plt, nrow=1)
.p

ggsave('fig-s2.png', .p, width=15, height=12)
