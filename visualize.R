library(tidyverse)
library(wtl)

.indirs = wtl::command_args()$args
if (length(.indirs) < 1L) {
    .indirs = list.dirs(full.names=FALSE, recursive=FALSE)
}

.infiles = file.path(.indirs, 'activities.tsv') %>%
    purrr::keep(file.exists)

extract_params = function(filename, params=c('lambda', 'xi', 'nu')) {
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

.data %>%
    dplyr::mutate(copies= copies / 500) %>%
    ggplot(aes(time, copies, group=activity))+
    geom_area(aes(fill=activity))+
    facet_grid(desc(nu) ~ desc(xi) + lambda, labeller=label_both)+
    # wtl::scale_fill_gradientb('Spectral', reverse=TRUE)+
    scale_fill_gradientn(colours=rev(head(rainbow(15L), 12L)), breaks=c(0, 0.5, 1))+
    wtl::theme_wtl()+
    theme(legend.position='bottom')
