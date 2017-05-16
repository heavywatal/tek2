library(tidyverse)
library(pipeR)
library(wtl)

.tsv = read_tsv('activities.tsv')
.tsv

.p = .tsv %>>%
    dplyr::mutate(copies= copies / 500) %>>%
    ggplot(aes(time, copies, group=activity))+
    geom_area(aes(fill=activity))+
    wtl::scale_fill_gradientb('Spectral', reverse=TRUE)+
    wtl::theme_wtl()
.p
ggsave('activity.png', .p, width=10, height=7)
