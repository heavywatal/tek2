library(tidyverse)

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
