library(tidyverse)

.tbl_fitness = read_tsv("fitness.tsv.gz") %>% print()

N = max(dplyr::count(.tbl_fitness, generation)$n)
df = .tbl_fitness %>% dplyr::filter(generation >= 400L)
df_max = .max = df %>% dplyr::group_by(generation) %>% summarise(fitness = max(fitness))

.p = ggplot(df, aes(fitness, group = generation, colour = generation)) +
  geom_line(aes(y = stat(count)), stat = "bin", binwidth = 0.02, size = 1, alpha = 0.5) +
  geom_point(data = df_max, y = -10, size = 2, alpha = 0.5) +
  scale_colour_distiller(palette = "Spectral", direction = 1) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, N)) +
  labs(title = basename(getwd())) +
  theme_bw() + theme(legend.position = "bottom")
.p
ggsave("fitness_history.pdf", .p, width = 5, height = 5)
