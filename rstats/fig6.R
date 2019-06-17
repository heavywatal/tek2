# @ILP22% ~/git/teaposon/run.py -p2 -j6 -r6 te2fig6
# rsync -auv ilp22.local:~/working/te2-fig6/ ~/working/tek/te2-fig6/

source("~/git/teaposon/rstats/read.R")
.indirs = fs::dir_ls(regexp = "^r\\d+", type = "directory")
.metadata = read_metadata(.indirs) %>% print()

df6act = .metadata %>%
  dplyr::mutate(
    label = sprintf("n=%d c=%d l=%d u=%d repl=%d", n, coexist, lower, upper, repl),
    adata = purrr::map(indir, read_activity),
    indir = NULL
  ) %>%
  tidyr::unnest() %>%
  dplyr::mutate(copy_number = copy_number / n) %>%
  print()

fig6candidates = df6act %>%
  dplyr::filter((generation %% 2000) == 0) %>%
  ggplot_activity() +
  facet_wrap(~label, ncol = 6L) +
  theme_minimal()
fig6candidates
ggsave("fig6candidates.pdf", fig6candidates, width = 7, height = 9.9, scale = 2)

fig6candidates_c8 = df6act %>%
  dplyr::filter(coexist == 8L) %>%
  dplyr::filter((generation %% 500) == 0) %>%
  ggplot_activity() +
  facet_wrap(~label, ncol = 6L) +
  theme_minimal()
fig6candidates_c8
ggsave("fig6candidates-c8.pdf", fig6candidates_c8, width = 9.9, height = 7, scale = 2)

.fig6repl = tibble::tribble(
     ~n, ~coexist, ~lower, ~upper, ~repl,
  # 1000L,       8L,     6L,    18L,    1L,
  1000L,       8L,     6L,    18L,    3L,
  1000L,       8L,     6L,    30L,    1L,
  1000L,       8L,     9L,    18L,    5L,
  1000L,       8L,     9L,    30L,    1L,
) %>% print()

df6 = .metadata %>%
  dplyr::right_join(.fig6repl, by = names(.fig6repl)) %>%
  dplyr::mutate(
    label = sprintf("n=%d l=%d u=%d repl=%d", n, lower, upper, repl),
    adata = purrr::map(indir, read_activity),
    indir = NULL
  ) %>%
  tidyr::unnest() %>%
  dplyr::mutate(copy_number = copy_number / n) %>%
  print()

fig6 = df6 %>%
  dplyr::rename(lowerbound = lower, upperbound = upper) %>%
  ggplot_activity() +
  facet_grid(upperbound ~ lowerbound, labeller = label_both) +
  theme(strip.background = element_blank())
fig6
ggsave("fig6.pdf", fig6, width = 10, height = 4)
