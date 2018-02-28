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
