if (FALSE) {

.indirs = fs::dir_ls(regexp = "^n\\d+", type = "directory")
.metadata = read_metadata(.indirs) %>% print()

fig1df = .metadata %>%
  dplyr::filter(n == 1000L) %>%
  dplyr::mutate(
    adata = wtl::mcmap(indir, read_activity),
    tdata = wtl::mcmap(indir, ~{
      message(.x)
      read_fastas(.x, interval = 100L) %>%
        add_phylo() %>%
        eval_treeshape()
    })
  ) %>%
  print()

saveRDS(fig1df, "fig1.rds")
# rsync -auv ilp15.local:~/working/te2-fig1/ ~/working/tek/te2-fig1/

}
# #######1#########2#########3#########4#########5#########6#########7#########

source("~/git/teaposon/rstats/read.R")

fig1df = readRDS("fig1.rds")

fig1df$tdata[[1]] %>% ggplot_evolution()

fig1plts = fig1df %>%
  dplyr::arrange(desc(xi)) %>%
  dplyr::mutate(
    tplot = purrr::imap(tdata, ~{
      .xmin = min(.x$generation)
      .p = ggplot_evolution(.x %>% dplyr::filter(generation > 500L), only_bi = TRUE) +
        coord_cartesian(xlim = c(.xmin, max(.x$generation))) +
        theme(
          axis.title = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank()
        )
      if (.y == 1L) {
        .p = .p +
          annotate("rect", xmin = 25000L, xmax = 45000L, ymin = 0, ymax = 1, fill = NA, colour = "#666666", linetype = "dashed") +
          annotate("text", label = "Fig. 2", x = 40500, y = 0.05, hjust = 0, vjust = 0)
      }
      .p
    }),
    aplot = purrr::map2(adata, n, ~{
      .p = ggplot_activity(.x, .y)
      if (max(.x$copy_number) / .y < 100) {
        .p = .p + theme(axis.title.y = element_text(margin = margin(r = 11.5)))
      }
      .p
    }),
    plt = purrr::map2(tplot, aplot, ~{
      cowplot::plot_grid(.x, .y + theme(legend.position = "none", axis.title.x = element_blank()),
        ncol = 1, rel_heights = c(1, 1), align = "v", axis = "lr")
    })
  ) %>%
  print()
fig1plts$plt[[1]]

fig1 = cowplot::plot_grid(
  cowplot::plot_grid(
    cowplot::plot_grid(plotlist = fig1plts$plt, ncol = 1, labels = c("A", "B", "C"), scale = 0.95),
    cowplot::get_legend(fig1plts$aplot[[1]]),
    nrow = 1L, rel_widths = c(1, 0.1)
  ),
  cowplot::ggdraw() + cowplot::draw_label("Generation"),
  ncol = 1L, rel_heights = c(3, 0.1)
)
ggsave("fig1.pdf", fig1, width = 10, height = 10)
