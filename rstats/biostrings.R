library(Biostrings)

parse_fasta_header = function(x) {
  str_match_all(x, '(\\w+)=(\\S+)') %>%
  purrr::map_dfr(~{
    tibble::tibble(key = .x[,2], value = .x[,3]) %>%
    tidyr::spread(key, value, convert=FALSE)
  }) %>%
  dplyr::mutate_at(vars(activity, dn, ds), as.double) %>%
  dplyr::mutate_at(vars(copy_number, indel, species), as.integer)
}

read_tek_fasta = function(file, metadata=FALSE, nrec=-1L, skip=0L) {
  .dss = Biostrings::readDNAStringSet(file, nrec=nrec, skip=skip)
  .names = names(.dss)
  if (metadata) {
    mcols(.dss) = parse_fasta_header(.names)
    names(.dss) = mcols(.dss)$te
  } else {
    names(.dss) = str_extract(.names, "(?<=te=)\\S+")
  }
  .dss
}

read_fastas = function(dir, interval = 1000L) {
  .fastas = fs::dir_ls(dir, regexp='generation_\\d+\\.fa\\.gz$')
  tibble::tibble(
    path = .fastas,
    infile = fs::path_file(path),
    generation = as.integer(readr::parse_number(infile))) %>%
  dplyr::filter((generation %% interval) == 0L) %>%
  dplyr::transmute(
    generation,
    seqs = purrr::map(path, read_tek_fasta)
  )
}
# .tbl = read_fastas('lower10_upper30_20180130T172206_00') %>% print()


origin_seq = Biostrings::DNAStringSet(c('0x0' = str_dup("A", 300L)))
origin_name = names(origin_seq)
sample_size = 10L

#' @param x BStringSet with metadata
tidy_mcols = function(x) {
  mcols(x) %>% as.data.frame() %>% as_tibble() %>%
    dplyr::select(label=te, everything()) %>%
    rename_origin()
}

rename_origin = function(.mcols) {
  old_label = .mcols %>%
    dplyr::distinct(label, dn, ds, indel, species) %>%
    dplyr::filter(dn == 0, ds == 0, indel == 0L, species == 0L) %>%
    purrr::pluck('label')
  if (is.null(old_label)) {
    .mcols
  } else {
    .mcols %>%
      dplyr::mutate(label = dplyr::recode(label, !!old_label := origin_name))
  }
}

summarise_mcols = function(.mcols, .id = 'total') {
  .mcols %>%
    group_by(label, activity, dn, ds, indel, species) %>%
    summarise(copy_number = sum(copy_number)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(individual = .id)
}

count_holders = function(.mcols, .id = 'holders') {
  .mcols %>%
    dplyr::group_by(label, activity, dn, ds, indel, species) %>%
    dplyr::summarise(copy_number = sum(copy_number > 0L)) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(desc(copy_number), label) %>%
    dplyr::mutate(individual = .id)
}

freq_in_samples = function(.mcols) {
  .mcols %>%
    count_holders() %>%
    dplyr::transmute(
      label,
      is_major = (copy_number > 5),
      is_fixed = (copy_number == sample_size)
    )
}

tidy_metadata = function(dss, add_root = TRUE) {
  .mcols = tidy_mcols(dss)
  # validate uniqueness of TE address
  .mcols %>%
    dplyr::select(-individual, -copy_number) %>%
    dplyr::distinct() %>%
    {stopifnot(!any(.$label %>% duplicated()))}
  .inds_with_origin = dplyr::filter(.mcols, label == origin_name)$individual
  .inds_wo_origin = unique(.mcols$individual) %>% {.[!. %in% .inds_with_origin]}
  if (add_root && length(.inds_wo_origin) > 0L) {
    .mcols = .mcols %>%
      add_row(label=origin_name, activity=1.0, copy_number=0L, dn=0.0, ds=0.0, indel=0L, individual=.inds_wo_origin, species=0L)
  }
  .mcols_total = summarise_mcols(.mcols)
  .freq_cols = freq_in_samples(.mcols)
  .mcols %>%
    dplyr::bind_rows(.mcols_total) %>%
    dplyr::left_join(.freq_cols, by='label')
}

unique_dss = function(dss) {
  dss = c(dss, origin_seq) %>% {.[!duplicated(names(.))]}
  mcols(dss) = NULL
  dss
}
