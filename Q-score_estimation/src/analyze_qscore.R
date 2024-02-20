library(tidyverse)
library(Biostrings)
read_stats <- function(path){
  kit <- str_extract(path, "V14|R9")
  algorithm <- str_extract(path, "HAC|SUP")
  barcode <- str_extract(path, "barcode\\d\\d")
  
  col_types <- readr::cols(.default = readr::col_character())
  
  stats <- read_tsv(path, col_types = col_types) %>% 
    mutate(
      kit = kit, 
      algorithm = algorithm, 
      barcode = barcode
    )
}
read_fastqs <- function(path){
  kit <- str_extract(path, "V14|R9")
  algorithm <- str_extract(path, "HAC|SUP")
  barcode <- str_extract(path, "barcode\\d\\d")
  reads <- readQualityScaledDNAStringSet(path, use.names = TRUE)
  quals <- as(quality(reads), "IntegerList")
  
  fastq <- lapply(seq_along(reads), convert_to_df, reads, quals) %>% 
    bind_rows() %>% 
    mutate(
      kit = kit, 
      algorithm = algorithm, 
      barcode = barcode
    )
}
read_fastqs_per_pos <- function(path, error_stats_filtered){
  kit <- str_extract(path, "V14|R9")
  algorithm <- str_extract(path, "HAC|SUP")
  barcode <- str_extract(path, "barcode\\d\\d")
  reads <- readQualityScaledDNAStringSet(path, use.names = TRUE)
  quals <- as(quality(reads), "IntegerList")
  
  fastq <- lapply(seq_along(reads), convert_to_df_per_pos, reads, quals, error_stats_filtered) %>% 
    bind_rows() %>% 
    mutate(
      kit = kit, 
      algorithm = algorithm, 
      barcode = barcode,
      cluster_id = str_sub(cluster_id, end = -3)
    ) %>% 
    inner_join(error_stats_filtered, by = c("barcode", "kit", "algorithm", "cluster_id", "pos" = "read_position"))
}
convert_to_df_per_pos <- function(i, reads, quals, error_stats_filtered){
  read <- reads[i]
  df <- tibble(
    pos = 1:width(read),
    base = str_split(as.character(reads[i]), pattern = "")[[1]],
    qual = quals[[i]],
    cluster_id = names(read)
  )
}
convert_to_df <- function(i, reads, quals){
  read <- reads[i]
  df <- tibble(
    qual = mean(quals[[i]]),
    cluster_id = names(read)
  )
}

qscore_paths <- list.files("stats_context/",
                            pattern = "cluster_qscore_stats.tsv", 
                            recursive = TRUE, 
                            full.names = TRUE)

error_paths <- list.files("stats_context/",
                            pattern = "cluster_error_stats.tsv", 
                            recursive = TRUE, 
                            full.names = TRUE)


cluster_stats_paths <- list.files("stats/",
                                  pattern = "split_cluster_stats.tsv",
                                  recursive = TRUE, 
                                  full.names = TRUE)

fastq_paths <- list.files("fastq/",
                          pattern = ".fastq",
                          full.names = TRUE,
                          recursive = TRUE
                          )

max_qscore = 40
qscore_paths_filtered = qscore_paths[str_detect(qscore_paths, "consensus")]
error_paths_filtered = error_paths[str_detect(error_paths, "consensus")]
fastq_paths_filtered = fastq_paths[str_detect(fastq_paths, "consensus")]


qscore_stats <- lapply(qscore_paths_filtered, read_stats) %>% 
  bind_rows() %>% 
  type_convert() %>% 
  mutate(
    cluster_id = str_sub(cluster, end = -3),
    qscore = ifelse(qscore == 60, max_qscore, qscore)
  )

cluster_stats_raw <-
  lapply(cluster_stats_paths, read_stats) %>%
  bind_rows() %>% 
  type_convert()

cluster_stats <- 
  cluster_stats_raw %>% 
  filter(cluster_written == 1) %>% 
  mutate(
    reads_written = reads_written_fwd + reads_written_rev
  )

error_stats <- lapply(error_paths_filtered, read_stats) %>% 
  bind_rows() %>% 
  type_convert() %>% 
  mutate(
    cluster_id = str_sub(cluster, end = -3),
    qscore = ifelse(qscore == 60, max_qscore, qscore))

fastqs <- lapply(fastq_paths_filtered, read_fastqs) %>% 
  bind_rows() %>% 
  mutate(
    cluster_id = str_sub(cluster_id, end = -3),
  ) %>% 
  dplyr::rename(medaka_quality = qual)

qscore_stats_merged <- qscore_stats %>% 
  inner_join(cluster_stats, by = c("barcode", "kit", "algorithm", "cluster_id")) %>% 
  inner_join(fastqs, by = c("barcode", "kit", "algorithm", "cluster_id"))

qscore_stats_merged_filtered <- qscore_stats_merged %>% 
  filter(qscore != max_qscore)
  
error_stats_merged <- error_stats %>% 
  inner_join(cluster_stats, by = c("barcode", "kit", "algorithm", "cluster_id"))

systematic_errors <- error_stats_merged %>%
  filter(error_type != "-") %>% 
  group_by(kit, algorithm, fragment, sample_type, position, error_type, read_base, reference_base) %>% 
  summarise(error_count = n()) %>% 
  filter(error_count > 1)

systematic_errors_per_cluster_size <- error_stats_merged %>% 
  group_by(kit, algorithm, fragment, reads_written, position, error_type, read_base, reference_base) %>% 
  summarise(error_count = n()) %>% 
  filter(error_count > 1)

systematic_error_types <- 
  error_stats_merged %>%
  filter(reads_written > 4) %>% 
  filter(str_detect(error_type, "-", negate = TRUE)) %>%
  mutate(
    base_change = paste0(toupper(reference_base), toupper(read_base))
  ) %>% 
  group_by(kit, algorithm, fragment, error_type, base_change) %>% 
  summarise(n = n())

systematic_error_types_pos <- 
  error_stats_merged %>% 
  filter(reads_written > 8) %>% 
  filter(kit == "V14") %>% 
  filter(str_detect(error_type, "-", negate = TRUE)) %>%
  mutate(
    base_change = paste0(toupper(reference_base), toupper(read_base)),
    position = as.integer(str_remove(position, "\\.\\d+"))
  ) %>%
  group_by(kit, algorithm, fragment, sample_type, position, base_change) %>% 
  summarise(n = n())

error_stats_filtered <- 
  error_stats_merged %>% 
  mutate(
    base_change = paste0(toupper(reference_base), toupper(read_base))
  ) %>% 
  filter(kit == "V14" & algorithm == "SUP" & str_detect(base_change, "CA|GT")) 

n_total_errors_filtered <- 
  error_stats_filtered %>% 
  group_by(kit, algorithm, base_change) %>% 
  summarise(n_total = n())

error_stats_filtered_annotated <- 
  error_stats_filtered %>% 
  mutate(context = paste0(rev_context, fwd_context)) %>% 
  group_by(kit, algorithm, base_change, context, rev_context, fwd_context) %>%
  summarise(n = n()) %>% 
  inner_join(n_total_errors_filtered) %>% 
  mutate(n_relative = n / n_total)

error_stats_filtered_annotated %>% 
  group_by(kit, algorithm, base_change) %>% 
  summarise(percent = sum(n_relative))


# sliding through context
test <- lapply(seq(1, 5, 1), function(window_size){
  lapply(seq(1, window_size, 1), function(i, window_size){
    start = 7 - i
    end = start - 1 + window_size
    error_stats_filtered_annotated %>% 
      mutate(context = str_sub(context, start = start, end = end)) %>% 
      group_by(kit, algorithm, base_change, context) %>% 
      summarise(n = sum(n) / window_size, 
                n_relative = sum(n_relative) / window_size) %>%
      mutate(start = start, 
             end = end, 
             window_size = window_size)
  }, window_size) %>% 
    bind_rows()
}) %>% bind_rows()

test %>%
  filter(base_change == "CA") %>% 
  filter(context != "A") %>% 
  filter(n_relative > 0.03) %>%
  mutate(
    facet_cols = nchar(context),
  ) %>% 
  ggplot(aes(context, n_relative, fill = base_change)) +
  facet_grid(
    cols = vars(facet_cols), 
    # rows = vars(facet_rows), 
    scales = "free_x") +
  geom_bar(position="dodge", stat="identity") + 
  theme(
    axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), 
    legend.position = "bottom"
  )

test %>%
  filter(base_change == "GT") %>% 
  filter(context != "T") %>% 
  filter(n_relative > 0.03) %>%
  mutate(
    facet_cols = nchar(context),
  ) %>% 
  ggplot(aes(context, n_relative, fill = base_change)) +
  facet_grid(
    cols = vars(facet_cols), 
    # rows = vars(facet_rows), 
    scales = "free_x") +
  geom_bar(position="dodge", stat="identity") + 
  theme(
    axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), 
    legend.position = "bottom"
  )

window_size <- 1
test_fwd <- lapply(seq(1, 5 - window_size, 1), function(start, window_size){
  lapply(seq(start + window_size, 5, 1), function(end, start, window_size){
    error_stats_filtered_annotated %>% 
      mutate(context = str_sub(fwd_context, start = start, end = end)) %>% 
      group_by(kit, algorithm, base_change, context) %>% 
      summarise(n = sum(n), 
                n_relative = sum(n_relative)) %>%
      mutate(start = start, 
             end = end, 
             window_size = window_size, 
             direction = "fwd")
  }, start, window_size) %>% 
    bind_rows()
}, window_size) %>% bind_rows()
test_rev <- lapply(seq(5, 1 + window_size, -1), function(end, window_size){
  lapply(seq(end - window_size, 1, -1), function(start, end, window_size){
    error_stats_filtered_annotated %>% 
      mutate(context = str_sub(rev_context, start = start, end = end)) %>% 
      group_by(kit, algorithm, base_change, context) %>% 
      summarise(n = sum(n), 
                n_relative = sum(n_relative)) %>%
      mutate(start = start, 
             end = end, 
             window_size = window_size, 
             direction = "rev")
  }, end, window_size) %>% 
    bind_rows()
}, window_size) %>% bind_rows()

# increasing size of sequence context
test_fwd <- lapply(seq(1, 5, 1), function(end){
    error_stats_filtered_annotated %>% 
      mutate(context = str_sub(fwd_context, end = end)) %>% 
      group_by(kit, algorithm, base_change, context) %>% 
      summarise(n = sum(n), 
                n_relative = sum(n_relative)) %>%
      mutate(direction = "fwd")
}) %>% bind_rows()
test_rev <- lapply(seq(5, 1, -1), function(start){
    error_stats_filtered_annotated %>% 
      mutate(context = str_sub(rev_context, start = start)) %>% 
      group_by(kit, algorithm, base_change, context) %>% 
      summarise(n = sum(n), 
                n_relative = sum(n_relative)) %>%
      mutate(direction = "rev")
}) %>% bind_rows()

test <- 
  test_rev %>%
  bind_rows(test_fwd)



test %>%
  # filter(n_relative > 0.02) %>% 
  filter(abs(start - end) == 1) %>% 
  mutate(
    facet_cols = base_change,
  ) %>% 
  ggplot(aes(context, n_relative, fill = direction)) +
  facet_grid(cols = vars(facet_cols), scales = "free_x") +
  geom_bar(position="dodge", stat="identity") + 
  theme(
    axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), 
    legend.position = "bottom"
    )

test %>%
  filter(n_relative > 0.02) %>%
  mutate(
    facet_rows = nchar(context),
    facet_cols = base_change,
  ) %>% 
  ggplot(aes(context, n_relative, fill = direction)) +
  facet_grid(cols = vars(facet_cols), rows = vars(facet_rows), scales = "free_x") +
  geom_bar(position="dodge", stat="identity") + 
  theme(
    axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), 
    legend.position = "bottom"
  )

test %>%
  filter(base_change == "CA") %>% 
  filter(context != "A") %>% 
  filter(n_relative > 0.02) %>%
  mutate(
    facet_cols = nchar(context),
    # facet_cols = base_change,
  ) %>% 
  ggplot(aes(context, n_relative, fill = direction)) +
  facet_grid(
    cols = vars(facet_cols),
    # rows = vars(facet_rows),
    scales = "free_x") +
  geom_bar(position="dodge", stat="identity") + 
  theme(
    axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), 
    legend.position = "bottom"
    )

test %>%
  filter(n_relative > 0.04) %>%
  filter(abs(start - end) == 2) %>% 
  mutate(
    facet_cols = base_change,
  ) %>% 
  ggplot(aes(context, n_relative, fill = direction)) +
  facet_grid(cols = vars(facet_cols), scales = "free_x") +
  geom_bar(position="dodge", stat="identity") + 
  theme(
    axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), 
    legend.position = "bottom"
  )


error_stats_filtered_annotated %>% 
  filter(base_change == "GT") %>%
  mutate(context = str_sub(context, start = 1, end = 3)) %>%
  group_by(kit, algorithm, base_change, context) %>% 
  summarise(n = sum(n), 
            n_relative = sum(n_relative)) %>% 
  ggplot(aes(context, n_relative, fill = base_change)) +
  geom_bar(position="dodge", stat="identity", show.legend = FALSE) + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

error_stats_filtered %>% 
  filter(base_change == "CA") %>% 
  mutate(context = paste0(rev_context, fwd_context)) %>% 
  group_by(kit, algorithm, base_change, context) %>%
  summarise(n = n(), 
            n_relative = n / n_total_errors_filtered_CA * 100) %>%
  filter(n_relative > 0.1) %>%
  ggplot(aes(context, n_relative, fill = base_change)) +
  facet_grid(rows = vars(base_change)) +
  geom_bar(position="dodge", stat="identity", show.legend = FALSE) + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

error_stats_filtered %>% 
  mutate(fwd_context = str_sub(fwd_context, end = 2)) %>% 
  group_by(kit, algorithm, base_change, fwd_context) %>%
  summarise(n = n()) %>%
  ggplot(aes(fwd_context, n_relative, fill = base_change)) +
  geom_bar(position="dodge", stat="identity") + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

error_stats_filtered %>% 
  group_by(kit, algorithm, base_change, rev_context) %>%
  summarise(n = n(), 
            n_relative = n / n_total_errors_filtered * 100) %>% 
  ggplot(aes(rev_context, n_relative, fill = base_change)) +
  facet_grid(rows = vars(base_change)) +
  geom_bar(position="dodge", stat="identity", show.legend = FALSE) + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))


stats_all_clusters <- lapply(seq(2,100,2), FUN = function(min_cluster_size, qscore_stats_merged){
  qscore_stats_merged %>% 
    filter(reads_written >= min_cluster_size) %>% 
    group_by(kit, algorithm, fragment, sample_type) %>% 
    summarise(errors = sum(n_errors), 
              n_cluster = n()) %>% 
    mutate(
      n_nucleotides = n_cluster * fragment,
      error_probability = errors / n_nucleotides,
      qscore = -10 * log10(error_probability),
      qscore = ifelse(is.infinite(qscore ), max_qscore, qscore),
      min_cluster_size = min_cluster_size
    )
}, qscore_stats_merged) %>% 
  bind_rows()
stats_clusters_categories <- lapply(seq(2,100,2), FUN = function(min_cluster_size, qscore_stats_merged){
  qscore_stats_merged %>% 
    filter(reads_written >= min_cluster_size) %>% 
    group_by(kit, algorithm, fragment, sample_type) %>% 
    mutate(
      min_cluster_size = min_cluster_size
    )
}, qscore_stats_merged) %>% 
  bind_rows()

error_stats_categorized <- lapply(seq(2,100,2), FUN = function(cluster_threshold, stats){
  stats %>% 
    filter(reads_written >= cluster_threshold) %>% 
    mutate(
      cluster_threshold = cluster_threshold
    )
}, error_stats_merged) %>% 
  bind_rows()

clusters_without_error <- 
  error_stats_categorized %>% 
  filter(str_detect(error_type, "-")) %>%
  select(kit, algorithm, sample_type, fragment, cluster_threshold, cluster) %>% 
  mutate(
    n_errors_per_cluster = 0,
    n_errorrs_per_cluster_categorized = factor(0)
  )

clusters_with_error <- 
  error_stats_categorized %>% 
  filter(str_detect(error_type, "-", negate = TRUE)) %>%
  group_by(kit, algorithm, sample_type, fragment, cluster_threshold, cluster) %>%
  summarise(n_errors_per_cluster = n()) %>% 
  mutate(
    n_errorrs_per_cluster_categorized =
      cut(
        n_errors_per_cluster, 
        breaks = c(-Inf, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, Inf), 
        labels = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, "10+"))
  )

n_errors_per_cluster_categorized <- 
  clusters_with_error %>% 
  bind_rows(clusters_without_error)

n_errors_per_pos_categorized <- 
  error_stats_categorized %>% 
  filter(str_detect(error_type, "-", negate = TRUE)) %>% 
  mutate(
    position = as.integer(str_remove(position, "\\.\\d+"))
  ) %>% 
  group_by(kit, algorithm, sample_type, fragment, cluster_threshold, position) %>% 
  summarise(n_errors_per_pos = n()) %>% 
  mutate(
    n_errorrs_per_pos_categorized =
      cut(
        n_errors_per_pos, 
        breaks = c(-Inf, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, Inf), 
        labels = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, "10+"))
  )

  
n_errors_per_cluster_threshold_categorized <- 
  error_stats_categorized %>% 
  filter(str_detect(error_type, "-", negate = TRUE)) %>% 
  group_by(kit, algorithm, sample_type, fragment, cluster_threshold) %>% 
  summarise(total_errors = n()) 

n_errors_per_error_type_categorized <- 
  error_stats_categorized %>% 
  filter(str_detect(error_type, "-", negate = TRUE)) %>% 
  mutate(
    base_change = paste0(toupper(reference_base), read_base)
  ) %>% 
  group_by(kit, algorithm, sample_type, fragment, cluster_threshold, error_type, base_change) %>%
  summarise(n_errors_per_error_type = n())

n_errors_per_error_type_categorized_annotated <- 
  n_errors_per_error_type_categorized %>% 
  inner_join(n_errors_per_cluster_threshold_categorized) %>% 
  group_by(kit, algorithm, sample_type, fragment, cluster_threshold, error_type, base_change) %>% 
  mutate(
    relative_error_rate_per_error_type = n_errors_per_error_type / total_errors
  )

n_errors_per_cluster_size <- 
  error_stats_merged %>% 
  filter(str_detect(error_type, "-", negate = TRUE)) %>% 
  group_by(kit, algorithm, sample_type, fragment, reads_written) %>% 
  summarise(total_errors = n()) 

n_errors_per_error_type_size <- 
  error_stats_merged %>% 
  filter(str_detect(error_type, "-", negate = TRUE)) %>% 
  mutate(
    base_change = paste0(toupper(reference_base), read_base)
  ) %>% 
  group_by(kit, algorithm, sample_type, fragment, reads_written, error_type, base_change) %>%
  summarise(n_errors_per_error_type = n())

n_errors_per_error_type_size <- 
  n_errors_per_error_type_size %>% 
  inner_join(n_errors_per_cluster_size) %>% 
  group_by(kit, algorithm, sample_type, fragment, reads_written, error_type, base_change) %>% 
  mutate(
    relative_error_rate_per_error_type = n_errors_per_error_type / total_errors
  )

n_errors_per_error_type_size %>% 
  group_by(kit, algorithm, sample_type, fragment, reads_written) %>% 
  summarise(percent = sum(relative_error_rate_per_error_type))

error_stats_categorized_annotated <- 
  n_errors_per_cluster_threshold_categorized %>% 
  inner_join(n_errors_per_cluster_categorized)


stats_summary <- 
  qscore_stats_merged %>% 
  group_by(kit, algorithm, fragment, reads_written, sample_type) %>% 
  summarise(
    errors = sum(n_errors), 
    n_cluster = n(), 
    medaka_quality = mean(medaka_quality)) %>% 
  mutate(
    n_nucleotides = n_cluster * fragment,
    error_probability = errors / n_nucleotides,
    qscore = -10 * log10(error_probability),
  )


stats_summary_cluster_threshold  <- 
  stats_clusters_categories %>% 
  group_by(kit, algorithm, fragment, min_cluster_size) %>% 
  summarise(
    errors = sum(n_errors), 
    n_cluster = n(), 
    medaka_quality = mean(medaka_quality)) %>% 
  mutate(
    n_nucleotides = n_cluster * fragment,
    error_probability = errors / n_nucleotides,
    qscore = -10 * log10(error_probability),
  )


write_tsv(qscore_stats_merged, file = "out/tables/qscore_cluster_medaka_qual_stats.tsv")
write_tsv(error_stats_merged, file = "out/tables/error_cluster_stats.tsv")
write_tsv(stats_summary, "out/tables/qualities_per_cluster_size.tsv")
write_tsv(stats_summary_cluster_threshold, "out/tables/qualities_per_cluster_threshold.tsv")


errors_per_pos <- error_stats_merged %>%
  filter(str_detect(error_type, "-", negate = TRUE)) %>% 
  mutate(
    position = as.integer(str_remove(position, "\\.\\d+"))
  ) %>% 
  group_by(kit, algorithm, fragment, sample_type, position, reads_written) %>% 
  summarise(n_errors = n())

errors_per_pos <- error_stats_merged %>%
  filter(kit == "V14") %>% 
  filter(reads_written > 6) %>% 
  filter(str_detect(error_type, "-", negate = TRUE)) %>% 
  mutate(
    position = as.integer(str_remove(position, "\\.\\d+"))
  ) %>% 
  group_by(kit, algorithm, fragment, sample_type, position) %>% 
  summarise(n_errors = n())

systematic_errors_per_pos <- error_stats_merged %>%
  filter(reads_written > 8) %>% 
  filter(str_detect(error_type, "-", negate = TRUE)) %>% 
  mutate(
    position = as.integer(str_remove(position, "\\.\\d+"))
  ) %>% 
  group_by(kit, algorithm, fragment, sample_type, position) %>% 
  summarise(n_errors = n()) %>% 
  filter(n_errors > 1)

systematic_errors <- 
  systematic_errors_per_pos %>% 
  group_by(position, fragment, sample_type) %>% 
  summarise(occurence = n()) %>% 
  filter(occurence > 1)

# cluster analysis
cluster_stats_raw %>% 
  group_by(kit, algorithm, barcode) %>% 
  summarise(n = n()) %>% 
  ggplot(aes(paste(kit, algorithm, barcode), n, fill = paste(kit, barcode))) + 
  geom_bar(position="dodge", stat="identity") + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

ggsave(
  filename = "number of clusters.jpeg", 
  device ="jpeg",
  path ="out/plots/",
  width = 12, 
  height = 9, 
  dpi = 400
)

cluster_stats_raw %>%
  filter(cluster_written == 1) %>% 
  group_by(kit, algorithm, barcode) %>% 
  summarise(n = n()) %>% 
  ggplot(aes(paste(kit, algorithm, barcode), n, fill = paste(kit, barcode))) + 
  geom_bar(position="dodge", stat="identity") + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

ggsave(
  filename = "number of clusters written.jpeg", 
  device ="jpeg",
  path ="out/plots/",
  width = 12, 
  height = 9, 
  dpi = 400
)

# cluster size
{
  stats_summary %>% 
    filter(!is.infinite(qscore)) %>% 
    mutate(
      facet_cols = paste(kit, algorithm),
      facet_rows = paste(fragment),
      color = paste(kit, algorithm)
    ) %>% 
    ggplot() + 
    facet_grid(rows = vars(facet_rows), cols = vars(facet_cols)) +
    geom_point(aes(reads_written, qscore, color = color), size = 0.8) +
    scale_x_continuous(breaks = seq(0,100,8)) +
    scale_y_continuous(breaks = seq(0,80,5)) +
    labs(
      x = "cluster size",
      y = "Q-score"
    ) +
    theme(
      legend.position = "bottom", 
      text = element_text(size = 14),
      axis.text = element_text(size = 12)
    )
  
  ggsave(
    filename = "Q_score_per_cluster_size_without_medaka.jpg",
    path = "out/plots/",
    dpi = 600,
    width = 12, 
    height = 9
  )
}
{
  stats_summary %>% 
    filter(!is.infinite(qscore)) %>% 
    mutate(
      facet_cols = paste(kit, algorithm),
      facet_rows = paste(fragment),
      color = paste(kit, algorithm)
    ) %>% 
    ggplot() + 
    facet_grid(rows = vars(facet_rows), cols = vars(facet_cols)) +
    geom_point(aes(reads_written, qscore, color = color), size = 0.8) +
    geom_point(aes(reads_written, medaka_quality), color = "grey", size = 0.8) +
    scale_x_continuous(breaks = seq(0,100,8)) +
    scale_y_continuous(breaks = seq(0,80,5)) +
    labs(
      x = "cluster size",
      y = "Q-score"
    ) +
    theme(
      legend.position = "bottom", 
      text = element_text(size = 14),
      axis.text = element_text(size = 12)
    )
  
  ggsave(
    filename = "Q_score_per_cluster_size_with_medaka.jpg",
    path = "out/plots/",
    dpi = 600,
    width = 12, 
    height = 9
  )
}
{
  qscore_stats_merged %>% 
    mutate(
      facet_cols = paste(kit, algorithm),
      facet_rows = paste(fragment, sample_type),
      color = paste(kit, algorithm)
    ) %>% 
    ggplot() + 
    facet_grid(rows = vars(facet_rows), cols = vars(facet_cols)) +
    geom_point(aes(reads_written, qscore, color = color), size = 0.8) +
    geom_point(aes(reads_written, medaka_quality), color = "grey", size = 0.8) + 
    scale_x_continuous(breaks = seq(0,100,8)) +
    scale_y_continuous(breaks = seq(0,80,5)) +
    labs(
      x = "cluster size",
      y = "Q-score"
    ) +
    theme(
      legend.position = "bottom", 
      text = element_text(size = 14),
      axis.text = element_text(size = 12)
    )
  
  ggsave(
    filename = "Q_score_per_cluster_size_per_sample_type_with_medaka.jpg",
    path = "out/plots/",
    dpi = 600,
    width = 12, 
    height = 9
  )
}
{
  qscore_stats_merged %>% 
  mutate(
    facet_cols = paste(kit, algorithm),
    facet_rows = paste(fragment),
    color = paste(kit, algorithm)
  ) %>% 
  ggplot() + 
  facet_grid(rows = vars(facet_rows), cols = vars(facet_cols)) +
  geom_point(aes(reads_written, qscore, color = color), size = 0.8) +
  geom_point(aes(reads_written, medaka_quality), color = "grey", size = 0.8) + 
  scale_x_continuous(breaks = seq(0,100,8)) +
  scale_y_continuous(breaks = seq(0,80,5)) +
  labs(
    x = "cluster size",
    y = "Q-score"
  ) +
  theme(
    legend.position = "bottom", 
    text = element_text(size = 14),
    axis.text = element_text(size = 12)
  )

ggsave(
  filename = "Q_score_per_cluster_size_with_medaka.jpg",
  path = "out/plots/",
  dpi = 600,
  width = 12, 
  height = 9
)
}
{
  qscore_stats_merged %>% 
  filter(reads_found < 20) %>% 
  mutate(
    facet_cols = paste(kit, algorithm),
    facet_rows = paste(fragment),
    color = paste(kit, algorithm)
  ) %>% 
  ggplot() + 
  facet_grid(rows = vars(facet_rows), cols = vars(facet_cols)) +
  geom_point(aes(reads_found, qscore, color = color), size = 0.8) +
  geom_point(aes(reads_found, medaka_quality), color = "grey", size = 0.8) + 
  scale_x_continuous(breaks = seq(0,100,2)) +
  scale_y_continuous(breaks = seq(0,80,5)) +
  labs(
    x = "cluster size",
    y = "Q-score"
  ) +
  theme(
    legend.position = "bottom",
    text = element_text(size = 14)
  )

ggsave(
  filename = "Q_score_per_cluster_size_with_medaka_zoomed.jpg",
  path = "out/plots/",
  dpi = 600,
  width = 12, 
  height = 9
)
}
{
  qscore_stats_merged %>% 
  filter(reads_found < 20) %>% 
  mutate(
    facet_cols = paste(kit, algorithm),
    facet_rows = paste(fragment),
    color = paste(kit, algorithm),
    group = reads_found
  ) %>% 
  ggplot(aes(reads_found, qscore)) + 
  facet_grid(rows = vars(facet_rows), cols = vars(facet_cols)) +
  geom_boxplot(aes(color = color, group = group), size = 0.5) + 
  scale_x_continuous(breaks = seq(0,100,2)) +
  scale_y_continuous(breaks = seq(0,80,5)) +
  labs(
    x = "cluster size",
    y = "Q-score"
  ) +
  theme(
    legend.position = "bottom", 
    text = element_text(size = 14),
    axis.text = element_text(size = 12)
  )

ggsave(
  filename = "Q_score_boxplots_per_cluster_size_zoomed.jpg",
  path = "out/plots/",
  dpi = 600,
  width = 12, 
  height = 9
)
}
{
  qscore_stats_merged %>% 
  filter(reads_found < 20) %>% 
  mutate(
    facet_cols = paste(kit, algorithm),
    facet_rows = paste(fragment),
    color = paste(kit, algorithm),
    group = reads_found
  ) %>% 
  ggplot() + 
  facet_grid(rows = vars(facet_rows), cols = vars(facet_cols)) +
  geom_boxplot(aes(reads_found, qscore, color = color, group = group), size = 0.5) + 
  geom_boxplot(aes(reads_found, medaka_quality, group = group), color = "grey", size = 0.5) + 
  scale_x_continuous(breaks = seq(0,100,2)) +
  scale_y_continuous(breaks = seq(0,80,5)) +
  labs(
    x = "cluster size",
    y = "Q-score"
  ) +
  theme(
    legend.position = "bottom",
    text = element_text(size = 14)
  )

ggsave(
  filename = "Q_score_boxplots_per_cluster_size_with_medaka_zoomed.jpg",
  path = "out/plots/",
  dpi = 600,
  width = 12, 
  height = 9
)
}
{
  qscore_stats_merged %>% 
  mutate(
    facet_cols = paste(kit, algorithm),
    facet_rows = paste(fragment),
    color = paste(kit, algorithm),
    group = reads_written
  ) %>% 
  ggplot(aes(reads_written, qscore)) + 
  facet_grid(rows = vars(facet_rows), cols = vars(facet_cols)) +
  geom_boxplot(aes(color = color, group = group), size = 0.5) + 
  scale_x_continuous(breaks = seq(0,100,8)) +
  scale_y_continuous(breaks = seq(0,80,5)) +
  labs(
    x = "cluster size",
    y = "Q-score"
  ) +
  theme(
    legend.position = "bottom", 
    text = element_text(size = 14),
    axis.text = element_text(size = 12)
  )

ggsave(
  filename = "Q_score_boxplots_per_cluster_size.jpg",
  path = "out/plots/",
  dpi = 600,
  width = 12, 
  height = 9
)
}

# cluster threshold
{
  stats_summary_cluster_threshold %>% 
    filter(!is.infinite(qscore)) %>% 
    mutate(
      facet_cols = paste(kit, algorithm),
      facet_rows = paste(fragment),
      color = paste(kit, algorithm)
    ) %>% 
    ggplot() + 
    facet_grid(rows = vars(facet_rows), cols = vars(facet_cols)) +
    geom_point(aes(min_cluster_size, qscore, color = color), size = 0.8) +
    geom_point(aes(min_cluster_size, medaka_quality), color = "grey", size = 0.8)  + 
    scale_x_continuous(breaks = seq(0,100,8)) +
    scale_y_continuous(breaks = seq(0,80,5)) +
    labs(
      x = "cluster threshold",
      y = "Q-score"
    ) +
    theme(
      legend.position = "bottom", 
      text = element_text(size = 14),
      axis.text = element_text(size = 12)
    )
  
  ggsave(
    filename = "Q_score_per_cluster_threshold_with_medaka.jpg",
    path = "out/plots/",
    dpi = 600,
    width = 12, 
    height = 9
  )
}
{
  stats_summary_cluster_threshold %>% 
    filter(!is.infinite(qscore)) %>% 
    mutate(
      facet_cols = paste(kit, algorithm),
      facet_rows = paste(fragment),
      color = paste(kit, algorithm)
    ) %>% 
    ggplot() + 
    facet_grid(rows = vars(facet_rows), cols = vars(facet_cols)) +
    geom_point(aes(min_cluster_size, qscore, color = color), size = 0.8) +
    scale_x_continuous(breaks = seq(0,100,8)) +
    scale_y_continuous(breaks = seq(0,80,5)) +
    labs(
      x = "cluster threshold",
      y = "Q-score"
    ) +
    theme(
      legend.position = "bottom", 
      text = element_text(size = 14),
      axis.text = element_text(size = 12)
    )
  
  ggsave(
    filename = "Q_score_per_cluster_threshold_without_medaka.jpg",
    path = "out/plots/",
    dpi = 600,
    width = 12, 
    height = 9
  )  
}
{
  stats_clusters_categories %>% 
    mutate(
      facet_cols = paste(kit, algorithm),
      facet_rows = paste(fragment),
      color = paste(kit, algorithm),
      group = min_cluster_size
    ) %>% 
    ggplot(aes(min_cluster_size, qscore)) + 
    facet_grid(rows = vars(facet_rows), cols = vars(facet_cols)) +
    geom_boxplot(aes(color = color, group = group), size = 0.5) + 
    scale_x_continuous(breaks = seq(0,100,8)) +
    scale_y_continuous(breaks = seq(0,80,5)) +
    labs(
      x = "cluster threshold",
      y = "Q-score"
    ) +
    theme(
      legend.position = "bottom", 
      text = element_text(size = 14),
      axis.text = element_text(size = 12)
    )
  
  ggsave(
    filename = "Q_score_boxplots_per_cluster_threshold.jpg",
    path = "out/plots/",
    dpi = 600,
    width = 12, 
    height = 9
  )
}
{
  stats_clusters_categories_grouped <- 
    stats_clusters_categories  %>% 
    filter(min_cluster_size < 15) %>% 
    mutate(
      facet_cols = paste(kit, algorithm),
      facet_rows = paste(fragment),
      color = paste(kit, algorithm),
      group = min_cluster_size
    )  
  n_stats_clusters_categories_grouped <- 
    stats_clusters_categories_grouped %>% 
    group_by(kit, algorithm, fragment, min_cluster_size) %>% 
    summarise(number_of_cluster = n())
  
  stats_clusters_categories_grouped_annotated <- 
    stats_clusters_categories_grouped %>% 
    inner_join(n_stats_clusters_categories_grouped, by = c("kit", "algorithm", "fragment", "min_cluster_size"))
  
  n_stats_clusters_categories_grouped %>% 
    write_tsv("out/tables/qualities_per_cluster_threshold_zoomed.tsv")
  
  stats_clusters_categories_grouped_annotated %>% 
    ggplot(aes(min_cluster_size, qscore)) + 
    facet_grid(rows = vars(facet_rows), cols = vars(facet_cols)) +
    geom_point(aes(color = color, group = group), size = 0.01, position = "jitter") + 
    scale_x_continuous(breaks = seq(0,100,2)) +
    scale_y_continuous(breaks = seq(0,80,5)) +
    labs(
      x = "cluster threshold",
      y = "Q-score"
    ) +
    theme(
      legend.position = "bottom", 
      text = element_text(size = 14)
    )
  
  ggsave(
    filename = "Q_score_point_per_cluster_threshold_zoomed.jpg",
    path = "out/plots/",
    dpi = 600,
    width = 12, 
    height = 9
  )
  
  stats_clusters_categories_grouped_annotated %>% 
    ggplot(aes(min_cluster_size, qscore)) + 
    facet_grid(rows = vars(facet_rows), cols = vars(facet_cols)) +
    geom_boxplot(aes(color = color, group = group)) + 
    scale_x_continuous(breaks = seq(0,100,2)) +
    scale_y_continuous(breaks = seq(0,80,5)) +
    labs(
      x = "cluster threshold",
      y = "Q-score"
    ) +
    theme(
      legend.position = "bottom", 
      text = element_text(size = 14)
    )

  ggsave(
    filename = "Q_score_boxplots_per_cluster_threshold_zoomed.jpg",
    path = "out/plots/",
    dpi = 600,
    width = 12, 
    height = 9
  )
}

lapply(seq(2,30,2), FUN = function(minimal_cluster_size, stats, filter){
  lapply(unique(stats$kit), function(kit_filter, stats, minimal_cluster_size, filter){
    if (str_detect(filter, "threshold")) {
      stats_filtered <- 
        stats %>% 
        filter(reads_written >= minimal_cluster_size)
    } else {
      stats_filtered <- 
        stats %>% 
        filter(reads_written == minimal_cluster_size)
    }
    
    number_of_cluster <- 
      stats_filtered %>%
      filter(kit == kit_filter) %>% 
      nrow()
    
    stats_filtered_parsed <- 
      stats_filtered %>%
      dplyr::filter(kit == kit_filter) %>% 
      dplyr::filter(str_detect(error_type, "-", negate = TRUE)) %>% 
      dplyr::group_by(kit, algorithm, fragment, sample_type, position) %>% 
      dplyr::summarise(n_errors = n() / number_of_cluster) %>% 
      dplyr::ungroup() %>% 
      dplyr::mutate(
        facet_cols = paste(kit, algorithm),
        facet_rows = paste(fragment, sample_type),
        color = paste(kit, algorithm),
        position = as.integer(str_remove(position, "\\.\\d*"))
      )
    
    stats_filtered_parsed %>% 
      ggplot(aes(position, n_errors, color = color)) + 
      facet_grid(rows = vars(facet_rows), cols = vars(facet_cols)) +
      geom_point(size = 0.5) +
      scale_y_continuous(breaks = seq(0, 0.008, 0.001)) +
      scale_x_continuous(breaks = seq(0,5000, 400)) +
      coord_cartesian(ylim = c(0, 0.008))
      labs(
        x = "Position", 
        y = "Number of errors",
        title = paste("error count for", filter, minimal_cluster_size), 
      ) + 
      theme(
        legend.position = "bottom",
        text = element_text(size = 12)
      )
    
    ggsave(filename = paste0(kit_filter, "_", minimal_cluster_size, ".jpg"),
           path = paste0("out/plots/error_profile/", filter),
           device = "jpg",
           width = 12, 
           height = 9, 
           dpi = 400
    )
  }, stats, minimal_cluster_size, filter)
}, error_stats_merged, "threshold")


n_errors_per_error_type_size %>% 
  filter(reads_written < 30) %>%
  # filter(str_detect(base_change, ":") | base_change %in% c("CA", "GT", "CT", "GA")) %>% 
  mutate(
    facet_cols = paste(kit, algorithm),
    facet_rows = paste(fragment, sample_type),
    fill = factor(reads_written, levels = seq(2, 100, 2))
  ) %>% 
  ggplot(aes(base_change, relative_error_rate_per_error_type, fill = fill)) + 
  facet_grid(rows = vars(facet_rows), cols = vars(facet_cols), scales = "free") +
  geom_bar(position="dodge", stat="identity") +
  scale_y_continuous(breaks = seq(0, 1, 0.1)) +
  labs(
    x = "error type", 
    y = "relative occurence of error type",
    fill = "cluster size"
  ) +
  theme(
    legend.position = "bottom"
  )

ggsave(
  filename = "Relative abundance of error type per cluster size.jpg", 
  device = "jpg", 
  path = "out/plots/error_profile/size/", 
  width = 12, 
  height = 9, 
  dpi = 600)

n_errors_per_error_type_size %>% 
  filter(reads_written < 20) %>%
  filter(str_detect(base_change, ":") | base_change %in% c("CA", "GT", "CT", "GA")) %>% 
  mutate(
    facet_cols = paste(kit, algorithm),
    facet_rows = paste(fragment, sample_type),
    fill = factor(reads_written, levels = seq(2, 100, 2))
  ) %>% 
  ggplot(aes(base_change, relative_error_rate_per_error_type, fill = fill)) + 
  facet_grid(rows = vars(facet_rows), cols = vars(facet_cols), scales = "free") +
  geom_bar(position="dodge", stat="identity") +
  scale_y_continuous(breaks = seq(0, 1, 0.1)) +
  labs(
    x = "error type", 
    y = "relative occurence of error type",
    fill = "cluster size"
  ) +
  theme(
    legend.position = "bottom"
  )

ggsave(
  filename = "Relative abundance of error type per cluster size selected.jpg", 
  device = "jpg", 
  path = "out/plots/error_profile/size/", 
  width = 12, 
  height = 9, 
  dpi = 600)


n_errors_per_error_type_categorized_annotated %>% 
  filter(cluster_threshold <= 30) %>% 
  mutate(
    facet_cols = paste(kit, algorithm),
    facet_rows = paste(fragment, sample_type),
    fill = factor(cluster_threshold, levels = seq(2, 30, 2))
  ) %>% 
  ggplot(aes(base_change, relative_error_rate_per_error_type, fill = fill)) + 
  facet_grid(rows = vars(facet_rows), cols = vars(facet_cols), scales = "free") +
  geom_bar(position="dodge", stat="identity") +
  scale_y_continuous(breaks = seq(0, 1, 0.1)) +
  labs(
    x = "error type", 
    y = "relative occurence of error type",
    fill = "cluster threshold"
  ) +
  theme(
    legend.position = "bottom"
  )

ggsave(
  filename = "Relative abundance of error type per cluster threshold.jpg", 
  device = "jpg", 
  path = "out/plots/error_profile/threshold/", 
  width = 12, 
  height = 9, 
  dpi = 600)

n_errors_per_error_type_categorized_annotated %>% 
  filter(cluster_threshold <= 30) %>% 
  filter(str_detect(base_change, ":") | base_change %in% c("CA", "GT", "CT", "GA")) %>% 
  mutate(
    facet_cols = paste(kit, algorithm),
    facet_rows = paste(fragment, sample_type),
    fill = factor(cluster_threshold, levels = seq(2, 30, 2))
  ) %>% 
  ggplot(aes(base_change, relative_error_rate_per_error_type, fill = fill)) + 
  facet_grid(rows = vars(facet_rows), cols = vars(facet_cols), scales = "free") +
  geom_bar(position="dodge", stat="identity") +
  scale_y_continuous(breaks = seq(0, 1, 0.1)) +
  labs(
    x = "error type", 
    y = "relative occurence of error type",
    fill = "cluster threshold"
  ) +
  theme(
    legend.position = "bottom"
  )

ggsave(
  filename = "Relative abundance of error type per cluster threshold selected.jpg", 
  device = "jpg", 
  path = "out/plots/error_profile/threshold/", 
  width = 12, 
  height = 9, 
  dpi = 600)

error_stats_categorized_annotated %>% 
  # filter(kit == "V14") %>% 
  filter(cluster_threshold < 30) %>% 
  mutate(
    facet_cols = paste(kit, algorithm),
    facet_rows = paste(fragment, sample_type),
    color = paste(kit, algorithm)
  ) %>% 
  ggplot(aes(cluster_threshold, fill = factor(n_errorrs_per_cluster_categorized, levels = c(seq(0,10,1), "10+")))) + 
  facet_grid(rows = vars(facet_rows), cols = vars(facet_cols), scales = "free") +
  geom_bar(position = "fill", alpha = 0.8) +
  scale_x_continuous(breaks = seq(0, 40, 2)) +
  theme_bw() +
  theme(
    legend.position = "bottom"
  ) +
  labs(
    fill = "Number of errors",
    x = "Cluster threshold", 
    y = "relative abundance of errors per cluster"
  )

ggsave(
  filename = "Relative abundance of errors per cluster threshold.jpg", 
  device = "jpg", 
  path = "out/plots/error_profile/threshold/", 
  width = 12, 
  height = 9, 
  dpi = 600)

qscore_stats_merged %>% 
  filter(kit == "V14") %>% 
  filter(reads_written <= 100) %>% 
  mutate(
    n_errors_categorized =
      cut(
        n_errors, 
        breaks = c(-Inf, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, Inf), 
        labels = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, "10+"))
  ) %>% 
  mutate(
    facet_cols = paste(kit, algorithm),
    facet_rows = paste(fragment, sample_type),
    color = paste(kit, algorithm)
  ) %>% 
  ggplot(aes(reads_written, fill = factor(n_errors_categorized))) + 
  facet_grid(rows = vars(facet_rows), cols = vars(facet_cols), scales = "free") +
  geom_bar(position = "fill", alpha = 0.8) +
  scale_x_continuous(breaks = seq(2, 100, 4)) +
  theme_bw() +
  theme(
    legend.position = "bottom"
  ) +
  labs(
    fill = "Number of errors",
    x = "Cluster size", 
    y = "relative abundance of errors per cluster"
  )

ggsave(
  filename = "Relative abundance of errors per cluster size.jpg", 
  device = "jpg", 
  path = "out/plots/error_profile/size/", 
  width = 12, 
  height = 9, 
  dpi = 600)

n_errors_per_pos_categorized %>%
  filter(cluster_threshold <= 10 | cluster_threshold == 20 | cluster_threshold == 30) %>% 
  filter(kit == "R9" & algorithm == "HAC") %>% 
  mutate(
    facet_cols = paste(fragment, sample_type),
    facet_rows = factor(cluster_threshold, levels = seq(2, 30, 2)),
    color = paste(fragment, sample_type)
  ) %>% 
  ggplot(aes(position, n_errorrs_per_pos_categorized)) + 
  facet_grid(rows = vars(facet_rows), cols = vars(facet_cols), scale = "free_x") +
  geom_point(size = 0.2) +
  scale_x_continuous(breaks = seq(0,5000,1000))  +
  labs(
    y = "Number of errors", 
    x = "Position"
  )

ggsave(
  filename = "Relative abundance of errors per Pos per threshold R9 HAC.jpg", 
  device = "jpg", 
  path = "out/plots/error_profile/threshold/", 
  width = 12, 
  height = 9, 
  dpi = 600)

n_errors_per_pos_categorized %>%
  filter(cluster_threshold <= 10 | cluster_threshold == 20 | cluster_threshold == 30) %>% 
  filter(kit == "V14" & algorithm == "SUP") %>% 
  mutate(
    facet_cols = paste(fragment, sample_type),
    facet_rows = factor(cluster_threshold, levels = seq(2, 30, 2)),
    color = paste(fragment, sample_type)
  ) %>% 
  ggplot(aes(position, n_errorrs_per_pos_categorized)) + 
  facet_grid(rows = vars(facet_rows), cols = vars(facet_cols), scale = "free_x") +
  geom_point(size = 0.2) +
  scale_x_continuous(breaks = seq(0,5000,1000))  +
  labs(
    y = "Number of errors", 
    x = "Position"
  )

ggsave(
  filename = "Relative abundance of errors per Pos per threshold V14 SUP.jpg", 
  device = "jpg", 
  path = "out/plots/error_profile/threshold/", 
  width = 12, 
  height = 9, 
  dpi = 600)

n_errors_per_pos_categorized %>%
  filter(cluster_threshold <= 10 | cluster_threshold == 20 | cluster_threshold == 30) %>% 
  filter(kit == "V14" & algorithm == "HAC") %>% 
  mutate(
    facet_cols = paste(fragment, sample_type),
    facet_rows = factor(cluster_threshold, levels = seq(2, 30, 2)),
    color = paste(fragment, sample_type)
  ) %>% 
  ggplot(aes(position, n_errorrs_per_pos_categorized)) + 
  facet_grid(rows = vars(facet_rows), cols = vars(facet_cols), scale = "free_x") +
  geom_point(size = 0.2) +
  scale_x_continuous(breaks = seq(0,5000,1000)) +
  labs(
    y = "Number of errors", 
    x = "Position"
  )

ggsave(
  filename = "Relative abundance of errors per Pos per threshold V14 HAC.jpg", 
  device = "jpg", 
  path = "out/plots/error_profile/threshold/", 
  width = 12, 
  height = 9, 
  dpi = 600)

error_stats_merged %>%
  filter(str_detect(error_type, "-", negate = TRUE)) %>% 
  mutate(
    position = as.integer(str_remove(position, "\\.\\d+"))
  ) %>% 
  group_by(kit, algorithm, sample_type, fragment, reads_written, position) %>% 
  summarise(n_errors_per_pos = n()) %>% 
  mutate(
    n_errorrs_per_pos_categorized =
      cut(
        n_errors_per_pos, 
        breaks = c(-Inf, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, Inf), 
        labels = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, "10+"))
  ) %>% 
  filter(kit == "V14" & algorithm == "HAC") %>% 
  filter(reads_written <= 20) %>%
  mutate(
    facet_cols = paste(fragment, sample_type),
    facet_rows = factor(reads_written, levels = seq(2, 100, 2)),
    color = paste(fragment, sample_type)
  ) %>% 
  ggplot(aes(position, n_errorrs_per_pos_categorized)) + 
  facet_grid(rows = vars(facet_rows), cols = vars(facet_cols), scale = "free_x") +
  geom_point(size = 0.2) +
  scale_x_continuous(breaks = seq(0,5000,400)) +
  labs(
    y = "Number of errors", 
    x = "Position"
  )

ggsave(
  filename = "Relative abundance of errors per Pos per cluster size V14 HAC.jpg", 
  device = "jpg", 
  path = "out/plots/error_profile/size/", 
  width = 12, 
  height = 9, 
  dpi = 600)

error_stats_merged %>%
  filter(str_detect(error_type, "-", negate = TRUE)) %>% 
  mutate(
    position = as.integer(str_remove(position, "\\.\\d+"))
  ) %>% 
  group_by(kit, algorithm, sample_type, fragment, reads_written, position) %>% 
  summarise(n_errors_per_pos = n()) %>% 
  mutate(
    n_errorrs_per_pos_categorized =
      cut(
        n_errors_per_pos, 
        breaks = c(-Inf, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, Inf), 
        labels = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, "10+"))
  ) %>% 
  filter(kit == "V14" & algorithm == "SUP") %>% 
  filter(reads_written <= 20) %>%
  mutate(
    facet_cols = paste(fragment, sample_type),
    facet_rows = factor(reads_written, levels = seq(2, 100, 2)),
    color = paste(fragment, sample_type)
  ) %>% 
  ggplot(aes(position, n_errorrs_per_pos_categorized)) + 
  facet_grid(rows = vars(facet_rows), cols = vars(facet_cols), scale = "free_x") +
  geom_point(size = 0.2) +
  scale_x_continuous(breaks = seq(0,5000,400)) +
  labs(
    y = "Number of errors", 
    x = "Position"
  )

ggsave(
  filename = "Relative abundance of errors per Pos per cluster size V14 SUP.jpg", 
  device = "jpg", 
  path = "out/plots/error_profile/size/", 
  width = 12, 
  height = 9, 
  dpi = 600)

qscore_stats_merged %>% 
  mutate(
    facet_cols = paste(kit, algorithm),
    facet_rows = paste(fragment),
    color = paste(kit, algorithm)) %>% 
  ggplot(aes(qscore, medaka_quality, color = color)) + 
  facet_grid(cols = vars(facet_cols)) +
  scale_x_continuous(breaks = seq(0,70,5)) +
  scale_y_continuous(breaks = seq(0,70,5)) +
  coord_cartesian(xlim = c(0,70), ylim = c(0,70)) +
  theme_bw() +
  geom_abline(alpha = 0.2) +
  geom_point( position = "jitter", size = 0.4) +
  geom_smooth(method = "lm", color = "black", se = FALSE, show.legend = FALSE) +
  theme(
    legend.position = "bottom",
    text = element_text(size = 12)
  ) + 
  labs(
    x = "Q-score real",
    y = "Q-score Medaka", 
  )

ggsave(
  filename = "Correlation real vs Medaka Q-score.jpg", 
  device = "jpg", 
  path = "out/plots/Q-score/correlation_medaka/", 
  width = 12, 
  height = 9, 
  dpi = 600)

stats_summary %>%
  filter(!is.infinite(qscore)) %>%
  mutate(
    # qscore = ifelse(is.infinite(qscore), 50, qscore),
    facet_cols = paste(kit, algorithm),
    color = paste(kit, algorithm)) %>% 
  ggplot(aes(qscore, medaka_quality, color = color)) + 
  facet_grid(cols = vars(facet_cols)) +
  scale_x_continuous(breaks = seq(0,70,5)) +
  scale_y_continuous(breaks = seq(0,70,5)) +
  coord_cartesian(xlim = c(0,70), ylim = c(0,70)) +
  theme_bw() +
  geom_abline(alpha = 0.2) +
  geom_point( position = "jitter", size = 0.4) +
  geom_smooth(method = "lm", color = "black", se = FALSE, show.legend = FALSE) +
  theme(
    legend.position = "bottom",
    text = element_text(size = 12)
  ) + 
  labs(
    x = "Q-score real",
    y = "Q-score Medaka", 
  )

ggsave(
  filename = "Correlation real vs Medaka Q-score per cluster size.jpg", 
  device = "jpg", 
  path = "out/plots/Q-score/correlation_medaka/", 
  width = 12, 
  height = 9, 
  dpi = 600)

stats_summary_cluster_threshold %>%
  filter(!is.infinite(qscore)) %>%
  mutate(
    # qscore = ifelse(is.infinite(qscore), 50, qscore),
    facet_cols = paste(kit, algorithm),
    color = paste(kit, algorithm)) %>% 
  ggplot(aes(qscore, medaka_quality, color = color)) + 
  facet_grid(cols = vars(facet_cols)) +
  scale_x_continuous(breaks = seq(0,70,5)) +
  scale_y_continuous(breaks = seq(0,70,5)) +
  coord_cartesian(xlim = c(0,70), ylim = c(0,70)) +
  theme_bw() +
  geom_abline(alpha = 0.2) +
  geom_point( position = "jitter", size = 0.4) +
  geom_smooth(method = "lm", color = "black", se = FALSE, show.legend = FALSE) +
  theme(
    legend.position = "bottom",
    text = element_text(size = 12)
  ) + 
  labs(
    x = "Q-score real",
    y = "Q-score Medaka", 
  )

ggsave(
  filename = "Correlation real vs Medaka Q-score per cluster threshold.jpg", 
  device = "jpg", 
  path = "out/plots/Q-score/correlation_medaka/", 
  width = 12, 
  height = 9, 
  dpi = 600)


