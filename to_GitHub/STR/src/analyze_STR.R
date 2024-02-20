library(tidyverse)
library(Biostrings)
library(jsonlite)
library(ggridges)
library(viridis)
library(ggh4x)

read_extacted_STR <- function(path){
  barcode <- str_extract(path, "barcode\\d{2}")
  kit <- str_extract(path, "V14|R9")
  algorithm <- str_extract(path, "SUPDUP|SUP|HAC")
  
  extracted_STR <- 
    read_tsv(path) %>% 
    mutate(
      barcode = barcode, 
      kit = kit, 
      algorithm = algorithm
    )
}
read_sample_sheet <- function(path){
  sample_sheet <- 
    as_tibble(read_json(sample_sheet_5104_path, simplifyVector = TRUE)) %>% 
    drop_na() %>% 
    mutate(barcode = str_replace(Barcode, "NB", "barcode"), 
           sample = str_sub(Sample, end = -6),
           sample = ifelse(str_detect(sample, "SAPHIR"), sample, ifelse(sample == "A_B_1000_0", "100% KIV-2A", "100% KIV-2B")),
           fragment = as.character(5104)
    ) %>% 
    select(!c(Barcode, Sample)) 
}
get_annotated_reference_data <- function(reference_data){
  type_b <- reference_data %>% 
    filter(pos == 666 | pos == 621 | pos == 594) %>% 
    group_by(id) %>% 
    summarize(n = n()) %>% 
    filter(n == 3) %>% 
    mutate(type_B_ngs = TRUE) %>% 
    select(id, type_B_ngs)
  
  type_b_c <- reference_data %>% 
    filter(pos == 666 | pos == 621 | pos == 594) %>%
    select(id, pos, variant_level) %>% 
    pivot_wider(names_from = pos, values_from = variant_level) %>%
    mutate(
      is_type_c = (as.numeric(`594`) - as.numeric(`666`) > 0.01)
    ) 
  
  not_type_b <- reference_data %>%
    anti_join(type_b) %>%
    group_by(id) %>% 
    summarize() %>% 
    mutate(type_B_ngs = FALSE)
  
  type_b_annotated <- type_b %>% 
    bind_rows(not_type_b) %>% 
    mutate(
      sample = ifelse(str_detect(id, "AK"), str_extract(id, "AK\\d\\d"), id))
  
  reference_data_parsed <- reference_data %>%
    filter(!is.na(id)) %>% 
    group_by(id) %>%
    mutate(
      sample = ifelse(str_detect(id, "AK"), str_extract(id, "AK\\d\\d"), id),
      expected_KIV2_repeats = ceiling(ifelse(heterozygous == 1, as.numeric(isosum) - 18, kiv2_qpcr_total)),
      expected_KIV2_repeats_source = ifelse(heterozygous == 1, "WB", "qPCR")
    ) %>% 
    group_by(sample, heterozygous, lmw, expected_KIV2_repeats, expected_KIV2_repeats_source, iso1, iso2, isosum, kiv2_qpcr_total, lpa) %>%
    summarize() %>% 
    inner_join(type_b_annotated, by = "sample") %>% 
    select(!id)
  
  return(reference_data_parsed)
}
read_ddPCR_sheet <- function(path){
  experiment <- str_extract(path, "first|repeat")
  replicates <- read_tsv(path) %>% 
    filter(str_detect(Well, "M", negate = TRUE)) %>%
    filter(Sample != "NTC") %>% 
    group_by(Sample) %>% 
    summarise(replicates = n())
  
  read_tsv(path) %>% 
    filter(str_detect(Well, "M")) %>% 
    mutate(
      experiment = experiment,
      confidence_interval = TotalCNVMax - TotalCNVMin, 
      sample = paste0("SAPHIR_", Sample)
    ) %>% 
    drop_na() %>% 
    inner_join(replicates) %>% 
    mutate(
      expected_KIV2_repeats_ddPCR = CNV - 2,
      TotalCNVMax = TotalCNVMax - 2, 
      TotalCNVMin = TotalCNVMin - 2
    )
}

extracted_STRs_path <- "20231121_STR/out/V14/SUP/barcode03/extracted_STR.tsv"
extracted_STRs_paths <- list.files(
  path = "20231121_STR/out/", 
  pattern = "extracted_STR.tsv", 
  recursive = TRUE, 
  full.names = TRUE
)

ddPCR_reference_paths <- 
  list.files(
    path = "20231121_STR/data/ddPCR/", 
    pattern = "ddPCR", 
    recursive = TRUE, 
    full.names = TRUE
  )

reference_data_path <- "20231121_STR/data/20221122_201803_pop_variants_noBAQ.csv"

sample_sheet_5104_path <- "20231121_STR/data/run_SAPHIR_5104/lib/Barcode_Sample_overview.js"

ddPCR_reference <- lapply(ddPCR_reference_paths, read_ddPCR_sheet) %>% 
  bind_rows() %>% 
  filter( replicates == 3 ) %>% 
  filter(!(sample == "SAPHIR_5254" & experiment == "first")) %>% 
  add_row(
    sample = "100% KIV-2A", 
    expected_KIV2_repeats_ddPCR = 1 
  ) %>% 
  add_row(
    sample = "100% KIV-2B", 
    expected_KIV2_repeats_ddPCR = 1
  ) %>% 
  select(sample, expected_KIV2_repeats_ddPCR)

reference_data <- read_csv(reference_data_path, na = c("", "NA", "#N/A")) %>% 
  filter(str_detect(id, "SAPHIR"))

reference_data_parsed <- 
  get_annotated_reference_data(reference_data) %>% 
  ungroup() %>%
  select(sample, type_B_ngs, lmw, lpa) %>%
  add_row(
    sample = "100% KIV-2A",
    type_B_ngs = FALSE, 
    lmw = NA, 
    lpa = NA
  ) %>%
  add_row(
    sample = "100% KIV-2B",
    type_B_ngs = TRUE,
    lmw = NA, 
    lpa = NA
  ) %>% 
  inner_join(ddPCR_reference)

extracted_STRs <- lapply(extracted_STRs_paths, read_extacted_STR) %>% 
  bind_rows()

sample_sheet <- 
  read_sample_sheet(sample_sheet_5104_path)

extracted_STRs_parsed <- 
  extracted_STRs %>% 
  inner_join(sample_sheet) %>% 
  rowwise() %>% 
  mutate(
    CA_occurence = length(str_extract_all(STR_sequence, "CA", simplify = TRUE)), 
    STR_length = nchar(STR_sequence) / 2,
    STR_end = str_locate_all(STR_sequence, "CA")[[1]][length(str_locate_all(STR_sequence, "CA")[[1]])],
    STR_start = str_locate(STR_sequence, "CA")[1],
    STR_sequence_trimmed = str_sub(STR_sequence, start = STR_start, end = STR_end), 
    STR_length_trimmed = nchar(STR_sequence_trimmed) / 2,
    mutations = str_remove_all(STR_sequence_trimmed, "CA"),
    mutation_position = ifelse(mutations == "", 0, str_locate(STR_sequence_trimmed, mutations)[1])
  ) %>% 
  ungroup() %>%
  reframe(
    sample = sample, 
    fragment = fragment,
    kit = kit, 
    algorithm = algorithm, 
    STR_sequence = STR_sequence,
    STR_length_trimmed = STR_length_trimmed,
    STR_sequence_trimmed = STR_sequence_trimmed,
    mutations = mutations, 
    mutation_position = mutation_position,
    CA_occurence = CA_occurence,
    occurence = n(),
    .by = c(STR_length_trimmed, sample, fragment, kit, algorithm), 
  ) %>% 
  reframe(
    sample = sample, 
    fragment = fragment,
    kit = kit, 
    algorithm = algorithm, 
    STR_sequence = STR_sequence,
    STR_sequence_trimmed = STR_sequence_trimmed, 
    STR_length_trimmed = STR_length_trimmed,
    mean_STR_length_trimmed = mean(STR_length_trimmed), 
    median_STR_length_trimmed = median(STR_length_trimmed), 
    quartile_1_STR_length_trimmed = quantile(STR_length_trimmed, 1/4),
    quartile_3_STR_length_trimmed = quantile(STR_length_trimmed, 3/4),
    mutations = mutations, 
    mutation_position = mutation_position,
    CA_occurence = CA_occurence,
    occurence = occurence,
    total_sequences = n(),
    .by = c(sample, fragment, kit, algorithm), 
  ) %>% 
  inner_join(reference_data_parsed) %>%
  filter( occurence > total_sequences * 0.0085)

degenerated_STR_per_sample <- 
  extracted_STRs_parsed %>%
  unique() %>% 
  select(sample, mutations, mutation_position, lpa, lmw, type_B_ngs) %>% 
  reframe(
    sequences = n(), 
    .by = c(sample, mutations, mutation_position)
  )

degenerated_STR_summary <- 
  degenerated_STR_per_sample %>% 
  reframe(
    samples = n(), 
    sequences = sum(sequences), 
    .by = c(mutations, mutation_position)
  )
  
diversity_STR <-
  extracted_STRs_parsed %>% 
  unique() %>% 
  group_by(STR_sequence_trimmed, STR_length_trimmed) %>% 
  summarise(
    occurence = n()
    )

total_STR_variants <- 
  diversity_STR %>% 
  nrow()

min_STR_length_trimmed <- min(extracted_STRs_parsed$STR_length_trimmed )
max_STR_length_trimmed = max(extracted_STRs_parsed$STR_length_trimmed)

p.STR_length_trimmed_per_sample_ridgeline_barplot <-
extracted_STRs_parsed %>% 
  ggplot(aes(
    x = STR_length_trimmed, 
    y = reorder(sample, mean_STR_length_trimmed),
    # y = reorder(sample, lpa),
    fill = sample
    )) +
  # facet_wrap(vars(type_B_ngs)) +
  geom_density_ridges(
    alpha = 0.9, 
    stat = "binline", 
    binwidth = 1, 
    rel_min_height = 0.001, 
    scale = 0.99
    ) +
  scale_x_continuous(breaks = seq(0, 40, 1), expand = expansion(add = c(0.6, 0.6))) +
  scale_y_discrete(expand = expansion( add = c(0.1, 1))) +
  coord_cartesian(xlim = c(min_STR_length_trimmed, max_STR_length_trimmed)) +
  labs(
    x = "STR length [bp]",
    y = ""
    ) +
  theme_bw() +
  theme(
    legend.position="none",
    panel.spacing = unit(0.1, "lines"),
    strip.text.x = element_text(size = 14), 
    panel.grid.minor.x = element_blank(),
    text = element_text(size = 14)
  )

p.STR_length_trimmed_per_sample_ridgeline_barplot

ggsave(
  plot = p.STR_length_trimmed_per_sample_ridgeline_barplot,
  filename = "STR_length_trimmed_per_sample_ridgeline_barplot.jpeg", 
  path = "20231121_STR/out/V14/SUP/plots/", 
  device = "jpeg", 
  dpi = 600, 
  width = 12, 
  height = 9
)

p.STR_length_trimmed_per_sample_relative_barplots <-
extracted_STRs_parsed %>%
  mutate(
    relative_occurence = occurence / total_sequences
  ) %>% 
  ggplot(aes(
    x = STR_length_trimmed,
    y = relative_occurence,
    fill = sample, 
  )) +
  facet_wrap(vars(sample), 
             nrow = 4,
             # ncol = 5
             ) +
  geom_bar(width = 1, stat = "identity", position = "identity", color = "black", linewidth = 0.2) +
  scale_x_continuous(breaks = seq(0, 40, 2)) +
  scale_y_continuous(breaks = seq(0, 1, 0.1)) +
  labs(
    x = "STR length [bp]",
    y = "Relative occurence"
  ) +
  theme_bw() +
  theme(
    text = element_text(size = 12),
    axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1, size = 12, margin = margin(t = 1, r = 2, b = 2, l = 1, "pt")),
    axis.text.y = element_text(size = 12, margin = margin(t = 1, r = 2, b = 2, l = 1, "pt")),
    axis.title.x = element_text(size = 14, margin = margin(t = 1, r = 2, b = 2, l = 1, "pt")),
    axis.title.y = element_text(angle = 90, size = 14, margin = margin(t = 1, r = 2, b = 2, l = 1, "pt")),
    panel.background = element_rect(fill = "white", color = NA),
    panel.border = element_rect(fill = NA, color = "grey20"),
    panel.grid.minor = element_line(linewidth = rel(0.5), color = "grey92"),
    panel.grid = element_line(color = "grey92"),
    strip.text = element_text(face = "bold", size = 12, margin = margin(t = 1, r = 2, b = 2, l = 1, "pt")),
    strip.background = element_rect(fill="white", color = "grey20"),
    legend.key = element_rect(fill = "white", color = NA),
    legend.key.width = unit(4, 'pt'),
    legend.margin = margin(0, 0, 0, 0),
    plot.margin = margin(0, 0, 0, 0),
    plot.title = element_blank(),
    legend.position="none", 
    panel.grid.minor.x = element_blank()
  )

p.STR_length_trimmed_per_sample_relative_barplots

ggsave(
  plot = p.STR_length_trimmed_per_sample_relative_barplots,
  filename = "STR_length_trimmed_per_sample_relative_barplots.jpeg", 
  path = "20231121_STR/out/V14/SUP/plots/", 
  device = "jpeg", 
  dpi = 600, 
  width = 12, 
  height = 9
)