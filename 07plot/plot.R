#!/usr/bin/env Rscript
library(ggplot2)
library(tidyr)
library(dplyr)
library(stringr)
library(colorspace)
library(cowplot)

# ---------------------------------------------------------------------------
# 1. Load Data
# ---------------------------------------------------------------------------
# Load the data
# E:\paper\07 costatae evolution\data_analysis\alleles_geographic_distribution\01total
raw_data_original <- read.table("merged_loci_with_genes.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE, quote = "'", check.names = FALSE)

# Check if there are duplicate arab_gene entries (multiple loci mapping to same arab_gene)
# For duplicate arab_genes, average the frequency values across all populations
cat("Original data has", nrow(raw_data_original), "rows\n")

# Count rows with NA annotation for reporting
na_annotation_count <- sum(is.na(raw_data_original$annotation) |
                           raw_data_original$annotation == "NA" |
                           raw_data_original$annotation == "")
cat("Rows with NA annotation (will be ignored):", na_annotation_count, "\n")

# Identify rows with valid arab_gene AND valid annotation (not NA)
# Filter out rows where annotation is NA - these should be ignored
valid_arab_gene <- raw_data_original %>%
  filter(!is.na(arab_gene) & arab_gene != "NA" & arab_gene != "") %>%
  filter(!is.na(annotation) & annotation != "NA" & annotation != "")

cat("Rows with valid arab_gene AND annotation:", nrow(valid_arab_gene), "\n")

# Get population columns (numeric frequency data)
population_cols <- c("LJS", "LP", "LJA", "JS", "WZ", "DQ", "DCX", "YJA", "JLX", "JD", "SLS", "MEK", "LFG", "EBY", "SDX", "BX", "MYL", "RTX", "BMX", "TZZ", "DTX", "BYS", "XYB", "SND", "TTH", "QLL", "TBX", "TBS", "HP", "HHG", "LYC", "LHS", "SWP", "LS", "PQG", "BA", "XLA", "LYL", "LLZ", "FHSJ", "CBS", "BJF")

# For each arab_gene, average the frequencies across all loci
# Keep the first annotation (name and GO annotation) for each arab_gene
aggregated_data <- valid_arab_gene %>%
  group_by(arab_gene) %>%
  summarise(
    name = first(name),
    annotation = first(annotation),
    across(all_of(population_cols), ~ mean(.x, na.rm = TRUE)),
    n_loci = n()  # Track how many loci were averaged
  ) %>%
  ungroup()

# Report genes with multiple loci
multi_loci_genes <- aggregated_data %>% filter(n_loci > 1)
if (nrow(multi_loci_genes) > 0) {
  cat("\nFound", nrow(multi_loci_genes), "arab_genes with multiple loci (averaged):\n")
  print(multi_loci_genes %>% select(arab_gene, n_loci, name))
}

# Prepare final data structure for downstream analysis
raw_data <- aggregated_data %>%
  select(-n_loci) %>%
  rename(Gene_Symbol = arab_gene, GO_Keywords = annotation)

cat("\nAfter aggregation:", nrow(raw_data), "unique genes\n")

# ---------------------------------------------------------------------------
# 2. Define GO Terms, Categories, Colors, and Orders
# ---------------------------------------------------------------------------

go_terms_list <- c(
  "root development", "flower development", "reproductive phase transition",
  "vernalization response", # Development & Reproduction
  "DNA repair", "DNA damage response", # DNA Processes
  "response to abscisic acid", "response to auxin", "response to salicylic acid",
  "response to jasmonic acid", "response to ethylene", "response to cytokinin", "response to gibberellin",# Hormone Signaling
  "response to bacterium", "response to fungus", "response to virus",
  "response to oomycetes", "response to nematode", # Biotic Stress
  "response to cold", "response to heat", "response to wounding",
  "response to water deprivation", "response to salt stress", "response to UV", # Abiotic Stress
  "response to absence of light","circadian rhythm", "photomorphogenesis",
  "photoperiodism, flowering"# Light & Circadian
)


go_term_definitions <- data.frame(
  GO_term = go_terms_list,
  Category = factor(c(
    rep("Development & Reproduction", 4),
    rep("DNA Processes", 2),
    rep("Hormone Signaling", 7),
    rep("Biotic Stress", 5),
    rep("Abiotic Stress", 6),
    rep("Light & Circadian", 4)
  ), levels = c("Development & Reproduction", "DNA Processes", "Hormone Signaling",
                "Biotic Stress", "Abiotic Stress", "Light & Circadian")),
  stringsAsFactors = FALSE
)

category_colors <- c(
  "Development & Reproduction" = "#FAB9AC",
  "DNA Processes" = "#7BBC53",            
  "Hormone Signaling" = "#DE6736",
  "Biotic Stress" = "#67C1EC",            
  "Abiotic Stress" = "#dcbeff",      
  "Light & Circadian" = "#E6B90D"         
)
go_term_definitions$Color <- category_colors[go_term_definitions$Category]

level_order_go_user <- go_terms_list 

# Gene order for the plots
level_order_gene <- unique(raw_data$Gene_Symbol)

# Population columns & groups for population
all_population_columns <- c("LJS", "LP", "LJA", "JS", "WZ", "DQ", "DCX", "YJA", "JLX", "JD", "SLS", "MEK", "LFG", "EBY", "SDX", "BX", "MYL", "RTX", "BMX", "TZZ", "DTX", "BYS", "XYB", "SND", "TTH", "QLL", "TBX", "TBS", "HP", "HHG", "LYC", "LHS", "SWP", "LS", "PQG", "BA", "XLA", "LYL", "LLZ", "FHSJ", "CBS", "BJF")
population_groups <- list(
  utilis = c("LJS", "LJA", "DCX", "JS", "WZ", "DQ", "LP"),
  uti_albo = c("RTX", "MYL", "BMX", "MEK", "SLS", "BX", "LFG", "JLX", "EBY", "YJA", "SDX","JD"),
  albo = c("SND", "XYB", "TBX", "DTX", "HP", "HHG", "TBS", "TZZ", "QLL", "LYC", "SWP", "BYS", "TTH", "LHS"),
  albo_erman = c("LS" , "PQG", "BA", "XLA"),
  ermanii = c("LYL", "LLZ", "FHSJ", "CBS", "BJF")
)
level_order_location <- rev(names(population_groups)) # For bottom-to-top display (e.g., utilis at top)


# ---------------------------------------------------------------------------
# 3. Data Preparation
# ---------------------------------------------------------------------------

# GO Annotation
go_plot_data_raw <- raw_data %>%
  select(Gene_Symbol, GO_Keywords) %>%
  mutate(GO_Keywords = str_split(GO_Keywords, ";\\s*")) %>%
  unnest(GO_Keywords) %>%
  rename(Gene = Gene_Symbol, GO_term = GO_Keywords)

# unique genes and GO terms
genes_before_filter <- unique(go_plot_data_raw$Gene)
go_terms_before_filter <- unique(go_plot_data_raw$GO_term)

go_plot_data <- go_plot_data_raw %>%
  filter(Gene %in% level_order_gene, GO_term %in% level_order_go_user) %>%
  distinct(Gene, GO_term) %>%
  mutate(Flag = 1)

genes_after_filter <- unique(go_plot_data$Gene)
go_terms_after_filter <- unique(go_plot_data$GO_term)

# Identify filtered out genes &  GO terms
filtered_out_genes <- setdiff(genes_before_filter, genes_after_filter)
if (length(filtered_out_genes) > 0) {
  warning(paste("Gene not in 'level_order_gene': ",
                paste(filtered_out_genes, collapse = ", ")))
}
filtered_out_go_terms <- setdiff(go_terms_before_filter, go_terms_after_filter)

if (length(filtered_out_go_terms) > 0) {
  warning(paste("GO not in 'go_terms_list': ",
                paste(filtered_out_go_terms, collapse = ", ")))
}


go_plot_data <- go_plot_data %>%
  left_join(go_term_definitions %>% select(GO_term, Color), by = "GO_term") %>%
  rename(mycolor = Color)
go_plot_data$mycolor[is.na(go_plot_data$mycolor)] <- "grey50" 

go_plot_data$mydotsize <- 1.7

# --- B. For Allele Fraction Plot ---
fraction_data_individual <- raw_data %>%
  select(Gene_Symbol, all_of(all_population_columns)) %>%
  pivot_longer(cols = all_of(all_population_columns), names_to = "Location_Individual", values_to = "Fraction_Individual") %>%
  rename(Gene = Gene_Symbol)

#  calculate the average fraction for each group
fraction_plot_data <- fraction_data_individual %>%
  mutate(
    Location = case_when(
      Location_Individual %in% population_groups$utilis ~ "utilis",
      Location_Individual %in% population_groups$uti_albo ~ "uti_albo",
      Location_Individual %in% population_groups$albo ~ "albo",
	  Location_Individual %in% population_groups$albo_erman ~ "albo_erman",
      Location_Individual %in% population_groups$ermanii ~ "ermanii",
      TRUE ~ NA_character_ 
    )
  ) %>%
  filter(!is.na(Location)) %>% 
  group_by(Gene, Location) %>%
  summarise(Fraction = mean(Fraction_Individual, na.rm = TRUE), .groups = 'drop') %>%
  filter(Gene %in% level_order_gene, Location %in% names(population_groups)) # Ensure only desired genes and groups are kept


fraction_plot_data <- fraction_plot_data %>%
  group_by(Gene) %>%
  filter(any(Fraction >= 0.3)) %>% # filter the average fraction for each group
  ungroup()
level_order_gene <- fraction_plot_data %>%
  group_by(Gene) %>%
  summarise(mean_fraction = mean(Fraction)) %>%
  arrange(desc(mean_fraction)) %>%
  pull(Gene)

go_plot_data <- go_plot_data %>%
  filter(Gene %in% level_order_gene)

base_dot_size_fraction <- 1.7 # Point size
fraction_plot_data$mydotsize <- base_dot_size_fraction * fraction_plot_data$Fraction

fraction_colors <- c("#e6194B", "#ffe119", "#3cb44b")
fraction_color_palette_func <- colorRampPalette(fraction_colors, space = "Lab")
fraction_plot_data$aux <- round(fraction_plot_data$Fraction * 99) + 1
fraction_plot_data$myplotcolors <- fraction_color_palette_func(100)[fraction_plot_data$aux]


# ---------------------------------------------------------------------------
# 4. Plotting
# ---------------------------------------------------------------------------

GOplotFUN_user <- function(data_plot_Ngo, go_order, gene_order, plot_margins, plot_title_text) {
  data_plot_Ngo$GO_term <- factor(data_plot_Ngo$GO_term, levels = go_order)
  data_plot_Ngo$Gene <- factor(data_plot_Ngo$Gene, levels = gene_order)
  data_plot_Ngo <- data_plot_Ngo[!is.na(data_plot_Ngo$GO_term) & !is.na(data_plot_Ngo$Gene), ]
  
  plotN1go <- ggplot(data_plot_Ngo, aes(x = GO_term, y = Gene)) +
    geom_point(stat = 'identity', size = data_plot_Ngo$mydotsize, colour = data_plot_Ngo$mycolor, na.rm = TRUE) + 
    ggtitle(plot_title_text) +
    theme_minimal(base_size = 10) +
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text.x = element_text(size = 5.5, angle = 90, hjust = 1, vjust = 0.5,color="#000000"),
      axis.text.y = element_text(size = 5.5, angle = 0,color="#000000"),
      panel.grid.major = element_line(linewidth = 0.3, linetype = 'solid', colour = "#D1D1D1"),
      panel.background = element_rect(fill = 'white', colour = 'black'),
      plot.margin = unit(plot_margins, 'cm'),
      plot.title = element_text(hjust = 0.01, size = 12)
    ) +
    coord_flip()
  return(plotN1go)
}

plotDots_user <- function(data_plot_dots, location_order, gene_order, plot_margins, base_dot_sz) {
  data_plot_dots$Location <- factor(data_plot_dots$Location, levels = location_order)
  data_plot_dots$Gene <- factor(data_plot_dots$Gene, levels = gene_order)
  data_plot_dots <- data_plot_dots[!is.na(data_plot_dots$Location) & !is.na(data_plot_dots$Gene), ]
  
  plotN1_dots <- ggplot(data_plot_dots, aes(x = Location, y = Gene)) +
    geom_point(stat = 'identity', size = data_plot_dots$mydotsize, colour = data_plot_dots$myplotcolors, na.rm = TRUE) + 
    theme_minimal(base_size = 10) +
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text.x = element_text(size = 5.5,angle = 90, hjust = 1, vjust = 0.5,color="#000000"),
      axis.text.y = element_text(size = 5.5,color="#000000"),
      panel.grid.major = element_line(linewidth = 0.3, linetype = 'solid', colour = "#D1D1D1"),
      panel.background = element_rect(fill = 'white', colour = 'black'),
      plot.margin = unit(plot_margins, 'cm')
    ) +
    coord_flip()
  return(plotN1_dots)
}

# ---------------------------------------------------------------------------
# 5. Generate and Combine Plots
# ---------------------------------------------------------------------------
plot_title_main <- "B. ashburneri origin"
go_plot_margins <- c(0.5, 0.5, 0.1, 0.5) # Increased left margin for GO terms if needed

final_go_order_for_plot <- rev(level_order_go_user)

go_plot_data$Gene <- factor(go_plot_data$Gene, levels = level_order_gene)
go_plot_data$GO_term <- factor(go_plot_data$GO_term, levels = final_go_order_for_plot) # Use reversed order for plotting

final_go_plot <- GOplotFUN_user(go_plot_data,
                                go_order = final_go_order_for_plot, # Use reversed
                                gene_order = level_order_gene,
                                plot_margins = go_plot_margins,
                                plot_title_text = plot_title_main)

fraction_plot_margins <- c(0.1, 0.5, 0.5, 0.5)
fraction_plot_data$Gene <- factor(fraction_plot_data$Gene, levels = level_order_gene)
fraction_plot_data$Location <- factor(fraction_plot_data$Location, levels = level_order_location) 

final_fraction_plot <- plotDots_user(fraction_plot_data,
                                     location_order = level_order_location,
                                     gene_order = level_order_gene,
                                     plot_margins = fraction_plot_margins,
                                     base_dot_sz = base_dot_size_fraction)

combined_plot <- plot_grid(
  final_go_plot,
  final_fraction_plot,
  ncol = 1,
  align = 'v',
  axis = "lr",
  rel_heights = c(0.75, 0.25) # GO plot takes more space
)

#combined_plot


#ggsave("my_custom_GO_fraction_plot.pdf", combined_plot, width = 300, height = 150, units = "mm")


# ---------------------------------------------------------------------------
# 6. Calculate Category Frequencies
# ---------------------------------------------------------------------------
# Create a mapping from each gene to its functional categories
gene_to_category_map <- raw_data %>%
  select(Gene_Symbol, GO_Keywords) %>%
  mutate(GO_term = strsplit(gsub("'", "", GO_Keywords), ";\\s*")) %>%
  unnest(GO_term) %>%
  left_join(go_term_definitions, by = "GO_term") %>%
  # Keep only the essential columns and remove duplicates
  # (e.g., if a gene has two "Abiotic Stress" terms, we only need one link)
  select(Gene = Gene_Symbol, Category) %>%
  filter(!is.na(Category)) %>%
  distinct(Gene, Category)

# Calculate the average allele fraction for each gene in each population group
# This is similar to code for 'fraction_plot_data'
group_allele_fractions <- raw_data %>%
  select(Gene_Symbol, all_of(all_population_columns)) %>%
  pivot_longer(
    cols = -Gene_Symbol,
    names_to = "Location_Individual",
    values_to = "Fraction"
  ) %>%
  mutate(
    Group = case_when(
      Location_Individual %in% population_groups$utilis     ~ "utilis",
      Location_Individual %in% population_groups$uti_albo   ~ "uti_albo",
      Location_Individual %in% population_groups$albo       ~ "albo",
      Location_Individual %in% population_groups$albo_erman ~ "albo_erman",
      Location_Individual %in% population_groups$ermanii    ~ "ermanii"
    )
  ) %>%
  group_by(Gene = Gene_Symbol, Group) %>%
  summarise(Mean_Gene_Fraction = mean(Fraction, na.rm = TRUE), .groups = 'drop')


# Combine and calculate the final category frequency
category_frequency_by_group <- group_allele_fractions %>%
  left_join(gene_to_category_map, by = "Gene") %>%
  filter(!is.na(Category)) %>%
  group_by(Group, Category) %>%
  summarise(Mean_Category_Fraction = mean(Mean_Gene_Fraction, na.rm = TRUE), .groups = 'drop')

# ---------------------------------------------------------------------------
# 3. View the Result
# ---------------------------------------------------------------------------

print(category_frequency_by_group)
category_frequency_wide <- category_frequency_by_group %>%
  pivot_wider(names_from = Group, values_from = Mean_Category_Fraction)

print(category_frequency_wide)

## Create a heatmap to visualize the results
category_frequency_by_group$Group <- factor(category_frequency_by_group$Group, levels = names(population_groups))
category_frequency_by_group$Category <- factor(category_frequency_by_group$Category, levels = levels(go_term_definitions$Category))

heatmap_plot <- ggplot(category_frequency_by_group, aes(x = Group, y = reorder(Category, desc(Category)), fill = Mean_Category_Fraction)) +
  geom_tile(color = "white", lwd = 1.5) +
  geom_text(aes(label = round(Mean_Category_Fraction, 2)), color = "black", size = 3) +
  scale_fill_gradient(low = "#e6f5ff", high = "#3869c1", name = "Mean Allele\nFraction") +
  theme_minimal(base_size = 12) +
  labs(
    x = "Population Group",
    y = "Functional Category"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  )


#low = "#ffe3dc", high = "#ea4d1c"

#low = "#fff6e9", high = "#e89425"


print(heatmap_plot)

######## pairing frequencies for five pop

library(ggpubr)
library(tidyr)
library(dplyr)
library(stringr)
library(tibble)
library(grid)

# NOTE: Update this path to your local environment
# setwd("E://paper/07 costatae evolution/data_analysis/polyAncestor20250329/")
# For Linux/Mac, use forward slashes, e.g.:
# setwd("/path/to/your/polyAncestor20250329/")

cols_keep <- c("pendula","platyphylla","populifolia","nana","occidentalis",
"cos","ash","cos_ash","bug","mca","len","basal")

species_levels <- c("len","nana","occidentalis","populifolia",
"pendula","platyphylla","mca","bug","cos","ash","cos_ash","basal")

species_labels <- c("B. lenta","B. nana","B. occidentalis","B. populifolia",
"B. pendula","B. platyphylla","B. mcallisteri","B. buggsii",
"B. costata","B. ashburneri","cos_ash","basal")

title_from_file <- function(file) {
  f <- tolower(basename(file))
  if (str_detect(f, "ash"))   return("ash ref")
  if (str_detect(f, "bugg"))  return("bugg ref")
  if (str_detect(f, "cos"))   return("cos ref")
  if (str_detect(f, "lenta")) return("lenta ref")
  "reference"
}

prep_plot_for_group <- function(file, prefixes, group_label, ylim_max = 40) {
  df <- read.delim(file, header = TRUE, row.names = 1, check.names = FALSE)
  cols_present <- intersect(cols_keep, colnames(df))
  
  if (length(cols_present) == 0) {
    stop(sprintf(
      "Required columns not found in %s. Expected any of: %s",
      file, paste(cols_keep, collapse = ", ")
    ))
  }
  
  df <- df[, cols_present, drop = FALSE]
  re <- paste0("^(", paste(prefixes, collapse = "|"), ")")
  df_sub <- df[grepl(re, rownames(df)), , drop = FALSE]
  
  if (nrow(df_sub) == 0) {
    stop(sprintf(
      "No samples matched the given prefixes in %s. Prefixes: %s",
      file, paste(prefixes, collapse = ", ")
    ))
  }
  
  df_sub <- as.data.frame(df_sub * 100) %>% rownames_to_column(var = "sample")
  
  df_long <- df_sub %>%
    pivot_longer(
      cols = -sample,
      names_to = "species",
      values_to = "value"
    ) %>%
    mutate(
      species = factor(species, levels = species_levels, labels = species_labels)
    ) %>%
    filter(!is.na(species))
  
  ggbarplot(
    df_long, x = "species", y = "value",
    add = "mean_se",
    xlab = "", ylab = sprintf("pairing with %s (%%)", group_label),
    ylim = c(0, ylim_max)
  ) +
    stat_summary(
      aes(label = round(after_stat(y), 1)),
      geom = "text", fun = mean, vjust = -0.7, size = 3
    ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "italic")) +
    ggtitle(title_from_file(file))
}

file_ash   <- "pairing_matrix_all_species_ash.txt"
file_bugg  <- "pairing_matrix_all_species_bugg.txt"
file_cos   <- "pairing_matrix_all_species_cos.txt"
file_lenta <- "pairing_matrix_all_species_lenta.txt"

prefix_group1 <- c("LJS", "LJA", "DCX", "JS", "WZ", "DQ", "LP")
prefix_group2 <- c("RTX", "MYL", "BMX", "MEK", "SLS", "BX", "LFG", "JLX", "EBY", "YJA", "SDX","JD")
prefix_group3 <- c("SND", "XYB", "TBX", "DTX", "HP", "HHG", "TBS", "TZZ", "QLL", "LYC", "SWP", "BYS", "TTH", "LHS")
prefix_group4 <- c("LS" , "PQG", "BA", "XLA")
prefix_group5 <- c("LYL", "LLZ", "FHSJ", "CBS", "BJF")

label_group1 <- "pure utilis"
label_group2 <- "contact zone between uti&albo"
label_group3 <- "pure albosinensis"
label_group4 <- "contact zone between albo&erman"
label_group5 <- "pure ermanii"

## ===== Group 1: pure utilis=====
g1_ash   <- prep_plot_for_group(file_ash,   prefixes = prefix_group1, group_label = label_group1)
g1_bugg  <- prep_plot_for_group(file_bugg,  prefixes = prefix_group1, group_label = label_group1)
g1_cos   <- prep_plot_for_group(file_cos,   prefixes = prefix_group1, group_label = label_group1)
g1_lenta <- prep_plot_for_group(file_lenta, prefixes = prefix_group1, group_label = label_group1)

fig_group1 <- ggarrange(
  g1_ash, g1_bugg, g1_cos, g1_lenta,
  nrow = 1, ncol = 4, common.legend = TRUE, legend = "bottom", align = "hv"
)

## ===== Group 2: contact zone between utilis and albosinese=====
g2_ash   <- prep_plot_for_group(file_ash,   prefixes = prefix_group2, group_label = label_group2)
g2_bugg  <- prep_plot_for_group(file_bugg,  prefixes = prefix_group2, group_label = label_group2)
g2_cos   <- prep_plot_for_group(file_cos,   prefixes = prefix_group2, group_label = label_group2)
g2_lenta <- prep_plot_for_group(file_lenta, prefixes = prefix_group2, group_label = label_group2)

fig_group2 <- ggarrange(
  g2_ash, g2_bugg, g2_cos, g2_lenta,
  nrow = 1, ncol = 4, common.legend = TRUE, legend = "bottom", align = "hv"
)

## ===== Group 3: pure albosinese=====
g3_ash   <- prep_plot_for_group(file_ash,   prefixes = prefix_group3, group_label = label_group3)
g3_bugg  <- prep_plot_for_group(file_bugg,  prefixes = prefix_group3, group_label = label_group3)
g3_cos   <- prep_plot_for_group(file_cos,   prefixes = prefix_group3, group_label = label_group3)
g3_lenta <- prep_plot_for_group(file_lenta, prefixes = prefix_group3, group_label = label_group3)

fig_group3 <- ggarrange(
  g3_ash, g3_bugg, g3_cos, g3_lenta,
  nrow = 1, ncol = 4, common.legend = TRUE, legend = "bottom", align = "hv"
)

## ===== Group 4:contact zone between albosinese and ermanii=====
g4_ash   <- prep_plot_for_group(file_ash,   prefixes = prefix_group4, group_label = label_group4)
g4_bugg  <- prep_plot_for_group(file_bugg,  prefixes = prefix_group4, group_label = label_group4)
g4_cos   <- prep_plot_for_group(file_cos,   prefixes = prefix_group4, group_label = label_group4)
g4_lenta <- prep_plot_for_group(file_lenta, prefixes = prefix_group4, group_label = label_group4)

fig_group4 <- ggarrange(
  g4_ash, g4_bugg, g4_cos, g4_lenta,
  nrow = 1, ncol = 4, common.legend = TRUE, legend = "bottom", align = "hv"
)

## ===== Group 5:contact zone between albosinese and ermanii=====
g5_ash   <- prep_plot_for_group(file_ash,   prefixes = prefix_group5, group_label = label_group5)
g5_bugg  <- prep_plot_for_group(file_bugg,  prefixes = prefix_group5, group_label = label_group5)
g5_cos   <- prep_plot_for_group(file_cos,   prefixes = prefix_group5, group_label = label_group5)
g5_lenta <- prep_plot_for_group(file_lenta, prefixes = prefix_group5, group_label = label_group5)

fig_group5 <- ggarrange(
  g5_ash, g5_bugg, g5_cos, g5_lenta,
  nrow = 1, ncol = 4, common.legend = TRUE, legend = "bottom", align = "hv"
)

print(fig_group1)
print(fig_group2)
print(fig_group3)
print(fig_group4)
print(fig_group5)

combined <- ggarrange(
  fig_group1, fig_group2, fig_group3,
  fig_group4, fig_group5,
  ncol = 1, nrow = 5,
  labels = c("A", "B", "C", "D", "E"),
  common.legend = F,
  legend = "bottom",
  align = "hv"
)

ggsave("five_groups_combined.pdf", combined, width = 15, height = 30)

