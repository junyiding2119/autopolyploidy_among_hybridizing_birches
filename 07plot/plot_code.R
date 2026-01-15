######## genomic composition (Fig. S15)
library(ggplot2)
library(mgcv)
library(tidyr)
library(dplyr)
library(ggpubr)
library(scales)

# Read geographic distribution data (lenta reference)

uti_pure <- c("LJS", "LJA", "DCX", "JS", "WZ", "DQ", "LP")
H1       <- c("RTX", "MYL", "BMX", "MEK", "SLS", "BX", "LFG", "JLX", "EBY", "YJA", "SDX","JD")
alb_pure <- c("SND", "XYB", "TBX", "DTX", "HP", "HHG", "TBS", "TZZ", "QLL", "LYC", "SWP", "BYS", "TTH", "LHS")
H2       <- c("LS" , "PQG", "BA", "XLA")
erm_pure <- c("LYL", "LLZ", "FHSJ", "CBS", "BJF")


allInd_lentaRef <- read.delim("clipboard",header = T)

get_prefix <- function(id) {
  prefix <- gsub("[0-9_-].*$", "", id)
  return(prefix)
}

assign_group <- function(id) {
  prefix <- get_prefix(id)
  if (prefix %in% uti_pure) return("uti_pure")
  if (prefix %in% H1) return("H1")
  if (prefix %in% alb_pure) return("alb_pure")
  if (prefix %in% H2) return("H2")
  if (prefix %in% erm_pure) return("erm_pure")
  return("Unknown")
}

allInd_lentaRef$group <- sapply(allInd_lentaRef$ID, assign_group)

# Set factor levels for consistent ordering
allInd_lentaRef$group <- factor(allInd_lentaRef$group, 
                                levels = c("uti_pure", "H1", "alb_pure", "H2", "erm_pure"))

# Define colors and borders
fill_colors <- c(
  "uti_pure" = "#a0cc58",
  "H1"       = "#a0cc58",
  "alb_pure" = "#da86b5",
  "H2"       = "#da86b5",
  "erm_pure" = "#f5d239"
)
border_colors <- c(
  "uti_pure" = "#a0cc58",
  "H1"       = "black",
  "alb_pure" = "#da86b5",
  "H2"       = "black",
  "erm_pure" = "#f5d239"
)

# PCA
pca <- prcomp(allInd_lentaRef[, c("lat", "long")], scale = TRUE)
allInd_lentaRef$PC1 <- pca$x[, 1]

# Plot
pdf("Relative_genomic_composition.pdf", width = 6, height = 6)

# cos plot
model_pca <- lm(cos ~ poly(PC1, 4), data = allInd_lentaRef)
summary(model_pca)
cos <- ggplot(allInd_lentaRef, aes(PC1, cos)) +
  geom_point(aes(fill = group, color = group), alpha = 0.5, shape = 21, size = 2, stroke = 0.5) +
  scale_fill_manual(values = fill_colors) +
  scale_color_manual(values = border_colors) +
  ylim(0, 0.6) +
  geom_smooth(method = "lm",
              formula = y ~ poly(x, 4),
              color = "#8A9CC4",
              se = FALSE) +
  geom_text(label = paste("R² = ", round(summary(model_pca)$adj.r.squared, 2)), x = -1.35, y = 0.55) +
  geom_text(label = "p < 2.2e-16 ", x = -1.3, y = 0.5) +
  theme_classic() +
  theme(axis.text.y = element_text(color = "#000000"), 
        axis.text.x = element_text(color = "#000000"),
        legend.position = "none")

# ash plot
model_pca <- lm(ash ~ poly(PC1, 4), data = allInd_lentaRef)
summary(model_pca)
ash <- ggplot(allInd_lentaRef, aes(PC1, ash)) +
  geom_point(aes(fill = group, color = group), alpha = 0.5, shape = 21, size = 2, stroke = 0.5) +
  scale_fill_manual(values = fill_colors) +
  scale_color_manual(values = border_colors) +
  ylim(0, 0.6) +
  geom_smooth(method = "lm",
              formula = y ~ poly(x, 4),
              color = "#F08961",
              se = FALSE) +
  geom_text(label = paste("R² = ", round(summary(model_pca)$adj.r.squared, 2)), x = -1.35, y = 0.55) +
  geom_text(label = "p < 2.2e-16 ", x = -1.3, y = 0.5) +
  theme_classic() +
  theme(axis.text.y = element_text(color = "#000000"), 
        axis.text.x = element_text(color = "#000000"),
        legend.position = "none")

# bug plot
model_pca <- lm(bug ~ poly(PC1, 1), data = allInd_lentaRef)
summary(model_pca)
bug <- ggplot(allInd_lentaRef, aes(PC1, bug)) +
  geom_point(aes(fill = group, color = group), alpha = 0.5, shape = 21, size = 2, stroke = 0.5) +
  scale_fill_manual(values = fill_colors) +
  scale_color_manual(values = border_colors) +
  ylim(0, 0.6) +
  geom_smooth(method = "lm",
              formula = y ~ poly(x, 1),
              color = "#E0C092",
              se = FALSE) +
  geom_text(label = paste("R² = ", round(summary(model_pca)$adj.r.squared, 2)), x = -1.35, y = 0.55) +
  geom_text(label = "p < 2.2e-16 ", x = -1.3, y = 0.5) +
  theme_classic() +
  theme(axis.text.y = element_text(color = "#000000"), 
        axis.text.x = element_text(color = "#000000"),
        legend.position = "none")

# cos_ash plot
model_pca <- lm(cos_ash ~ poly(PC1, 2), data = allInd_lentaRef)
summary(model_pca)
cos_ash <- ggplot(allInd_lentaRef, aes(PC1, cos_ash)) +
  geom_point(aes(fill = group, color = group), alpha = 0.5, shape = 21, size = 2, stroke = 0.5) +
  scale_fill_manual(values = fill_colors) +
  scale_color_manual(values = border_colors) +
  ylim(0, 0.6) +
  geom_smooth(method = "lm",
              formula = y ~ poly(x, 2),
              color = "#A0A0A0",
              se = FALSE) +
  geom_text(label = paste("R² = ", round(summary(model_pca)$adj.r.squared, 2)), x = -1.35, y = 0.55) +
  geom_text(label = "p < 2.2e-16 ", x = -1.3, y = 0.5) +
  theme_classic() +
  theme(axis.text.y = element_text(color = "#000000"), 
        axis.text.x = element_text(color = "#000000"),
        legend.position = "none")
ggarrange(cos, ash, bug, cos_ash, ncol = 2, nrow = 2)

dev.off()

######## pairing frequencies (Fig. 2a)
library(ggpubr)
library(tidyr)
library(dplyr)

new_labels <- c("B. lenta", "B. nana", "B. occidentalis", "B. populifolia",
                "B. pendula", "B. platyphylla", "B. mcallisteri", 
                "B. buggsii", "B. costata", "B. ashburneri", "cos_ash")

new_levels <- c("len", "nana", "occidentalis", "populifolia", 
                "pendula", "platyphylla", "mca", 
                "bug", "cos", "ash", "cos_ash")

###albosinensis
albo_lentaRef <- read.delim("clipboard",header = T,row.names = 1)
albo_lentaRef <- albo_lentaRef * 100

albo_lentaRef_long <- albo_lentaRef %>%
  pivot_longer(
    cols = c(pendula, platyphylla, populifolia, nana, occidentalis, cos, ash, bug, mca, cos_ash , len),
    names_to = "species",
    values_to = "value"
  )

albo_lentaRef_long$species <- factor(albo_lentaRef_long$species, 
                                     levels = new_levels,
                                     labels = new_labels)

albo <- ggbarplot(albo_lentaRef_long, x = "species", y = "value",
          add = "mean_se",  lab.vjust = -.7,
          ylab = "pairing	with B. albosinensis (%)",
          xlab = "", ylim = c(0, 40)) +
stat_summary(aes(label = round(..y.., 1)),
             geom = "text", fun = mean, vjust = -.7) +
theme(axis.text.x = element_text(angle = 45, hjust = 1,face = "italic"))

###ermanii
erman_lentaRef <- read.delim("clipboard",header = T,row.names = 1)
erman_lentaRef <- erman_lentaRef * 100

erman_lentaRef_long <- erman_lentaRef %>%
  pivot_longer(
    cols = c(pendula, platyphylla, populifolia, nana, occidentalis, cos, ash, bug,cos_ash, mca, len),
    names_to = "species",
    values_to = "value"
  )

erman_lentaRef_long$species <- factor(erman_lentaRef_long$species, 
                                      levels = new_levels,
                                      labels = new_labels)

erman <- ggbarplot(erman_lentaRef_long, x = "species", y = "value",
          add = "mean_se",  lab.vjust = -.7,
          ylab = "pairing	with B. ermanii (%)",
          xlab = "", ylim = c(0, 40)) +
  stat_summary(aes(label = round(..y.., 1)),
               geom = "text", fun = mean, vjust = -.7) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1,face = "italic"))

###utilis
uti_lentaRef <- read.delim("clipboard",header = T,row.names = 1)
uti_lentaRef <- uti_lentaRef * 100

uti_lentaRef_long <- uti_lentaRef %>%
  pivot_longer(
    cols = c(pendula, platyphylla, populifolia, nana, occidentalis, cos, ash, bug,cos_ash, mca, len),
    names_to = "species",
    values_to = "value"
  )

uti_lentaRef_long$species <- factor(uti_lentaRef_long$species, 
                                    levels = new_levels,
                                    labels = new_labels)

uti <- ggbarplot(uti_lentaRef_long, x = "species", y = "value",
          add = "mean_se",  lab.vjust = -.7,
          ylab = "pairing	with B. utilis (%)",
          xlab = "", ylim = c(0, 40)) +
  stat_summary(aes(label = round(..y.., 1)),
               geom = "text", fun = mean, vjust = -.7) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1,face = "italic"))

ggarrange(
  uti, albo, erman,
  nrow = 1,        
  ncol = 3,  
  common.legend = TRUE, 
  legend = "bottom",
  align = "hv"
)

######## visualizing gene ontology annotations and allele fractions (Fig. s20)
library(ggplot2)
library(tidyr)
library(dplyr)
library(stringr)
library(colorspace)
library(cowplot)


# 1. Load Data
# E:\paper\07 costatae evolution\data_analysis\alleles_geographic_distribution\01total
raw_data <- read.table("clipboard", header = TRUE, sep = "\t", stringsAsFactors = FALSE, quote = "'", check.names = FALSE)
colnames(raw_data)[1] <- "Gene_Symbol"
colnames(raw_data)[2] <- "GO_Keywords"

# 2. Define GO Terms, Categories, Colors
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


# 3. Data Preparation
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


# 4. Plotting
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

# 5. Generate and Combine Plots
plot_title_main <- "B. ashburneri origin"
go_plot_margins <- c(0.5, 0.5, 0.1, 0.5)

final_go_order_for_plot <- rev(level_order_go_user)

go_plot_data$Gene <- factor(go_plot_data$Gene, levels = level_order_gene)
go_plot_data$GO_term <- factor(go_plot_data$GO_term, levels = final_go_order_for_plot)

final_go_plot <- GOplotFUN_user(go_plot_data,
                                go_order = final_go_order_for_plot,
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


# 6. Calculate Category Frequencies
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

#View the Result
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

######## pairing frequencies for five pop (Fig. S17)

library(ggpubr)
library(tidyr)
library(dplyr)
library(stringr)
library(tibble)
library(grid)

setwd("E://paper/07 costatae evolution/data_analysis/polyAncestor20250329/")

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


######## Generate summary tables for five groups

library(dplyr)
library(tibble)
library(stringr)

setwd("E://paper/07 costatae evolution/data_analysis/polyAncestor20250329/")

cols_to_extract <- c("pendula", "platyphylla", "pend_platy", "populifolia", "pop_out", 
                     "nana", "nana_out", "occidentalis", "occ_out", 
                     "cos", "ash", "cos_ash", "bug", "bug_out", 
                     "mca", "mca_out", "len", "basal")

file_ash   <- "pairing_matrix_all_species_ash.txt"
file_bugg  <- "pairing_matrix_all_species_bugg.txt"
file_cos   <- "pairing_matrix_all_species_cos.txt"

ref_files <- list(
  "ash ref" = file_ash,
  "bugg ref" = file_bugg,
  "cos ref" = file_cos
)

groups <- list(
  list(
    name = "pure_utilis",
    label = "pure utilis",
    prefixes = c("LJS", "LJA", "DCX", "JS", "WZ", "DQ", "LP")
  ),
  list(
    name = "contact_zone_uti_albo",
    label = "contact zone between uti&albo",
    prefixes = c("RTX", "MYL", "BMX", "MEK", "SLS", "BX", "LFG", "JLX", "EBY", "YJA", "SDX", "JD")
  ),
  list(
    name = "pure_albosinensis",
    label = "pure albosinensis",
    prefixes = c("SND", "XYB", "TBX", "DTX", "HP", "HHG", "TBS", "TZZ", "QLL", "LYC", "SWP", "BYS", "TTH", "LHS")
  ),
  list(
    name = "contact_zone_albo_erman",
    label = "contact zone between albo&erman",
    prefixes = c("LS", "PQG", "BA", "XLA")
  ),
  list(
    name = "pure_ermanii",
    label = "pure ermanii",
    prefixes = c("LYL", "LLZ", "FHSJ", "CBS", "BJF")
  )
)

calculate_group_means <- function(file_path, prefixes) {
  df <- read.delim(file_path, header = TRUE, row.names = 1, check.names = FALSE)
  cols_present <- intersect(cols_to_extract, colnames(df))

  re <- paste0("^(", paste(prefixes, collapse = "|"), ")")
  df_sub <- df[grepl(re, rownames(df)), cols_present, drop = FALSE]
  
  if (nrow(df_sub) == 0) {
    warning(sprintf("No samples matched prefixes in %s", file_path))
    # Return a row of NAs
    result <- rep(NA, length(cols_to_extract))
    names(result) <- cols_to_extract
    return(result)
  }
  
  means <- colMeans(df_sub, na.rm = TRUE)
  
  result <- rep(NA, length(cols_to_extract))
  names(result) <- cols_to_extract
  result[names(means)] <- means
  
  return(result)
}

for (group in groups) {
  cat(sprintf("\nProcessing group: %s\n", group$label))
  
  result_matrix <- matrix(NA, nrow = length(ref_files), ncol = length(cols_to_extract))
  rownames(result_matrix) <- names(ref_files)
  colnames(result_matrix) <- cols_to_extract
  
  for (i in seq_along(ref_files)) {
    ref_name <- names(ref_files)[i]
    ref_file <- ref_files[[i]]
    
    cat(sprintf("  Calculating means for %s...\n", ref_name))
    means <- calculate_group_means(ref_file, group$prefixes)
    result_matrix[i, ] <- means
  }
  
  result_df <- as.data.frame(result_matrix)
  
  # Save to file
  output_file <- sprintf("summary_%s.txt", group$name)
  write.table(result_df, file = output_file, sep = "\t", 
              quote = FALSE, row.names = TRUE, col.names = NA)
  
  cat(sprintf("  Saved to: %s\n", output_file))
}

#### Linkage disequilibrium decay plot (Fig. S1)
library(ggplot2)
library(dplyr)

ash <- read.table('E://paper/07 costatae evolution/data_analysis/admixture/LD/ash_LD.AvgR2.txt', 
                  col.names = c('Sample', 'Ploidy', 'Start', 'End', 'R2', 'Variance'))
cos <- read.table('E://paper/07 costatae evolution/data_analysis/admixture/LD/cos_LD.AvgR2.txt', 
                  col.names = c('Sample', 'Ploidy', 'Start', 'End', 'R2', 'Variance'))
bug <- read.table('E://paper/07 costatae evolution/data_analysis/admixture/LD/bug_LD.AvgR2.txt', 
                  col.names = c('Sample', 'Ploidy', 'Start', 'End', 'R2', 'Variance'))

data <- rbind(ash, cos, bug)

data$Distance <- (data$Start + data$End) / 2 / 1000

data$Sample <- gsub("_LD", "", data$Sample)

p <- ggplot(data, aes(x = Distance, y = R2, color = Sample, group = Sample)) +
  geom_line(size = 1, alpha = 0.8) +
  geom_point(size = 1.5, alpha = 0.6) +
  scale_x_log10(breaks = c( 1, 10, 100, 1000, 10000),
                labels = c( "1", "10", "100", "1,000", "10,000")) +
  scale_color_manual(values = c("ash" = "#e6855f", 
                                "cos" = "#8798bd", 
                                "bug" = "#d4b78b")) +
  labs(x = "Distance (kb)",
       y = expression(paste("Genotypic correlation (", r^2, ")")),
       color = "Sample") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        legend.position = c(0.85, 0.85),
        legend.background = element_rect(fill = "white", color = "black"),
        axis.text = element_text(size = 11),
        axis.title = element_text(size = 12, face = "bold"),
        legend.title = element_text(size = 11, face = "bold"))

print(p)

erman <- read.table('E://paper/07 costatae evolution/data_analysis/admixture/LD/erman_LD.AvgR2.txt', 
                    col.names = c('Sample', 'Ploidy', 'Start', 'End', 'R2', 'Variance'))
albo <- read.table('E://paper/07 costatae evolution/data_analysis/admixture/LD/albo_LD.AvgR2.txt', 
                   col.names = c('Sample', 'Ploidy', 'Start', 'End', 'R2', 'Variance'))
uti <- read.table('E://paper/07 costatae evolution/data_analysis/admixture/LD/uti_LD.AvgR2.txt', 
                  col.names = c('Sample', 'Ploidy', 'Start', 'End', 'R2', 'Variance'))
erman$Sample <- "erman"
albo$Sample <- "albo"
uti$Sample <- "uti"


data <- rbind(erman, albo, uti)

data$Distance <- (data$Start + data$End) / 2 / 1000

p <- ggplot(data, aes(x = Distance, y = R2, color = Sample, group = Sample)) +
  geom_line(size = 1, alpha = 0.8) +
  geom_point(size = 1.5, alpha = 0.6) +
  scale_x_log10(breaks = c(0.1, 1, 10, 100, 1000, 10000),
                labels = c("0.1", "1", "10", "100", "1,000", "10,000")) +
  scale_color_manual(values = c("erman" = "#ebcb42", 
                                "albo" = "#da86b5", 
                                "uti" = "#9bc459")) +
  labs(x = "Distance (kb)",
       y = expression(paste("Genotypic correlation (", r^2, ")")),
       color = "Sample") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        legend.position = c(0.85, 0.85),
        legend.background = element_rect(fill = "white", color = "black"),
        axis.text = element_text(size = 11),
        axis.title = element_text(size = 12, face = "bold"),
        legend.title = element_text(size = 11, face = "bold"))
print(p)

# Matrix of relative (FST) and absolute (dxy) divergence for pairwise species comparisons (Fig 1)
library(tidyr)
library(dplyr)
library(ggplot2)

raw_data <- read.delim("clipboard", header = TRUE)
#Species1	Species2	Fst	Dxy
#2ash	1bug	0.354499789	0.010538055
#5albo	1bug	0.383346607	0.011484311
#5albo	2ash	0.157789818	0.008924637
#1bug	3cos	0.432031141	0.011400987
#2ash	3cos	0.235463755	0.009139384
#5albo	3cos	0.264025804	0.009802776
#1bug	4uti	0.422196999	0.011310054
#2ash	4uti	0.18868109	0.008745591
#3cos	4uti	0.312861605	0.00986099
#5albo	4uti	0.102383475	0.008216603
#6erm	4uti	0.184479375	0.009817451
#1bug	6erm	0.394596229	0.012460632
#2ash	6erm	0.239526963	0.010519845
#3cos	6erm	0.196305342	0.009579332
#5albo	6erm	0.146961164	0.009652669


species <- c("1bug", "2ash", "3cos", "4uti", "5albo", "6erm")


target <- data.frame()
for (j in 1:length(species)) {     
  for (i in 1:length(species)) {  
    if (i == j) {
      pos <- "diagonal"
    } else if (i < j) {
      pos <- "Dxy"
    } else {
      pos <- "Fst"
    }
    
    target <- rbind(target, data.frame(
      Species1 = species[i],
      Species2 = species[j],
      value = NA,
      pos = pos
    ))
  }
}

print(target)

fill_value <- function(sp1, sp2, pos, raw_data) {
  if (pos == "diagonal") return(NA)
  match1 <- raw_data[(raw_data$Species1 == sp1 & raw_data$Species2 == sp2), ]
  match2 <- raw_data[(raw_data$Species1 == sp2 & raw_data$Species2 == sp1), ]
  if (nrow(match1) > 0) {
    return(ifelse(pos == "Fst", match1$Fst, match1$Dxy))
  } else if (nrow(match2) > 0) {
    return(ifelse(pos == "Fst", match2$Fst, match2$Dxy))
  } else {
    return(NA)
  }
}

for (i in 1:nrow(target)) {
  if (target$pos[i] != "diagonal") {
    target$value[i] <- fill_value(target$Species1[i], target$Species2[i], 
                                  target$pos[i], raw_data)
  }
}

combined_long <- read.delim("clipboard",header=T)
#plot
ggplot(combined_long, aes(x = Species1, y = Species2, fill = value)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.3f", value)),
            color = "black", size = 3.5, na.rm = TRUE) +
  scale_fill_gradient2(low = "#C5DAF6FF", mid = "#6996E3FF", high = "#1A318BFF",
                       midpoint = median(combined_long$value, na.rm = TRUE),
                       na.value = "grey90",
                       name = "Value") +
  scale_y_discrete(limits = rev) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 11, face = "bold"),
    axis.text.y = element_text(size = 11, face = "bold"),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
  ) +
  coord_fixed()

##############
###Fig2 e
##############

library(ggplot2)
library(dplyr)
library(tidyr)

raw_data <- read.delim("E://paper/07 costatae evolution/data_analysis/05s_ABC_geographic_distribution/Fig2_e_data.txt")

long_data <- raw_data %>%
  pivot_longer(cols = c(ashburneri, costata, buggsii, ancestor),
               names_to = "component",
               values_to = "value") %>%
  mutate(value = value * 100) 

# mean se
summary_data <- long_data %>%
  group_by(Population, component) %>%
  summarise(
    mean = mean(value),
    sd = sd(value),
    se = sd(value) / sqrt(n()), 
    .groups = 'drop'
  )

summary_data$component <- factor(summary_data$component,
                                 levels = c("ashburneri", "costata", 
                                            "buggsii", "ancestor"))

# plot
p <- ggplot(summary_data, aes(x = Population, y = mean, fill = component)) +
  geom_bar(stat = "identity", position = position_dodge(0.9), width = 0.8) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se),
                position = position_dodge(0.9),
                width = 0.25) +
  geom_text(aes(label = round(mean, 1)),
            position = position_dodge(0.9),
            vjust = -0.9, size = 3) + 
  scale_fill_manual(values = c("ashburneri" = "#E07856",
                               "costata" = "#E8C547",
                               "buggsii" = "#D4B896",
                               "ancestor" = "#A1A1A1"),
                    labels = c("B. ashburneri", "B. costata", 
                               "B. buggsii", "Ancestor")) +
  labs(x = NULL, y = "Genomic component (%)", fill = "component") +
  theme_classic() +
  theme(
    legend.position = "right",
    axis.text = element_text(size = 11),
    axis.title = element_text(size = 12)
  ) +
  ylim(0, 105) 

ggsave("plot.pdf", p, width = 10, height = 4)
