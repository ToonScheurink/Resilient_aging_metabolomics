### Packages
library(data.table)
library(dplyr)
library(tidyr)
library(tidyverse)
library(purrr)
library(pheatmap)
library(gridExtra)
library(ComplexUpset)
library(ggtext)

### Loading Data
#MicrobeMASST
df <- fread("microbes-CR-taxonomy.csv")
df$Scan <- as.character(df$Scan)

#Features from CR
model_metab <- readRDS("ranked_metabolites.rds")
model_metab <-dplyr::rename(model_metab, Scan = Metabolite)

#annotations
metabanno  <- fread("annotated_molecules_edit.csv")
metabanno <-dplyr::rename(metabanno, Scan = scanID)
metabanno$Scan <- as.character(metabanno$Scan)

df_clean <- df %>%
  filter(QC != "Yes" & Blank != "Yes")

### Counting microbial features, exporting unique scans for later use
# Removing samples that are either QC or Blank
unique_scans <- unique(df_clean$Scan)
unique_scans

#write.csv(unique_scans, "unique_microbial_scans.csv", row.names = FALSE)

# Counts per taxonomic level
tax_cols <- c("kingdom", "phylum", "class", "order", "family", "genus", "species")
sapply(df_clean[, ..tax_cols], function(col) length(unique(col[!is.na(col)])))

### Counts per phylum
df_count <- df_clean %>%
  group_by(Scan, phylum) %>%
  summarise(count = n(), .groups = "drop")

df_wide_count <- df_count %>%
  tidyr::pivot_wider(
    names_from = phylum,
    values_from = count,
    values_fill = 0
  ) %>%
  mutate(total_matches = rowSums(select(., -Scan)))

### Counts per class
df_count_class <- df_clean %>%
  group_by(Scan, Direction, class) %>%
  summarise(count = n(), .groups = "drop")

df_count_class_for_wide <- df_clean %>%
  group_by(Scan, class) %>%
  summarise(count = n(), .groups = "drop")

df_wide_count_class <- df_count_class_for_wide %>%
  tidyr::pivot_wider(
    names_from = class,
    values_from = count,
    values_fill = 0
  ) %>%
  mutate(total_matches = rowSums(select(., -Scan)))

### Plot (Figure XC)
#Obtain phylym info
df_phylum <- df_clean %>%
  select(phylum, class) 

df_phylum <- df_phylum %>% 
  distinct(phylum, class)

# Create df with absence/prescence counts of scans per per microbe
df_class_direction <- df_clean %>%
  distinct(Scan, class, Direction) %>%  # keep only unique ScanID–class–Direction combos
  group_by(class, Direction) %>%
  summarise(count = n(), .groups = "drop")

# Join phylum information to your plotting dataframe
df_plot <- df_class_direction %>%
  left_join(df_phylum, by = "class") %>%
  # remove cases where class is NA but phylum is not
  filter(!(is.na(class) & !is.na(phylum))) %>%
  mutate(phylum = ifelse(is.na(phylum), "Z Unclassified class", phylum)) %>%
  arrange(phylum, class) %>%
  mutate(class = factor(class, levels = unique(class)))

# Create a factor for class ordered by phylum, then class (both alphabetically)
df_plot <- df_plot %>%
  arrange(phylum, class) %>%
  mutate(class = factor(class, levels = unique(class)))

# Numeric x positions for class factor
df_plot <- df_plot %>%
  mutate(class_num = as.numeric(class))

# numeric x positions for class factor
df_plot <- df_plot %>%
  mutate(class = factor(class, levels = unique(class)),
         class_num = as.numeric(class))

# compute phylum ranges in the same order as they appear on the x-axis
phylum_ranges <- df_plot %>%
  group_by(phylum) %>%
  summarize(xmin = min(class_num) - 0.5,
            xmax = max(class_num) + 0.5,
            .groups = "drop") %>%
  # order by factor levels as in plot
  arrange(match(phylum, unique(df_plot$phylum))) %>%
  mutate(bg_fill = rep(c("grey85", "white"), length.out = n()))

# Plot
microbes_counts_class_2<-ggplot() +
  # background shading
  geom_rect(data = phylum_ranges, 
            aes(xmin = xmin, xmax = xmax, ymin = 0, ymax = Inf),
            fill = phylum_ranges$bg_fill,
            alpha = 0.6,
            inherit.aes = FALSE) +
  # bars
  geom_bar(data = df_plot, aes(x = class, y = count, fill = Direction), stat = "identity") +
  # horizontal grid lines
  scale_fill_manual(values = c("Positive" = "#583F8F", "Negative" = "#E1642F")) +
  labs(
    title = "",
    x = "",
    y = "Number of Matches"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    panel.grid = element_blank(),         
    axis.title.y = element_blank()        
  )

microbes_counts_class_2

