### Packages
library(data.table)
library(dplyr)
library(tibble)
library(ggpubr)
library(stringr)
library(tidyr)
library(caret)
library(tidyverse)
library(igraph)
library(ggraph)
library(visNetwork)
library(scales)
library(ggrepel)
library(ComplexUpset)
library(patchwork)
library(cowplot)
library(ggplotify)
library(RColorBrewer)
library(colorspace)
library(ggforce)
library(purrr)
library(ggnewscale)

### Data
matches <- read_tsv("/ENTER YOUR PATH HERE/output/data/Batch_MASST_matches.txt")
ReDu <-fread("/ENTER YOUR PATH HERE/all_sampleinformation.tsv")
metab <- readRDS("ranked_metabolites.rds")
metab <-rename(metab, Scan = Metabolite)
metabanno  <- fread("annotated_molecules_edit.csv")
datasets <- read_tsv("/ENTER YOUR PATH HERE/output/data/Batch_MASST_datasets.txt")

matches %>% summarise(n_unique = n_distinct(USI))
datasets %>% summarise(n_unique = n_distinct(Dataset))

matches$USI <- str_remove(matches$USI, ":scan:\\d+$")

df1 <- matches %>%
  left_join(ReDu, by = "USI")
df1$Scan<-as.character(df1$Scan)

df <- df1 %>%
  left_join(metab, by = "Scan")

##### Disease state analysis
### Dataframe for disease plots, which excludes all feature matches to A) non-humans, B) healthy samples, C) missing disease classifiers D) Not blood samples
df_clean <- df %>%
  filter(
    !is.na(DOIDCommonName),
    DOIDCommonName != "missing value",
    !str_detect(NCBITaxonomy, "10090|Mus musculus|10088|Mus|missing value"),
    str_detect(UBERONBodyPartName, "blood plasma|blood serum|blood"),
    !str_detect(HealthStatus, "healthy")
  ) %>%
  group_by(Scan) %>%
  filter(n_distinct(USI) >= 2) %>%   # Keep only scans with ≥2 unique USIs
  ungroup() %>%
  dplyr::select(
    USI, Scan, Cosine, ATTRIBUTE_DatasetAccession, DOIDCommonName,
    HealthStatus, Direction, Rank, Weight, NCBITaxonomy, UBERONBodyPartName
  )

#Checking number of matches
length(unique(df_clean$USI))  
length(unique(df_clean$Scan)) 

### Creating disease groups
# Create a disease-group mapping
disease_groups <- data.frame(
  DOIDCommonName = c("acquired immunodeficiency syndrome", 
                     "human immunodeficiency virus infectious disease",
                     "COVID-19",
                     "Chagas disease",
                     "primary bacterial infectious disease",
                     "toxoplasmosis",
                     "Crohn's disease",
                     "inflammatory bowel disease",
                     "Kawasaki disease",
                     "psoriasis",
                     "rheumatoid arthritis",
                     "ulcerative colitis",
                     "diabetes mellitus",
                     "hypertension",
                     "metabolic dysfunction-associated steatotic liver disease",
                     "obesity",
                     "chronic obstructive pulmonary disease",
                     "cystic fibrosis",
                     "pulmonary fibrosis",
                     "Alzheimer's disease",
                     "schizophrenia",
                     "dental caries",
                     "mucopolysaccharidosis",
                     "osteoarthritis"
  ),
  Group = c("Immunodeficiency", 
            "Immunodeficiency",
            "Infection", 
            "Infection",
            "Infection",
            "Infection",
            "Inflammatory/Autoimmune",
            "Inflammatory/Autoimmune",
            "Inflammatory/Autoimmune",
            "Inflammatory/Autoimmune",
            "Inflammatory/Autoimmune",
            "Inflammatory/Autoimmune",
            "Metabolic",
            "Metabolic",
            "Metabolic",
            "Metabolic",
            "Pulmonal",
            "Pulmonal",
            "Pulmonal",
            "Neurological",
            "Psychiatric",
            "Other",
            "Other",
            "Other"
  ))

### Counting and order cleaning
Disease_direction_clean <- df_clean %>%
  count(DOIDCommonName, Direction, name = "Count")

#Adding percentages
Disease_direction_clean <- Disease_direction_clean %>%
  group_by(DOIDCommonName) %>%
  mutate(Percent = Count / sum(Count) * 100) %>%
  ungroup()

# Ensure mapping is correct
Disease_direction_grouped_clean <- Disease_direction_clean %>%
  left_join(disease_groups, by = "DOIDCommonName") %>%
  filter(!is.na(Group))


### Plotting (Figure supplementary X)
# Order diseases within each group for consistent plotting
Disease_order <- Disease_direction_grouped_clean %>%
  group_by(Group) %>%
  summarise(Disease_order = list(unique(DOIDCommonName))) %>%
  deframe()

#Count total samples available per disease that are in our matched diseases
Disease_degree_all_clean <- ReDu %>%
  filter(
    !is.na(DOIDCommonName),
    DOIDCommonName != "missing value",
    str_detect(NCBITaxonomy, "9606|Homo sapiens"),
    str_detect(UBERONBodyPartName, "blood plasma|blood serum|blood"),
    !str_detect(HealthStatus, "healthy"),
    !str_detect(DOIDCommonName, "primary bacterial infectious disease|toxoplasmosis|inflammatory bowel disease|psoriasis|ulcerative colitis|diabetes mellitus|hypertension|obesity|chronic obstructive pulmonary disease|cystic fibrosis|pulmonary fibrosis|schizophrenia|dental caries|mucopolysaccharidosis")
  ) %>%
  count(DOIDCommonName, name = "Count") %>%
  left_join(disease_groups, by = "DOIDCommonName") %>%
  filter(!is.na(Group))

#Count number of libraries  
Library_degree_all_clean <- ReDu %>%
  filter(
    !is.na(DOIDCommonName),
    DOIDCommonName != "missing value",
    str_detect(NCBITaxonomy, "9606|Homo sapiens"),
    str_detect(UBERONBodyPartName, "blood plasma|blood serum|blood"),
    !str_detect(HealthStatus, "healthy"),
    !str_detect(DOIDCommonName, "primary bacterial infectious disease|toxoplasmosis|inflammatory bowel disease|psoriasis|ulcerative colitis|diabetes mellitus|hypertension|obesity|chronic obstructive pulmonary disease|cystic fibrosis|pulmonary fibrosis|schizophrenia|dental caries|mucopolysaccharidosis")
  ) %>%
  group_by(DOIDCommonName) %>%                   # group by disease
  summarise(LibraryCount = n_distinct(ATTRIBUTE_DatasetAccession)) %>%  # count unique libraries
  left_join(disease_groups, by = "DOIDCommonName") %>%
  filter(!is.na(Group))

#Code below is used to count features detected
Disease_degree_clean <- df_clean %>%
  count(DOIDCommonName, name = "Count") %>%
  left_join(disease_groups, by = "DOIDCommonName") %>%
  filter(!is.na(Group))


# --- HEATMAP ---
heatmap_disease <- ggplot(Disease_direction_grouped_clean,
                          aes(x = Direction, y = DOIDCommonName, fill = Percent)) +
  geom_tile(color = "white", width = 0.9, height = 0.9) +
  geom_text(aes(label = paste0(round(Percent, 1), "%")),
            color = "black", size = 5, fontface = "bold") +
  facet_grid(Group ~ ., scales = "free_y", space = "free_y", switch = "y") +
  scale_fill_gradient(low = "white", high = "#d95f02") +
  theme_minimal(base_size = 18) +
  labs(x = "Direction", y = NULL, fill = "Percent") +
  theme(
    axis.text.y = element_text(face = "bold"),
    axis.text.x = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    panel.grid = element_blank(),
    strip.placement = "outside",
    strip.background = element_blank(),
    strip.text.y.left = element_text(face = "bold", size = 18, angle = 0)
  )


# --- HISTOGRAM (flipped, grouped, and ordered properly) ---

# Reuse the same grouping and facet structure
histogram_disease <- ggplot(Disease_degree_all_clean, aes(y = DOIDCommonName, x = Count)) +
  geom_col(fill = "#d95f02") +
  geom_text(aes(label = Count),
            hjust = -0.2, color = "black", size = 5, fontface = "bold") +
  facet_grid(Group ~ ., scales = "free_y", space = "free_y") +
  theme_minimal(base_size = 18) +
  labs(x = "Count", y = NULL) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid = element_blank(),
    strip.text.y = element_blank(),
    axis.title = element_text(face = "bold"),
    plot.margin = margin(t = 5, r = 10, b = 5, l = 5)
  ) +
  coord_cartesian(clip = "off", xlim = c(0, max(Disease_degree_all_clean$Count) * 1.1))

#histogram for libraries
histogram_library <- ggplot(Library_degree_all_clean, aes(y = DOIDCommonName, x = LibraryCount)) +
  geom_col(fill = "#d95f02") +
  geom_text(aes(label = LibraryCount),
            hjust = -0.2, color = "black", size = 5, fontface = "bold") +
  facet_grid(Group ~ ., scales = "free_y", space = "free_y") +
  theme_minimal(base_size = 18) +
  labs(x = "Count", y = NULL) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid = element_blank(),
    strip.text.y = element_blank(),
    axis.title = element_text(face = "bold"),
    plot.margin = margin(t = 5, r = 10, b = 5, l = 5)
  ) +
  coord_cartesian(clip = "off", xlim = c(0, max(Library_degree_all_clean$LibraryCount) * 1.1))

#Individual plots
heatmap_disease
histogram_disease
histogram_library

# --- ALIGN PLOTS ---

# Convert ggplots to grobs
heatmap_grob <- as_grob(heatmap_disease)
histogram_grob <- as_grob(histogram_disease)
library_grob <- as_grob(histogram_library)

# Combine with cowplot, aligning by rows
combined_plot <- plot_grid(
  heatmap_grob,
  histogram_grob,
  library_grob,
  ncol = 3,
  align = "v",
  axis = "tb",
  rel_widths = c(3, 1, 1)  # adjust width ratio as needed
)

combined_plot

#### Body site analysis
### Body site data, excluding those not matched to human samples
df_bodypart <- df %>%
  filter(
    str_detect(NCBITaxonomy, "9606|Homo sapiens")) %>%
  group_by(Scan) %>%
  filter(n_distinct(USI) >= 2) %>%   # Keep only scans with ≥2 unique USIs
  ungroup() %>%
  select(
    USI, Scan, Cosine, ATTRIBUTE_DatasetAccession, DOIDCommonName,
    HealthStatus, Direction, Rank, Weight, NCBITaxonomy, UBERONBodyPartName, LifeStage
  )

#Prepping anatomy maps
sample_map_2 <- c(
  "arm skin" = "skin",
  "axilla skin" = "skin",
  "head or neck skin" = "skin",
  "skin of pes" = "skin",
  "skin of manus" = "skin",
  "skin of leg" = "skin",
  "skin of trunk" = "skin",
  "zone of skin" = "skin",
  "alveolar system" = "lung",
  "colostrum" = "milk",
  "corpus callosum" = "brain",
  "forebrain" = "brain",
  "proximal tubule" = "kidney",
  "renal tubule" = "kidney",
  "supragingival dental plaque" = "dental plaque",
  "venous blood" = "blood",
  "subcutaneous adipose tissue" = "adipose tissue"
  # Add all other detailed names here
)

anatomical_order <- rev(c(
  "brain", #1
  "cerebrospinal fluid", #2
  "oral cavity", #3
  "saliva", #4
  "dental plaque", #5
  "lung", #6
  "heart", #7
  "milk", #8
  "liver", #9
  "skin", #10
  "adipose tissue", #11
  "epithelial cell", #12
  "blood", #13
  "blood plasma", #14
  "blood serum", #14
  "kidney", #15
  "urine", #16
  "feces" #17
))

# Collapse UBERONBodyPartName
df_bodypart <- df_bodypart %>%
  mutate(UBERONBodyPartName = ifelse(
    UBERONBodyPartName %in% names(sample_map_2),
    sample_map_2[UBERONBodyPartName],
    UBERONBodyPartName
  ))

### Counting and order cleaning
#Counts
Body_direction_clean <- df_bodypart %>%
  count(UBERONBodyPartName, Direction, name = "Count")

#Percentages
Body_direction_clean <- Body_direction_clean %>%
  group_by(UBERONBodyPartName) %>%
  mutate(Percent = Count / sum(Count) * 100) %>%
  ungroup()

### Removing sample types with < 500 matches 
# Count sample types
Library_segments_body <- ReDu %>%
  filter(str_detect(NCBITaxonomy, "9606|Homo sapiens"),
         !str_detect(UBERONBodyPartName, "^tissue$"),
         !str_detect(UBERONBodyPartName,
                     "2 cell stage|B cell|T cell|adipocyte|bone marrow|bone tissue|breast|bronchial mucosa|cervical mucus|colon|epithelium of mammary gland|esophagus mucosa|eye|fibroblast|foreskin fibroblast|granulosa cell|hepatocyte|islet of Langerhans|keratinocyte|leaf|leukocyte|macrophage|mesenchymal stem cell|monocyte|mononuclear cell|myofibroblast cell|neuron|neutrophil|peripheral blood mononuclear cell|peritoneal fluid|pluripotent stem cell|prostate gland|rectal lumen|sebum|seminal fluid|skeletal muscle tissue|sputum|stem cell|stomach|sweat|synovial joint|thymus|thyroid gland|uterine cervix|vitreous humor|white adipose tissue|missing value")
  ) %>%
  count(UBERONBodyPartName, ATTRIBUTE_DatasetAccession, name = "Count") 


# Collapse UBERONBodyPartName
Library_segments_body <- Library_segments_body %>%
  mutate(UBERONBodyPartName = ifelse(
    UBERONBodyPartName %in% names(sample_map_2),
    sample_map_2[UBERONBodyPartName],
    UBERONBodyPartName
  ))

#Collapse per body part
Library_segments_body <- Library_segments_body %>%
  group_by(UBERONBodyPartName) %>%
  summarise(Count = sum(Count), .groups = "drop")%>%
  filter(Count >= 500) #Removing those with less than 500 samples


#count libraries
Library_body_all_clean <- ReDu %>%
  filter(str_detect(NCBITaxonomy, "9606|Homo sapiens"),
         !str_detect(UBERONBodyPartName, "^tissue$"),
         !str_detect(UBERONBodyPartName,
                     "2 cell stage|B cell|T cell|adipocyte|bone marrow|bone tissue|breast|bronchial mucosa|cervical mucus|colon|epithelium of mammary gland|esophagus mucosa|eye|fibroblast|foreskin fibroblast|granulosa cell|hepatocyte|islet of Langerhans|keratinocyte|leaf|leukocyte|macrophage|mesenchymal stem cell|monocyte|mononuclear cell|myofibroblast cell|neuron|neutrophil|peripheral blood mononuclear cell|peritoneal fluid|pluripotent stem cell|prostate gland|rectal lumen|sebum|seminal fluid|skeletal muscle tissue|sputum|stem cell|stomach|sweat|synovial joint|thymus|thyroid gland|uterine cervix|vitreous humor|white adipose tissue|missing value")
  ) %>%
  group_by(UBERONBodyPartName) %>%                   # group by body part
  summarise(LibraryCount = n_distinct(ATTRIBUTE_DatasetAccession))  # count unique libraries

# Collapse UBERONBodyPartName
Library_body_all_clean <- Library_body_all_clean %>%
  mutate(UBERONBodyPartName = ifelse(
    UBERONBodyPartName %in% names(sample_map_2),
    sample_map_2[UBERONBodyPartName],
    UBERONBodyPartName
  ))

#Collapse per body part
Library_body_all_clean <- Library_body_all_clean %>%
  group_by(UBERONBodyPartName) %>%
  summarise(LibraryCount = sum(LibraryCount), .groups = "drop")


# Keep only body parts that had >= 500 counts in Library_segments_body in the library count
Library_body_all_clean <- Library_body_all_clean %>%
  filter(UBERONBodyPartName %in% Library_segments_body$UBERONBodyPartName)

# Keep only body parts that had >= 500 counts in Library_segments_body in the heatmap
Body_direction_clean_filtered <- Body_direction_clean %>%
  filter(UBERONBodyPartName %in% Library_segments_body$UBERONBodyPartName)

### Ordering for plotting
#Ordering by anatomy
Body_direction_clean_filtered$UBERONBodyPartName <- factor(
  Body_direction_clean_filtered$UBERONBodyPartName,
  levels = anatomical_order
)
Library_segments_body$UBERONBodyPartName <- factor(
  Library_segments_body$UBERONBodyPartName,
  levels = anatomical_order
)
Library_body_all_clean$UBERONBodyPartName <- factor(
  Library_body_all_clean$UBERONBodyPartName,
  levels = anatomical_order
)

#Plotting (Figure XA and B)
#heatmap
heatmap_body_2 <- ggplot(Body_direction_clean_filtered, aes(x = Direction, y = UBERONBodyPartName, fill = Percent)) +
  geom_tile(color = "white", width = 0.9, height = 0.9) +
  geom_text(aes(label = paste0(round(Percent, 1), "%")),
            color = "black", size = 4.5, fontface = "bold") +
  scale_fill_gradient(low = "white", high = "#d95f02") +
  theme_minimal(base_size = 16) +
  labs(x = "Direction", y = NULL, fill = "Percent") +
  theme(
    axis.text.y = element_text(face = "bold"),
    axis.text.x = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    panel.grid = element_blank()
  )

# --- HISTOGRAM (flipped, grouped, and ordered properly) ---

# Reuse the same grouping and facet structure
histogram_sample_body <- ggplot(Library_segments_body, aes(y = UBERONBodyPartName, x = Count)) +
  geom_col(fill = "#d95f02") +
  geom_text(aes(label = Count),
            hjust = -0.2, color = "black", size = 4.5, fontface = "bold") +
  theme_minimal(base_size = 12) +
  labs(x = "Sample count", y = NULL) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid = element_blank(),
    strip.text.y = element_blank(),
    axis.title = element_text(face = "bold"),
    plot.margin = margin(t = 0, r = 0, b = 0, l = 0)
  ) +
  coord_cartesian(clip = "off", xlim = c(0, max(Library_segments_body$Count) * 1.1))

#histogram for libraries
histogram_library_body <- ggplot(Library_body_all_clean, aes(y = UBERONBodyPartName, x = LibraryCount)) +
  geom_col(fill = "#d95f02") +
  geom_text(aes(label = LibraryCount),
            hjust = -0.2, color = "black", size = 4.5, fontface = "bold") +
  theme_minimal(base_size = 12) +
  labs(x = "Number of Datasets", y = NULL) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid = element_blank(),
    strip.text.y = element_blank(),
    axis.title = element_text(face = "bold"),
    plot.margin = margin(t = 0, r = 0, b = 0, l = 0)
  ) +
  coord_cartesian(clip = "off", xlim = c(0, max(Library_body_all_clean$LibraryCount) * 1.1))



# --- ALIGN PLOTS ---

# Convert ggplots to grobs
heatmap_grob_body <- as_grob(heatmap_body_2)
histogram_grob_body <- as_grob(histogram_sample_body)
library_grob_body <- as_grob(histogram_library_body)

# Combine with cowplot, aligning by rows
combined_plot_body_2 <- plot_grid(
  heatmap_grob_body,
  histogram_grob_body,
  library_grob_body,
  ncol = 3,
  align = "v",
  axis = "tb",
  rel_widths = c(1, 1, 1)  # adjust width ratio as needed
)

combined_plot_body_2
