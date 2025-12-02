#Packages
library(data.table)
library(dplyr)
library(tibble)
library(ggpubr)
library(stringr)
library(tidyr)
library(caret)
library(mixOmics)
library(tidyverse)
library(vegan)
library(ggrepel)
library(scales)

#transposition function
transpose_and_format_df <- function(df) {
  # Transpose the dataframe
  df_t <- t(df)
  df_t <- as.data.frame(df_t)  # Ensure it's a dataframe
  
  # Use the first row to name the columns
  col_names <- as.character(unlist(df_t[1, ]))
  df_t <- df_t[-1, ]  # Remove the first row used for column names
  colnames(df_t) <- col_names
  
  # Make the row index the first column
  df_t <- cbind(filename = row.names(df_t), df_t)
  
  # Optionally, remove the row names
  row.names(df_t) <- NULL
  
  return(df_t)
}

###read the metadata from qiita, and reformat sample names to the metadata format###
qiita_metadata <- fread("Qiita_metadata.csv")
qiita_metadata$sample_name <- gsub("^\\d+\\.", "", qiita_metadata$sample_name)  # Remove the numeric prefix
qiita_metadata$sample_name <- gsub("\\.", "_", qiita_metadata$sample_name)      # Replace periods with underscores

#prepare the MS sequence files. This is the order in which the sequences were run. 
seq1 <- fread("Sequence_20241209_Wellcome_Plate1.csv")
seq2 <- fread("Sequence_20241210_Wellcome_Plate2.csv")
seq3 <- fread("Sequence_20241212_Wellcome_Plate3.csv")
seq4 <- fread("Sequence_20241214_Wellcome_Plate4.csv")
seq5 <- fread("Sequence_20241215_Wellcome_Plate5.csv")
seq6 <- fread("Sequence_20241216_Wellcome_Plate6.csv")
seq7 <- fread("Sequence_20241218_Wellcome_Plate7.csv")
seq8 <- fread("Sequence_20241220_Wellcome_Plate8.csv")
seq9 <- fread("Sequence_20241223_Wellcome_Plate9.csv")

##List of sequence dataframes
sequence_files <- list(seq1, seq2, seq3, seq4, seq5, seq6, seq7, seq8, seq9)

# List of plate names corresponding to each dataframe
plate_names <- c("Plate1", "Plate2", "Plate3", "Plate4", "Plate5", "Plate6", "Plate7", "Plate8", "Plate9")

# Add "Plate" column to each dataframe -> adding plate number to the sequence files 
for (i in seq_along(sequence_files)) {
  sequence_files[[i]]$Plate <- plate_names[i]
}
# Assign the updated dataframes back to their original variables
list2env(setNames(sequence_files, paste0("seq", 1:9)), envir = .GlobalEnv)

##Add another column called Sequence in all the seq files
# Add `Sample_order` column to each dataframe
for (i in seq_along(sequence_files)) {
  sequence_files[[i]]$Sample_order <- paste0(sequence_files[[i]]$Plate, "_", seq_len(nrow(sequence_files[[i]])))
}
# Assign the updated dataframes back to their original variables
list2env(setNames(sequence_files, paste0("seq", 1:9)), envir = .GlobalEnv)

##Cleaning up the columns not needed
# Keep only columns 2, 22, and 23 in each dataframe
for (i in seq_along(sequence_files)) {
  sequence_files[[i]] <- sequence_files[[i]][, c(2, 22, 23)]
}
# Assign the updated dataframes back to their original variables
list2env(setNames(sequence_files, paste0("seq", 1:9)), envir = .GlobalEnv)

# Combine all dataframes into one
stacked_seq <- do.call(rbind, sequence_files)

#read the feature table, removing columns we don't need
df_feature <- fread("20241225_Wellcome_MZmine_quant.csv") %>%
  dplyr::select(-c(4, 5, 6, 8:13))

#remove features that elute before 0.7 min. 
data <- df_feature %>%  
  dplyr::filter(!`row retention time` < 0.70)

# Replace ".mzML Peak area" in all column names
colnames(data) <- gsub(".mzML Peak area", "", colnames(data))

#remove columns that has the string "Wash", "Blank_midrun", "Blank_postrun". These are washes cleaning the instrument
data <- data[, !grepl("Wash|Blank_midrun|Blank_postrun|V969", names(data), ignore.case = TRUE), with = FALSE] %>% dplyr::select(-c(2:4))

#transposing data
data_t <- transpose_and_format_df(data)

### Investigating TIC
# Investigate total peak area
data_TIC <- data.frame(TIC = rowSums(data_t %>% column_to_rownames("filename"))) %>% 
  rownames_to_column("filename") %>% left_join(stacked_seq, by = c("filename" = "File Name"))

data_TIC %>% ggscatter("Sample_order", "TIC", add = "reg.line", label = "SampleID") +
  stat_cor() + ylim(0, 1e10) 

# Samples only
data_TIC %>% dplyr::filter(!(str_detect(pattern = "Blank|six|pool|blk", filename))) %>%
  ggscatter("Sample_order", "TIC", add = "reg.line", label = "filename") + ylim(0, 1e10) +
  stat_cor() # outliers with low TIC are the Sixmix_1 or Sixmix_3 for most of the plates

# Check Blank, QCpool and QCmix
data_TIC %>% dplyr::filter(str_detect(pattern = "pool", filename)) %>% 
  ggscatter("Sample_order", "TIC", add = "reg.line", label = "filename") +
  stat_cor() + ylim(0, 2e10) #no outliers

data_TIC %>% dplyr::filter(str_detect(pattern = "Sixmix", filename)) %>% 
  ggscatter("Sample_order", "TIC", add = "reg.line", label = "filename") +
  stat_cor() + ylim(0, 6e9) #Sixmix_2 has higher TIC due to contamination/leach from previous samples

data_TIC %>% dplyr::filter(str_detect(pattern = "Blank|blk", filename)) %>% 
  ggscatter("Sample_order", "TIC", add = "reg.line", label = "Blank") +
  stat_cor() + ylim(0, 1e9) #no outliers

### Library matches
df_lib <- fread("GNPS2_wellcome_libhit2.tsv", sep = "\t")  %>%
  dplyr::select(c(1,2,4,15))  

df_lib$`#Scan#` <- as.character(df_lib$`#Scan#`)

####All library matches
# Summarize count of scans per LibraryName. Shows where the annotations of the. molecules in our data came from. 
df_lib_summary <- df_lib %>%
  group_by(LibraryName) %>%
  summarise(ScanCount = n())

# Clean up LibraryName by removing ".mgf"
df_lib_summary <- df_lib_summary %>%
  mutate(LibraryName = gsub("\\.mgf$", "", LibraryName))

# Plot
libhits <- ggplot(df_lib_summary, aes(
  x = reorder(LibraryName, -ScanCount),
  y = ScanCount
)) +
  geom_col(fill = "#FFB3D9") +                    # single pastel pink
  geom_text(aes(label = comma(ScanCount)),
            vjust = -0.3, size = 3) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(
    x = "",
    y = "",
    title = ""
  ) +
  theme_bw(base_size = 18) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 60, hjust = 1, size = 12),
    axis.text.y = element_text(size = 8),
    axis.title   = element_text(size = 10),
    plot.title   = element_text(size = 13, face = "bold"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = margin(t = 10, r = 10, b = 10, l = 80)
  )

libhits

### Blank substraction
#creating the information table with feature scan and library annotation. Shows all the features in our dataset that were annotated in GNSP2
info_feature <- df_feature %>% dplyr::select(1:4)
colnames(info_feature) <- c("Feature", "mz", "RT", "Corr_ID")
info_feature$Feature <- as.character(info_feature$Feature)
info_feature_filter <- info_feature %>% dplyr::filter(!RT < 0.70)

info_feature_complete <- info_feature_filter %>% 
  left_join(df_lib, by = c("Feature" = "#Scan#")) 

#blank substraction on the feature table. 
blank_subtraction_edited <- function(df) {
  library(dplyr)
  
  # Replace NA values with 0
  df[is.na(df)] <- 0
  
  # Identify rows where "filename" contains "Blank" or "blk"
  blank_rows <- grepl("Blank|blk", df$filename, ignore.case = TRUE)
  
  # Calculate the mean of blank rows for each feature column
  blank_means <- df %>%
    filter(blank_rows) %>%
    dplyr::select(-filename) %>%
    summarise(across(everything(), ~ mean(.x, na.rm = TRUE)))
  
  # Identify columns where the mean of non-blank rows is > 5 times the blank means
  selected_columns <- df %>%
    filter(!blank_rows) %>% # Exclude blank rows
    dplyr::select(-filename) %>%
    summarise(across(
      everything(), 
      ~ ifelse(
        is.finite(blank_means[[cur_column()]]), 
        mean(.x, na.rm = TRUE) > 5 * blank_means[[cur_column()]], 
        FALSE # Default to FALSE for invalid cases
      )
    )) %>%
    dplyr::select(where(~ all(.))) %>% # Retain only columns satisfying the condition
    colnames()
  
  # Filter the original dataframe to keep the filename and selected columns
  df_filtered <- df %>%
    dplyr::select(filename, all_of(selected_columns))
  
  return(df_filtered)
}

data_t_clean <- blank_subtraction_edited(data_t)

###merge the metadata with the transposed feature table and exporting it for later use
#rename the filename to merge with the metadata
data_t_clean$sample_name <- ifelse(
  grepl("(D[0-9]{1}[A-Z]{2}[0-9]{5}_V[0-9]+|V[0-9]{1}[A-Z]{2}[0-9]{5}_V[0-9]+|BDV[0-9]{6}_V[0-9]+|[0-9]{2}[A-Z]{2}[0-9]{5}_V[0-9]+)", data_t_clean$filename),
  sub(".*(D[0-9]{1}[A-Z]{2}[0-9]{5}_V[0-9]+|V[0-9]{1}[A-Z]{2}[0-9]{5}_V[0-9]+|BDV[0-9]{6}_V[0-9]+|[0-9]{2}[A-Z]{2}[0-9]{5}_V[0-9]+).*", "\\1", data_t_clean$filename),
  data_t_clean$filename
)

data_meta_clean <- merge(data_t_clean, qiita_metadata, by = "sample_name", all.x = TRUE)

###Removing blank samples and other mixes
data_meta_clean <- data_meta_clean %>%
  filter(!grepl("Blank|pool|blk|Sixmix", filename, ignore.case = TRUE))

#seperating visit and participant
data_meta_clean <- data_meta_clean %>%
  separate(sample_name, into = c("record_ID", "visit"), sep = "_")

#saveRDS(data_meta_clean, "clean_data_wellcome.rds")

### Counting features
df_metab <- data_meta_clean %>% dplyr::select(4:22664)

df_metab %>% summarise(across(everything(), ~ any(. > 0))) %>% 
  summarise(total = sum(across(everything())))

sum(colSums(df_metab) > 0)

