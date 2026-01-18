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
library(RColorBrewer)
library(taxize)
library(purrr)
library(ape)
library(readr)
library(rlang)

### Loading data
#MicrobeMASST
df <- read_tsv("/ENTER YOUR PATH HERE/microbe_masst/output/data/Batch_MASST_microbe.txt")
df$Scan <- as.character(df$Scan)
df %>% summarise(n_unique = n_distinct(file_usi))
df %>% summarise(n_unique = n_distinct(Dataset))

#metabolite info
metab1 <- fread("20241225_Wellcome_MZmine_quant.csv", select = 1:13)
metab1 <-rename(metab1, Scan = `row ID`)
metab1$Scan <- as.character(metab1$Scan)

#Features from CR
model_metab <- readRDS("ranked_metabolites.rds")
model_metab <-rename(model_metab, Scan = Metabolite)

#merging CR model with metabolite info
metab <- model_metab %>%
  left_join(metab1, by = "Scan")       

#annotations
metabanno  <- fread("annotated_molecules_edit.csv")
metabanno <-rename(metabanno, Scan = scanID)
metabanno$Scan <- as.character(metabanno$Scan)

#Ranking taxons for microbes
rank <- fread("/ENTER YOUR PATH HERE/trees/microbe_masst_tree/list_ncbi_microbemasst.csv")

microbe_unranked <- df %>%
  left_join(metab, by = "Scan")       
microbe_unranked <-rename(microbe_unranked, name = Taxaname_file)

microbe1 <- microbe_unranked %>%
  left_join(rank, by = "name")

### Cleaning 
#Removing features also found in human (25705 matches before)
microbe <- microbe1 %>%
  group_by(Scan) %>%
  filter(!any(name == "Homo sapiens")) %>%
  ungroup()

### Preparing for taxonomy assignment
# Ensure NCBI IDs are numeric
microbe$Taxa_NCBI <- as.numeric(microbe$Taxa_NCBI)

# Get unique IDs (to avoid redundant API calls)
unique_ids <- unique(microbe$Taxa_NCBI)
length(unique_ids)

#Keep only unique microbes
df_unique <- microbe %>%
  distinct(Taxa_NCBI, .keep_all = TRUE) %>%  
  select(Taxa_NCBI, name)

# Filter valid entries (removes one)
df_unique <- df_unique[df_unique$Taxa_NCBI != "missing value" & !is.na(df_unique$Taxa_NCBI), ]
taxa_id <- as.numeric(df_unique$Taxa_NCBI)
names(taxa_id) <- df_unique$name

# save_file and final_file contain information extracted from NCBI for each batch
save_file <- "ncbi_classifications_species_partial.rds"   # intermediate save
final_file <- "ncbi_classifications_species_final.rds"

# batch size for processing samples - make it small to reduce errors
batch_size <- 5         

# number of retries for failed chunks
max_retries <- 3  

# Split into chunks and batch size
taxa_chunks <- split(taxa_id, ceiling(seq_along(taxa_id) / batch_size))

# total_chunks represents the total number of batches
total_chunks <- length(taxa_chunks)

# 3️⃣ Set your NCBI API key (for higher rate limits)
Sys.setenv(ENTREZ_KEY = "Enter key here") # <- change the part between "" to your key


### Loop to load taxonomy
# Load previous progress if it exists (e.g. some batches were already processed)
if (file.exists(save_file)) {
  message("Resuming from saved file...")
  # this line reads the RDS file (the type of file that save_file is) back into R
  # so that you can pick up where you left off
  all_results <- readRDS(save_file)
} else {
  # this line of code lists the number of batches so that you know it starts at beginning
  all_results <- vector("list", total_chunks)
}

# Function to process one batch safely
# taxa_chunks is chunk and i is each value or batch in it
process_chunk <- function(chunk, i) {
  message("Processing chunk ", i, " of ", total_chunks, " (", length(chunk), " taxa)")
  for (attempt in 1:max_retries) {
    # runs NCBI query to get information for each batch
    result <- tryCatch(
      classification(chunk, db = "ncbi"),
      error = function(e) {
        message("  Error in chunk ", i, " attempt ", attempt, ": ", e$message)
        return(NULL)
      }
    )
    if (!is.null(result)) {
      message("  Chunk ", i, " completed successfully")
      return(result)
    }
    Sys.sleep(2)  #  delay time before retry
  }
  message("  Skipping chunk ", i, " after ", max_retries, " attempts")
  return(NULL)
}

# Main loop for batches
for (i in seq_len(total_chunks)) {
  if (!is.null(all_results[[i]])) {
    message("Skipping chunk ", i, " (already done)")
    next
  }
  all_results[[i]] <- process_chunk(taxa_chunks[[i]], i)
  saveRDS(all_results, save_file)
  Sys.sleep(1)
}

# Retry failed batches at the end
skipped <- which(sapply(all_results, is.null))
if (length(skipped) > 0) {
  message("Retrying ", length(skipped), " skipped chunks...")
  for (i in skipped) {
    all_results[[i]] <- process_chunk(taxa_chunks[[i]], i)
    saveRDS(all_results, save_file)
    Sys.sleep(1)
  }
}

# Finalize results
all_results_combined <- do.call(c, all_results)
saveRDS(all_results_combined, final_file)
message("All results saved as ", final_file)

# Clean results & build tree
NCBIClass_clean <- all_results_combined[!sapply(all_results_combined, function(x) {
  is.null(x) || nrow(x) == 0
})]

# Remove duplicates
NCBIClass_nodup <- NCBIClass_clean[!duplicated(names(NCBIClass_clean))]

# Standard taxonomic ranks needed for all species
standard_ranks <- c("domain", "kingdom", "phylum", "class", 
                    "order", "family", "genus", "species")

# Keep only these ranks for all speices (some species contain others that will be removed for so that the tree can be created)
NCBIClass_filtered <- lapply(NCBIClass_nodup, function(df) {
  df[df$rank %in% standard_ranks, , drop = FALSE]
})

# Keep lineages with all standard ranks
NCBIClass_complete <- NCBIClass_filtered[sapply(NCBIClass_filtered, function(x) {
  all(standard_ranks %in% x$rank)
})]

cat("Original taxa:", length(NCBIClass_nodup), "\n")
cat("With all standard ranks:", length(NCBIClass_complete), "\n")

### Merge and export
rank <- NCBIClass_complete
microbe_df %>% microbe
select(1:35)

# Simplified get_lineage: returns only rank + name
get_lineage_names <- function(taxon_id, rank_list) {
  id_chr <- as.character(taxon_id)
  
  # 1. Direct match
  if (!is.null(rank_list[[id_chr]])) {
    lineage <- rank_list[[id_chr]]
  } else {
    # 2. Look downstream in child lineages
    found <- FALSE
    for (child in rank_list) {
      if (taxon_id %in% child$id) {
        cutoff <- which(child$id == taxon_id)
        lineage <- child[seq_len(cutoff), ]
        found <- TRUE
        break
      }
    }
    if (!found) return(tibble(rank = NA, name = NA))
  }
  
  # Keep only rank + name
  lineage %>% select(rank, name)
}

# Build wide dataframe with only taxon names
df_wide <- df %>%
  mutate(rank_info = purrr::map(Taxa_NCBI, get_lineage_names, rank)) %>%
  mutate(rank_info = purrr::map(rank_info, ~{
    pivot_wider(.x,
                names_from = rank,
                values_from = name,
                values_fn = \(x) paste(unique(x), collapse = "; "))
  })) %>%
  unnest(rank_info)

#write.csv(df_wide, "microbes-CR-taxonomy.csv", row.names = FALSE)




