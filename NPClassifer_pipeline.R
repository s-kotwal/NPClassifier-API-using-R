# Pipeline to use NP Classifier API
# Filters out rows where Superclass is NA


# Load required libraries
library(httr)
library(jsonlite)
library(readr)
library(dplyr)

# Input and output file paths
input_file <- "/Users/shashi/Desktop/BX_Final_Screen/pos_neg_ar_final.csv"
output_file <- "/Users/shashi/Desktop/BX_Final_Screen/test_pos_neg_ar_final_NPclassified.csv"

# Load the CSV
df <- read_csv(input_file)

# Get unique, non-NA SMILES
smiles_vec <- df$SMILES %>% na.omit() %>% unique()

# Function to classify SMILES with NPClassifier
classify_smiles <- function(smiles) {
  encoded <- URLencode(smiles, reserved = TRUE)
  url <- paste0("https://npclassifier.gnps2.org/classify?smiles=", encoded)
  
  res <- tryCatch(GET(url), error = function(e) return(NULL))
  if (is.null(res) || status_code(res) != 200) {
    return(data.frame(SMILES = smiles, Class = NA, Superclass = NA, Pathway = NA, IsGlycoside = NA))
  }
  
  raw_text <- content(res, as = "text", encoding = "UTF-8")
  json_text <- sub(".*<p>(\\{.*\\})</p>.*", "\\1", raw_text)
  
  classification <- tryCatch(fromJSON(json_text), error = function(e) return(NULL))
  if (is.null(classification)) {
    return(data.frame(SMILES = smiles, Class = NA, Superclass = NA, Pathway = NA, IsGlycoside = NA))
  }
  
  class_result <- if (!is.null(classification$class_results) && length(classification$class_results) > 0)
    classification$class_results[1] else NA
  superclass <- if (!is.null(classification$superclass_results) && length(classification$superclass_results) > 0)
    classification$superclass_results[1] else NA
  pathway <- if (!is.null(classification$pathway_results) && length(classification$pathway_results) > 0)
    classification$pathway_results[1] else NA
  glyco <- if (!is.null(classification$isglycoside) && length(classification$isglycoside) > 0)
    as.character(classification$isglycoside[1]) else NA
  
  return(data.frame(
    SMILES = smiles,
    Class = class_result,
    Superclass = superclass,
    Pathway = pathway,
    IsGlycoside = glyco,
    stringsAsFactors = FALSE
  ))
}

# Run (smiles_vec for all SMILES, or (smiles_vec[1:10].. for 1st 10
results_list <- lapply(smiles_vec, function(s) {
  Sys.sleep(2)  # Be respectful to the server
  classify_smiles(s)
})
results_df <- bind_rows(results_list)

# Merge classification results with original dataframe
df_merged <- df %>%
  left_join(results_df, by = "SMILES")

# Filter out rows with NA in Superclass or Pathway (but keep NA in Class)
df_filtered <- df_merged %>%
  filter(!is.na(Superclass) & !is.na(Pathway))

# Reorder columns: place "Class", "Superclass", "Pathway", "IsGlycoside" *after* "Compound"
new_cols <- c("Class", "Superclass", "Pathway", "IsGlycoside")
compound_index <- which(names(df_filtered) == "Compound")

# Build the new column order
reordered_cols <- append(
  names(df_filtered)[!names(df_filtered) %in% new_cols],  # all other columns
  values = new_cols,
  after = compound_index  # insert after Compound
)

df_ordered <- df_filtered[, reordered_cols]

# Save output
write_csv(df_ordered, output_file)
cat("✅ Successfully saved:", output_file, "\n")
