# Set your working directory
setwd("your_directory")

# These packages are part of CRAN repository
#install.packages("data.table", dependencies = TRUE)
#install.packages("tidyverse", dependencies = TRUE)
#install.packages("pheatmap", dependencies = TRUE)

## Load the packages required for the analysis
library(data.table)
library(tidyverse)
library(pheatmap)

# Specify the folder path - it should be the folder inside the working directory 
folder_path <- "the_FASST_search_directory"

## Download/Import the ReDU metadata file - it should be in the working directory folder and NOT be in the sub-folder with the csv files from the Fast Search

# Define the filename for the ReDU metadata
processed_redu_metadata <- "all_sampleinformation.tsv"

# Check if the pre-processed metadata file exists in the working directory
if (!file.exists(file.path(getwd(), processed_redu_metadata))) {
  redu_url <- "https://redu.gnps2.org/dump"
  options(timeout = 600) # If 10 min is not enough, add more time 
  download.file(redu_url, file.path(getwd(), processed_redu_metadata), mode = "wb")
  redu_metadata <- data.table::fread(processed_redu_metadata)
} else {
  redu_metadata <- data.table::fread(processed_redu_metadata)
}

# Optional: In lieu of the previous fread() command, if memory issues occur or the program crashes, we recommend commenting out the previous command, uncommenting the line with the read_tsv() command, and following the instructions below
## Within the col_select parameter, we recommend only specifying the columns needed for a given analysis to minimize problems with memory limitations; the filename and NCBITaxonomy columns are likely to be always used
### If memory issues persist, we recommend reading in the metadata in chunks, performing the analysis desired on each chunk, and appropriately recombining the results at the end
#### redu_metadata <- readr::read_tsv("all_sampleinformation.tsv", col_select = c("filename", "NCBITaxonomy", "UBERONBodyPartName", "DOIDCommonName", "HealthStatus", "BiologicalSex", "LifeStage"), show_col_types = FALSE)

# Get the list of all .csv files in the folder
file_list <- list.files(folder_path, pattern = "*.csv", full.names = TRUE)

# Read each .csv file and add the Compound column in each df
df_list <- lapply(file_list, function(file) {
  df <- read_csv(file)
  df$Compound <- tools::file_path_sans_ext(basename(file))
  return(df)
})

# Combine all dfs into a single df
molecules_interest <- bind_rows(df_list)

molecules_interest_filtered <- molecules_interest |> 
  dplyr::filter(`Delta Mass` >= -0.05 & `Delta Mass` <= 0.05)

# Prepare the data tables for merging
## Create a function to extract the desired segment from the USI column

MassiveID_filename <- function(USI) {
  USI <- gsub("/", ":", USI)
  USI <- sub("\\.[^\\.]*$", "", USI)
  parts <- unlist(strsplit(USI, ":"))
  combined <- paste(parts[2], parts[length(parts)], sep = ":")
  return(combined)
}

# Apply the function to each row of the USI column in the molecules_interest
molecules_interest_filtered$USI <- vapply(molecules_interest_filtered$USI, MassiveID_filename, FUN.VALUE = character(1))

# Prepare the ReDU metadata USI column for merging with FASST output table
## Create a function to extract the datasetID and the last segment (filename)

ReDU_USI <- function(USI) {
  USI <- gsub("/", ":", USI)
  USI <- sub("\\.[^\\.]*$", "", USI)
  parts <- unlist(strsplit(USI, ":"))
  combined <- paste(parts[2], parts[length(parts)], sep = ":")
  return(combined)
}

# Apply the function to each row of the fxlename column in the ReDU output table
redu_metadata$USI <- vapply(redu_metadata$USI, ReDU_USI, FUN.VALUE = character(1))

# Merge the ReDU metadata table and the FASST MASST output table
ReDU_MASST <- left_join(molecules_interest_filtered, redu_metadata, by = "USI", relationship = "many-to-many")

# Once both data tables are merged, ones can filter the table which based on the research question
## To note: not all publicly available files have associated metadata and we strongly encourage scientists to make 
### their data available with a very detailed metadata (sample information)
#### As more data are being deposited in repositories more matches will be uncovered and more results will be embedded in heatmaps 

# Standardize the body parts and Health Status
ReDU_MASST_standardize <- ReDU_MASST |> 
  dplyr::mutate(
    UBERONBodyPartName = str_replace_all(UBERONBodyPartName, 'skin of trunk|skin of pes|head or neck skin|axilla skin|skin of manus|arm skin|skin of leg', 'skin'),
    UBERONBodyPartName = str_replace_all(UBERONBodyPartName, 'blood plasma|blood serum', 'blood'),
    HealthStatus = str_replace(HealthStatus, 'Chronic Illness', 'chronic illness'),
    HealthStatus = str_replace(HealthStatus, 'Healthy', 'healthy')
  )

# Separate humans and rodents from the merged data table
df_humans <- ReDU_MASST_standardize |>  
  dplyr::filter(NCBITaxonomy == "9606|Homo sapiens")

# Define a list for rodents taxonomy IDs
list_rattus_mus <- c('10088|Mus', '10090|Mus musculus', '10105|Mus minutoides', '10114|Rattus', '10116|Rattus norvegicus')

# Separate rodents
df_rodents <- ReDU_MASST_standardize |>  
  dplyr::filter(NCBITaxonomy %in% list_rattus_mus)

analyze_counts <- function(df, column_interest) {
  
# Create a list of all unique entries in the column of interest and create a df
  df_body_parts <- df |>  distinct(across(all_of(column_interest)))
  
  # Count occurrences of each entry in the column of interest
  df_BodyPartName_counts <- df |> 
    count(across(all_of(column_interest)), name = "Counts_fastMASST")
  
  # Aggregate the number and list of unique Compounds for each entry
  compounds <- df |> 
    group_by(across(all_of(column_interest))) |> 
    summarise(Compounds = n_distinct(Compound),
              CompoundsList = toString(unique(Compound))) |> 
    ungroup()
  
  # Merge all the data into a single data frame
  combined <- df_body_parts |> 
    left_join(df_BodyPartName_counts, by = column_interest) |> 
    left_join(compounds, by = column_interest)
  
  return(combined)
}

# Get a glimpse of the number of counts per organ
body_counts_humans <- analyze_counts(df_humans, "UBERONBodyPartName")
head(body_counts_humans)

body_counts_rodents <- analyze_counts(df_rodents, "UBERONBodyPartName")
head(body_counts_rodents)

# Create a function to pivot the table for data visualization
prepare_pivot_table <- function(df, column_interest, compound) {
  
  grouped_df <- df |> 
    group_by(across(all_of(c(compound, column_interest)))) |> 
    summarise(Count = n(), .groups = 'drop')
  
  pivot_table <- grouped_df |> 
    pivot_wider(names_from = all_of(compound), values_from = Count, values_fill = list(Count = 0))
  
  return(pivot_table)
}

# Define the variables based on your research question 
## Here we are interesting in organ distribution in humans and rodents of the molecule of interest
variable <- 'UBERONBodyPartName'
pivot_table_humans <- prepare_pivot_table(df_humans, variable, 'Compound')
pivot_table_rodents <- prepare_pivot_table(df_rodents, variable, 'Compound')

# Prepare the table to be compatible with pheatmap package
humans_molecules_counts_by_bodypart <- pivot_table_humans |> 
  dplyr::arrange(UBERONBodyPartName) |> 
  tibble::column_to_rownames("UBERONBodyPartName")
# Prepare the table to be compatible with pheatmap package
rodents_molecules_counts_by_bodypart <- pivot_table_rodents |>
  dplyr::arrange(UBERONBodyPartName) |> 
  tibble::column_to_rownames("UBERONBodyPartName")

# Convert all columns to numeric for the humans df
humans_molecules_counts_by_bodypart <- humans_molecules_counts_by_bodypart |> 
  dplyr::mutate(across(everything(), as.numeric))
# Convert all columns to numeric for the rodents df
rodents_molecules_counts_by_bodypart <- rodents_molecules_counts_by_bodypart |> 
  dplyr::mutate(across(everything(), as.numeric))

# Define your chosen colors
colors_version <- c("#FFFFFF", "#C7D6F0", "#EBB0A6")
# Creating the gradient function
color_gradient <- colorRampPalette(colors_version)
# Generate 30 discrete colors from this gradient
gradient_colors <- color_gradient(30)

# The users can log scale or not the data
log_humans_molecules_counts_by_bodypart <- log2(1 + humans_molecules_counts_by_bodypart)
#write.csv(log_humans_molecules_counts_by_bodypart, 
#          file = "log_humans_molecules_counts_by_bodypart.csv", 
#          row.names = TRUE)
log_rodents_molecules_counts_by_bodypart <- log2(1 + rodents_molecules_counts_by_bodypart)
#write.csv(log_rodents_molecules_counts_by_bodypart, 
#         file = "log_rodents_molecules_counts_by_bodypart.csv", 
#         row.names = TRUE)
# Organ distribution in humans
## Use heatmap for data visualization or organ distribution - humans
### If one MS/MS spectrum is used in reverse metabolomics, the cluster_rows and cluster_cols should be set to FALSE
Organ_humans <- pheatmap(log_humans_molecules_counts_by_bodypart,
                         color = gradient_colors,
                         cluster_rows = FALSE,
                         cluster_cols = TRUE,
                         angle_col = 90,
                         main = "Organ distribution in humans",
                         fontsize = 10,
                         cellwidth = 15,
                         cellheight = 15,
                         treeheight_row = 100,
                         fontsize_row = 12,
                         fontsize_col = 12,
                         legend_fontsize = 10,
                         border_color = NA)
Organ_humans
ggsave("Organ_distribution_in_humans.pdf", plot = Organ_humans, width = 10, height = 10, dpi = 900)
getwd()

# Organ distribution in rodents
## Use heatmap for data visualization or organ distribution - rodents
### If one MS/MS spectrum is used in reverse metabolomics, the cluster_rows and cluster_cols should be set to FALSE
Organ_rodents <- pheatmap(log_rodents_molecules_counts_by_bodypart,
         color = gradient_colors,
         cluster_rows = FALSE,
         cluster_cols = TRUE,
         angle_col = 90,
         main = "Organ distribution in rodents",
         fontsize = 10,
         cellwidth = 15,
         cellheight = 15,
         treeheight_row = 100,
         fontsize_row = 12,
         fontsize_col = 12,
         legend_fontsize = 10,
         border_color = NA)
Organ_rodents
ggsave("Organ_distribution_in_rodents.pdf", plot = Organ_rodents, width = 10, height = 10, dpi = 900)
getwd()


# Health phenotype association
## Filter for human information in the ReDU metadata
df_redu_humans <- redu_metadata |>  
  dplyr::filter(NCBITaxonomy == "9606|Homo sapiens")

# Filter for rodent information in the ReDU metadata
df_redu_rodents <- redu_metadata |>  
  dplyr::filter(NCBITaxonomy %in% list_rattus_mus)

# Humans - filter for lifestage, DOIDCommonName, HealthStatus, BiologicalSex
## For LifeStage
human_ReDU_LifeStage <- df_redu_humans |> 
  dplyr::count(LifeStage) |> 
  dplyr::rename(LifeStage_counts = n, LifeStage = LifeStage)
human_ReDU_LifeStage$LifeStage_counts <- as.numeric(human_ReDU_LifeStage$LifeStage_counts)
# For DOIDCommonName
human_ReDU_DOIDCommonName <- df_redu_humans |> 
  dplyr::count(DOIDCommonName) |> 
  dplyr::rename(DOIDCommonName_counts = n, DOIDCommonName = DOIDCommonName)
human_ReDU_DOIDCommonName$DOIDCommonName_counts <- as.numeric(human_ReDU_DOIDCommonName$DOIDCommonName_counts)
  # For HealthStatus
human_ReDU_HealthStatus <- df_redu_humans |> 
  dplyr::count(HealthStatus) |> 
  dplyr::rename(HealthStatus_counts = n, HealthStatus = HealthStatus)
human_ReDU_HealthStatus$HealthStatus_counts <- as.numeric(human_ReDU_HealthStatus$HealthStatus_counts)
# For BiologicalSex
human_ReDU_BiologicalSex <- df_redu_humans |> 
  dplyr::count(BiologicalSex) |> 
  dplyr::rename(BiologicalSex_counts = n, BiologicalSex = BiologicalSex)
human_ReDU_BiologicalSex$BiologicalSex_counts <- as.numeric(human_ReDU_BiologicalSex$BiologicalSex_counts)

# Rodents - filter for lifestage, DOIDCommonName, HealthStatus, BiologicalSex
## For LifeStage
Rodents_ReDU_LifeStage <- df_redu_rodents |> 
  dplyr::count(LifeStage) |> 
  dplyr::rename(LifeStage_counts = n, LifeStage = LifeStage)
Rodents_ReDU_LifeStage$LifeStage_counts <- as.numeric(Rodents_ReDU_LifeStage$LifeStage_counts)
# For DOIDCommonName
Rodents_ReDU_DOIDCommonName <- df_redu_rodents |> 
  dplyr::count(DOIDCommonName) |> 
  dplyr::rename(DOIDCommonName_counts = n, DOIDCommonName = DOIDCommonName)
Rodents_ReDU_DOIDCommonName$DOIDCommonName_counts <- as.numeric(Rodents_ReDU_DOIDCommonName$DOIDCommonName_counts)
# For HealthStatus
Rodents_ReDU_HealthStatus <- df_redu_rodents |> 
  dplyr::count(HealthStatus) |> 
  dplyr::rename(HealthStatus_counts = n, HealthStatus = HealthStatus)
Rodents_ReDU_HealthStatus$HealthStatus_counts <- as.numeric(Rodents_ReDU_HealthStatus$HealthStatus_counts)
# For BiologicalSex
Rodents_ReDU_BiologicalSex <- df_redu_rodents |> 
  dplyr::count(BiologicalSex) |> 
  dplyr::rename(BiologicalSex_counts = n, BiologicalSex = BiologicalSex)
Rodents_ReDU_BiologicalSex$BiologicalSex_counts <- as.numeric(Rodents_ReDU_BiologicalSex$BiologicalSex_counts)


# Normalization of the FASST search output 
## Grouping and counting occurrences
### (Optional) remove the # symbol to activate the line and apply a minimum of 3 counts to be include before normalization 
#### Here we are interested in finding disease association
grouped_df_humans <- df_humans |>
  group_by(Compound, DOIDCommonName) |>  
  summarise(Count = n()) |> 
  ungroup() #|> 
  #dplyr::filter(Count >= 3)

grouped_df_humans_pivot_table <- grouped_df_humans |>
  pivot_wider(names_from = Compound, values_from = Count, values_fill = list(Count = 0))

# Merging 
merged_DOID_humans <- left_join(grouped_df_humans_pivot_table, human_ReDU_DOIDCommonName, by = "DOIDCommonName")
merged_DOID_humans$DOIDCommonName <- gsub("Crohn's disease", "crohn's disease", merged_DOID_humans$DOIDCommonName)

# Normalizing counts by the number of files available for each DOIDCommonName
columns_to_normalize <- setdiff(names(merged_DOID_humans), c("DOIDCommonName", "DOIDCommonName_counts"))
  
normalized_merged_DOID_humans <- merged_DOID_humans |>
  dplyr::filter(!DOIDCommonName == "missing value") |> 
  dplyr::mutate(across(all_of(columns_to_normalize), ~ .x / .data$DOIDCommonName_counts)) |> 
  dplyr::select(-DOIDCommonName_counts)
  
# Adding a sum row with correct names
sums <- colSums(dplyr::select(normalized_merged_DOID_humans, where(is.numeric)), na.rm = TRUE)
sums_df <- as.data.frame(t(sums))
sums_df$DOIDCommonName <- 'Sum'
sums_df <- sums_df[, names(normalized_merged_DOID_humans)]
merged_sum_humans_DOID <- bind_rows(normalized_merged_DOID_humans, sums_df)
merged_sum_humans_DOID <- merged_sum_humans_DOID |> 
  dplyr::filter(!is.na(DOIDCommonName)) |>
  dplyr::mutate(across(where(is.numeric), ~replace_na(.x, 0)))
  
# Normalizing each numeric column by the value in the 'Sum' row to get percentages
merged_sum_humans_DOID_percentage <- merged_sum_humans_DOID |> 
  dplyr::mutate(across(all_of(columns_to_normalize), ~ .x / .x[n()] * 100))
  
# Remove the 'Sum' row
merged_sum_humans_DOID_percentage_plot <- merged_sum_humans_DOID_percentage |>
  dplyr::filter(DOIDCommonName != "Sum") |> 
  dplyr::arrange(DOIDCommonName) |> 
  tibble::column_to_rownames("DOIDCommonName")

#write.csv(merged_sum_humans_DOID_percentage_plot, 
         #file = "merged_sum_humans_DOID_percentage_plot.csv", 
         #row.names = TRUE)

# Get health phenotype information
## If one MS/MS spectrum is used in reverse metabolomics, the cluster_rows and cluster_cols should be set to FALSE
Diseases_humans <-pheatmap(merged_sum_humans_DOID_percentage_plot,
         color = gradient_colors,
         cluster_rows = FALSE,
         cluster_cols = TRUE,
         angle_col = 90,
         main = "Health phenotype association",
         fontsize = 10,
         cellwidth = 15,
         cellheight = 15,
         treeheight_row = 100,
         fontsize_row = 12,
         fontsize_col = 12,
         legend_fontsize = 10,
         border_color = NA)
Diseases_humans

# Export the data 
ggsave("Diseases_humans.pdf", plot = Diseases_humans, width = 10, height = 10, dpi = 900)
getwd()


