#!/usr/bin/env Rscript
# Load libraries
library(plyr)
library(tidyverse)
library(ggrepel)
library(svglite)

PIPELINEDIR <- "MY_PIPELINE_DIR" # Replace with your actual pipeline directory

args = commandArgs(trailingOnly=TRUE)
# test if there is at least one argument: if not, return an error
if (length(args)!=3) {
  stop("Three argument must be provided", call.=FALSE)
} else {
  PCA_Folder <- args[1]
  sample_info_file<- args[2]
  vcf_file <- args[3]
}

# Clean duplicated sample names like "sample_sample" into "sample"
clean_sample_name <- function(id) {
  ifelse(grepl("^(.*)_\\1$", id), sub("^(.*)_\\1$", "\\1", id), id)
}

#Test code ajout échantillon#

# Charger la ligne d'en-tête contenant les noms d'échantillons depuis le VCF
#vcf_file <- "path/to/your_file.vcf" 
header_line <- readLines(vcf_file, n = 1000)
header_line <- header_line[grepl("^#CHROM", header_line)]

# Extraire les noms d'échantillons
sample_ids <- strsplit(header_line, "\t")[[1]][10:length(strsplit(header_line, "\t")[[1]])] # extrait de la 10e (commencement des sampleID) colonne à la fin

# Créer un dataframe pour ces nouveaux échantillons
new_samples <- data.frame(
  Sample = sample_ids,
  Family.ID = sample_ids,
  Population = "cohort",
  stringsAsFactors = FALSE
)

# Ajouter à kg uniquement ceux qui n’y sont pas déjà
#samples_to_add <- subset(new_samples, !(Sample %in% kg$Sample))
samples_to_add <- new_samples

# Ajouter et reconstruire SuperPopulation pour les nouveaux
#samples_to_add$SuperPopulation <- "cohort"

# Ajouter à `kg`


# Vérification (optionnel)
cat("Nombre de nouveaux échantillons ajoutés :", nrow(samples_to_add), "\n")

#===========================================OG Code===========================================#
# Read sample info file
sample_info_file<-file.path(PIPELINEDIR, "Software/Code/Ethnic_PCA/sample_info_file.csv")
kg <- read.csv(sample_info_file, header = TRUE, stringsAsFactors = FALSE, quote = "***", fill = TRUE)

print("Montre moi kg avant ajout")
print(head(kg))
print("Montre moi sample_to_add avant ajout")
print(head(samples_to_add))

kg <- rbind.fill(kg, samples_to_add)

print("Montre moi kg après ajout")
print(head(kg))
print("Montre moi kg cohort après ajout")
print(subset(kg, Population == "cohort"))
# Assign super-population labels
kg$SuperPopulation <- NA
kg$SuperPopulation[kg$Population %in% c("CHB", "JPT", "CHS", "CDX", "KHV")] <- "EAS"
kg$SuperPopulation[kg$Population %in% c("CEU", "TSI", "FIN", "GBR", "IBS")] <- "EUR"
kg$SuperPopulation[kg$Population %in% c("YRI", "ASW", "ACB", "LWK", "GWD", "MSL", "ESN")] <- "AFR"
kg$SuperPopulation[kg$Population %in% c("MXL", "PUR", "CLM", "PEL")] <- "AMR"
kg$SuperPopulation[kg$Population %in% c("GIH", "PJL", "BEB", "STU", "ITU")] <- "SAS"
kg$SuperPopulation[kg$Population == "cohort"] <- "cohort"

# Read PCA scores
pc1 <- read.table(file.path(PCA_Folder, "1kg-PC1-score.sscore"), header = TRUE, stringsAsFactors = FALSE, comment.char = "", fill = TRUE)
pc2 <- read.table(file.path(PCA_Folder, "1kg-PC2-score.sscore"), header = TRUE, stringsAsFactors = FALSE, comment.char = "", fill = TRUE)
pc3 <- read.table(file.path(PCA_Folder, "1kg-PC3-score.sscore"), header = TRUE, stringsAsFactors = FALSE, comment.char = "", fill = TRUE)


# Clean sample names consistently
kg$Sample <- clean_sample_name(kg$Sample)
pc1$IID <- clean_sample_name(pc1$IID)
pc2$IID <- clean_sample_name(pc2$IID)
pc3$IID <- clean_sample_name(pc3$IID)

# Filter to matching samples
samp <- kg$Sample
IID <- pc1$IID
not_in_pc <- setdiff(samp, IID)
not_in_kg <- setdiff(IID, samp)

pc1 <- pc1[!(pc1$IID %in% not_in_kg), ]
pc2 <- pc2[!(pc2$IID %in% not_in_kg), ]
pc3 <- pc3[!(pc3$IID %in% not_in_kg), ]
kg <- kg[!(kg$Sample %in% not_in_pc), ]

# Prepare for merging
kg$IID <- kg$Sample

# Merge data
# Print the header of pc1


# Extract relevant columns
df <- pc1[, c("X.FID", "IID", "SCORE1_AVG")]
df <- merge(df, pc2[, c("IID", "SCORE1_AVG")], by = "IID")
df <- merge(df, pc3[, c("IID", "SCORE1_AVG")], by = "IID")
df <- merge(df, kg[, c("IID", "SuperPopulation")], by = "IID")

# Rename columns
names(df)[names(df) == "SCORE1_AVG.x"] <- "PC1"
names(df)[names(df) == "SCORE1_AVG.y"] <- "PC2"
names(df)[names(df) == "SCORE1_AVG"]   <- "PC3"

# Clean sample names
df$CleanName <- clean_sample_name(df$IID)

# Subset only cohort samples
cohort_samples <- subset(df, SuperPopulation == "cohort")

# Define color palette (for consistency)
palette <- c("purple", "orange", "green", "blue", "red")


#Print de l'enfer du débuggage
print("Montre moi kg")
print(head(kg))
print("Montre moi kg cohort")
print(subset(kg, SuperPopulation == "cohort"))
print("Montre moi pc1")
print(head(pc1))
print("Montre moi df")
print(head(df))
print("Montre moi cohort_samples")
print(cohort_samples)


# Loop over each cohort sample and generate both plots
for (i in 1:nrow(cohort_samples)) {
  
  sample_row <- cohort_samples[i, ]
  clean_name <- make.names(sample_row$CleanName)  # sanitize for file name
  
  # Create output folder for this sample if it doesn't exist
  sample_results_dir <- file.path(PIPELINEDIR,"Results", clean_name)
  if (!dir.exists(sample_results_dir)) {
    dir.create(sample_results_dir, recursive = TRUE)
  }
  # Ajouter le nom de l'échantillon à SampleList.txt
  sample_list_file <- file.path(PIPELINEDIR,"Results", "SampleList.txt")
  write(clean_name, file = sample_list_file, append = TRUE)

  # --- PC1 vs PC2 ---
  p1 <- ggplot(data = subset(df, SuperPopulation != "cohort"),
               aes(x = PC1, y = PC2, color = SuperPopulation,)) +
    geom_point(size = 3) +
    scale_color_manual(values = palette) +
    geom_point(data = sample_row, aes(x = PC1, y = PC2), color = "black", size = 4) +
    geom_text_repel(data = sample_row, aes(x = PC1, y = PC2, label = CleanName),
                    size = 4, color = "black", force = 10, min.segment.length = 0) +
    ggtitle(paste("PCA PC1 vs PC2 - Sample:", sample_row$CleanName))
  
  
  ggsave(filename = file.path(sample_results_dir, paste0(clean_name, "_PC1_PC2.svg")), plot = p1, width = 8, height = 6, device = "svg")
  
  Sys.sleep(0.5)  # pause between plots

  # --- PC2 vs PC3 ---
  p2 <- ggplot(data = subset(df, SuperPopulation != "cohort"),
               aes(x = PC2, y = PC3, color = SuperPopulation)) +
    geom_point(size = 3) +
    scale_color_manual(values = palette) +
    geom_point(data = sample_row, aes(x = PC2, y = PC3), color = "black", size = 4) +
    geom_text_repel(data = sample_row, aes(x = PC2, y = PC3, label = CleanName),
                    size = 4, color = "black", force = 10, min.segment.length = 0) +
    ggtitle(paste("PCA PC2 vs PC3 - Sample:", sample_row$CleanName))
  
  
  ggsave(filename = file.path(sample_results_dir, paste0(clean_name, "_PC2_PC3.svg")), plot = p2, width = 8, height = 6, device = "svg")
  
  Sys.sleep(0.5)
}
