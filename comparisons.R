library("readr")
library("BiocManager") #cite
library("DESeq2") #cite
library("ggplot2")
library('openxlsx') 
library('EnhancedVolcano')
library('patchwork')
library(dplyr)
library(tibble)

#compare generated files with the base truth

base_truthdir <- "/home/kenmn/nfs_scratch/simulated_data_1/sim_counts_matrix.rda"

load(base_truthdir)
#loads in as counts_matrix

pipe_data <- read_delim('/home/kenmn/nfs_scratch/simulated_data_1/counts_uni.tsv', delim = "\t",
escape_double = FALSE, trim_ws = TRUE, skip = 1)

#the seqs have additional source info, remove
pipe_data$Geneid <- gsub("\\._.*","",pipe_data$Geneid)

counts_matrix_df <- as.data.frame(counts_matrix)

counts_rowcol <- tibble::rownames_to_column(counts_matrix_df, "Geneid")

#comparing the right answers
combined_data <- inner_join(counts_rowcol, pipe_data, by = 'Geneid')

n_base <- 2
n_exp <- 17

right_diff <- data.frame(matrix(ncol = 11, nrow = 1400))
right_diff[1] <- combined_data['Geneid']

for (i in 1:10) {
    right_diff[i+1] <- combined_data[n_exp] - combined_data[n_base]
    n_base <- n_base + 1
    n_exp <- n_exp + 1
}

n_base <- 2
n_exp <- 17

right_diff_per <- data.frame(matrix(ncol = 11, nrow = 1400))
right_diff_per[1] <- combined_data['Geneid']
for (i in 1:10) {
    right_diff_per[i+1] <- (combined_data[n_exp] - combined_data[n_base])/combined_data[n_base]
    n_base <- n_base + 1
    n_exp <- n_exp + 1
}


#false positives
fp_combined <- right_join(counts_rowcol, pipe_data, by = 'Geneid')