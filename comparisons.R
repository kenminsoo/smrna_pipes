library("readr")
library("BiocManager") #cite
library("DESeq2") #cite
library("ggplot2")
library(dplyr)
library(tibble)

#compare generated files with the base truth

base_truthdir <- "/home/kenmn/nfs_scratch/simulated_data_1/sim_counts_matrix.rda" # nolint
experimental_data_dir <- "/home/kenmn/nfs_scratch/sim_results/"

featurecounts_comparison <- function(true_dir, prefix, pipeline, exp_data_dir) {

load(true_dir)
#loads in as counts_matrix

pipe_data <- read_delim(paste(experimental_data_dir,pipeline,"/",prefix,"/",prefix,"-",pipeline,"-","counts.tsv",sep=""), delim = "\t",escape_double = FALSE, trim_ws = TRUE, skip = 1)

#the seqs might have additional source info, remove
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
fp_combined <- fp_combined[!complete.cases(fp_combined),]
fp_combined <- fp_combined[ rowSums(fp_combined[17:26]) > 0,]

write.csv(fp_combined, paste(prefix,pipeline,"fp.csv",sep="-"))
write.csv(right_diff_per, paste(prefix,pipeline,'per-error.csv',sep="-"))
write.csv(right_diff, paste(prefix,pipeline,'raw-error.csv',sep="-"))
write.csv(combined_data, paste(prefix,pipeline,'raw-counts.csv',sep="-"))

}
#we will analyze the csv files at a later date