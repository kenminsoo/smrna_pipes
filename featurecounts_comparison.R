#This script will compare featurecounts outputs
#with the ground truth.

#expected input will be a csv file:
#Geneid, sample1, ..., samplen, chr, start, end, strand, length, exp_sample1, ..., exp_samplen # nolint

#For now we will produce the following outputs:

#PCA of the samples from all pipes.
#Linear regression of ground truth with experimental.
#Sensitivity = TP / TP + FN # nolint

#if read found, it will be counted as TP 
#if read not found, it will be counted as FN (i.e. error >95%)

#Precision = TP / FP + TP # nolint

#if read is found in experimental, it will be counted as FP

#F-measure = 2 Precision X Sensitivity / Precision + Sensitivity  # nolint
#Jaccardi Similarity index
#Pearson correlation matrix with each ncRNA type

#============#

#import packages
library(ggplot2)
library(stringr)
library(dplyr)
library(reshape2)
library(ggpmisc)
library(plyr)

#============#
#import data, assuming running locally on computer
data_dir <- "/Users/kenminsoo/Desktop/benchmark_data/"
trial <- "base"
#for now, manually put in names
pipes <- list("aasra", "compsra", "manual")
data <- list()
counts <- "raw-counts.csv"
per <- "per-error.csv"
fp <- "fp.csv"
samples <- 10

#1 = raw counts
#2 = percent error
#3 = false positives
n <- 1

separated <- list()

bio_types <- list("5S", "rRNA", "miR", "piR", "tRNA")

#import data from specified dir
for (pipe in pipes){
    pipe_dir <- paste(data_dir, trial, "/", pipe, "/",sep = "")
    base_data_name <- paste(trial, "-", pipe, "-", sep = "")
    data[[pipe]] <- list()

    data[[pipe]][[1]] <- read.csv(paste(pipe_dir, 
    base_data_name, counts, sep = "")) 

    data[[pipe]][[2]] <- read.csv(paste(pipe_dir, 
    base_data_name, per, sep = "")) 

    data[[pipe]][[3]] <- read.csv(paste(pipe_dir, base_data_name, fp, sep = ""))
    n <- n + 1

    #separted data by biotypes

    separated[[pipe]] <- list()
    for (biotype in bio_types){
        separated[[pipe]][[biotype]] <- (data[[pipe]][[1]] %>%
        filter(str_detect(Geneid, biotype)))
    }
}


#subset the raw counts into different biotypes
#This subset method will only work for the first round of simulated data
#in the future let us keep the database source in synthetic data in order to have biotypes


#============#

#============#
#PCA



#============#

#============#
#linear regression
run_linear_regression <- FALSE

if (run_linear_regression){

organ_base <- paste("mkdir -p results/", trial, "/", sep = "")
organ_base_dir <- paste("results/", trial, "/", sep = "")

#organize files

system("mkdir -p results")


for (pipe in pipes){
    system(paste(organ_base, pipe, sep = ""))
    #will hold pdfs 
    system(paste(organ_base, pipe, "/", "overall", sep = ""))
    #will hold folders of biotypes, which hold pdfs
    system(paste(organ_base, pipe, "/", "biotypes", sep = ""))

    for (biotype in bio_types){
            system(paste(organ_base, pipe, "/", "biotypes",
            "/", biotype, sep = ""))
    }
}



for (pipe in pipes){

    #turn data loooong
    truth <- melt(data[[pipe]][[1]][2:(samples + 2)])
    experimental <- melt(data[[pipe]][[1]]
    [(samples + 8):(samples + 7 + samples)])

    tru_exp_df <- cbind(truth, experimental)

    colnames(tru_exp_df) <- c("geneid", "tru_sample",
    "tru_value", "exp_sample", "exp_value")

    #log2 transform the data
    tru_exp_df$tru_value <- log2(tru_exp_df$tru_value + 1)
    tru_exp_df$exp_value <- log2(tru_exp_df$exp_value + 1)

    overall_plot <- ggplot(tru_exp_df, aes(x = tru_value, y = exp_value)) + 
    geom_point() +
    stat_poly_line() +
    stat_poly_eq() +
    labs(x = "log2(Ground Truth Counts + 1)",
    y = paste(pipe, "log2(Experimental Counts)"),
    title = paste("Overall Experimental vs. Truth for", pipe))

    ggsave(paste(organ_base_dir, pipe, "/overall/", pipe, "_", trial,
    "_overall_linear", ".pdf", sep = ""), plot = overall_plot)
#make biotype specific plots
    for (biotype in bio_types){
        biotype_truth <- melt(separated[[pipe]][[biotype]][2:(samples + 2)])
        biotype_exp <- melt(separated[[pipe]][[biotype]]
        [(samples + 8):(samples + 7 + samples)])

#add names
        biotype_df <- cbind(biotype_truth, biotype_exp)
        colnames(biotype_df) <- c("geneid", "tru_sample",
        "tru_value", "exp_sample", "exp_value")
#log transform
        biotype_df$tru_value <- log2(biotype_df$tru_value + 1)
        biotype_df$exp_value <- log2(biotype_df$exp_value + 1)

        biotype_plot <- ggplot(biotype_df, aes(x = tru_value, y = exp_value)) + 
        geom_point() +
        stat_poly_line() +
        stat_poly_eq() +
        labs(x = "log2(Ground Truth Counts + 1)",
        y = paste(pipe, "log2(Experimental Counts)"),
        title = paste("Overall Experimental vs. Truth for", pipe, biotype))

         ggsave(paste(organ_base_dir, pipe, "/biotypes/",
         biotype, "/", pipe, "_", trial, "_", biotype,
        "_overall_linear", ".pdf", sep = ""), plot = biotype_plot)
    }
}

}
#============#
#Sensitivity and Precision Calculation
#14000 "right answers"

#sensitivity 
for (pipe in pipes){

    #turn data loooong
    truth <- melt(data[[pipe]][[1]][2:(samples + 2)])
    experimental <- melt(data[[pipe]][[1]]
    [(samples + 8):(samples + 7 + samples)])

    tru_exp_df <- cbind(truth, experimental)

    colnames(tru_exp_df) <- c("geneid", "tru_sample",
    "tru_value", "exp_sample", "exp_value")


#make biotype calculations
    for (biotype in bio_types){
        biotype_truth <- melt(separated[[pipe]][[biotype]][2:(samples + 2)])
        biotype_exp <- melt(separated[[pipe]][[biotype]]
        [(samples + 8):(samples + 7 + samples)])

#add names
        biotype_df <- cbind(biotype_truth, biotype_exp)
        colnames(biotype_df) <- c("geneid", "tru_sample",
        "tru_value", "exp_sample", "exp_value")
#log transform
        biotype_df$tru_value <- log2(biotype_df$tru_value + 1)
        biotype_df$exp_value <- log2(biotype_df$exp_value + 1)

        biotype_plot <- ggplot(biotype_df, aes(x = tru_value, y = exp_value)) + 
        geom_point() +
        stat_poly_line() +
        stat_poly_eq() +
        labs(x = "log2(Ground Truth Counts + 1)",
        y = paste(pipe, "log2(Experimental Counts)"),
        title = paste("Overall Experimental vs. Truth for", pipe, biotype))

         ggsave(paste(organ_base_dir, pipe, "/biotypes/",
         biotype, "/", pipe, "_", trial, "_", biotype,
        "_overall_linear", ".pdf", sep = ""), plot = biotype_plot)
    }
}

#============#
