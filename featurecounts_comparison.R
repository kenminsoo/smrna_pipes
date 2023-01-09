#This script will compare featurecounts outputs
#with the ground truth.

#expected input will be a csv file:
#Geneid, sample1, ..., samplen, chr, start, end, strand, length, exp_sample1, ..., exp_samplen # nolint

#For now we will produce the following outputs:

#PCA of the samples from all pipes.
#Linear regression of ground truth with experimental.
#Sensitivity = TP / TP + FN # nolint

#if read found, it will be counted as TP
#if read not found, it will be counted as FN

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
library(ggplot2)
library(viridis)
library(hrbrthemes)
library(ggsci)

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
fp_separated <- list()
per_separated <- list()

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
    fp_separated[[pipe]] <- list()
    per_separated[[pipe]] <- list()
    for (biotype in bio_types){
        separated[[pipe]][[biotype]] <- (data[[pipe]][[1]] %>%
        filter(str_detect(Geneid, biotype)))

        fp_separated[[pipe]][[biotype]] <- (data[[pipe]][[3]] %>%
        filter(str_detect(Geneid, biotype)))

        per_separated[[pipe]][[biotype]] <- (data[[pipe]][[2]] %>%
        filter(str_detect(X1, biotype)))
    }
}

#subset the raw counts into different biotypes
#This subset method will only work for the first round of simulated data
#in the future let us keep the database source
#in synthetic data in order to have biotypes

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

#create dir if does not exist
organ_base <- paste("mkdir -p results/", trial, "/", sep = "")
organ_base_dir <- paste("results/", trial, "/", sep = "")

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

#make empty csv files for sensitivity and precision + percent error calculations
#this isn't actually necessary
calc_base <- paste("touch results/", trial, "/", sep = "")
calc_base_dir <- paste("results/", trial, "/", sep = "")

for (pipe in pipes){
    #will hold overall metrics
    system(paste(calc_base, pipe, "/", "overall", "/", pipe,
    "_calculated.csv", sep = ""))

    for (biotype in bio_types){
            system(paste(calc_base, pipe, "/", "biotypes", "/", biotype,
            "/", pipe, "_", biotype, "_", "_calculated.csv", sep = ""))
    }
}

#sensitivity and precision and average percent error

system(paste("echo ", "measurement, value, pipeline, biotype", ">", " ",
organ_base_dir, "biotype_calculated.csv", sep = ""))

for (pipe in pipes){

    #turn data loooong
    truth <- melt(data[[pipe]][[1]][2:(samples + 2)])
    experimental <- melt(data[[pipe]][[1]]
    [(samples + 8):(samples + 7 + samples)])

    tru_exp_df <- cbind(truth, experimental)

    colnames(tru_exp_df) <- c("geneid", "tru_sample",
    "tru_value", "exp_sample", "exp_value")

    false_neg <- colSums(tru_exp_df["exp_value"] == 0)
    true_pos <- colSums(tru_exp_df["exp_value"] != 0)

    sensitivity <- ((true_pos) / (true_pos + false_neg))

    #precision
    fp_df <- melt(data[[pipe]][[3]][(samples + 8):(samples + 7 + samples)])
    false_pos <- colSums(fp_df["value"] != 0) 

    precision <- ((true_pos) / (false_pos + true_pos))

    #percent error
    per_df <- melt(data[[pipe]][[2]][2:(samples + 2)])
    per_df_non_zero <- per_df[3][per_df[3] != 0]
    per_avg0_pos <- mean(per_df_non_zero[per_df_non_zero > 0])
    per_sd0_pos <- sd(per_df_non_zero[per_df_non_zero > 0])
    per_avg0_neg <- mean(per_df_non_zero[per_df_non_zero < 0 & per_df_non_zero != 0])
    per_sd0_neg <- sd(per_df_non_zero[per_df_non_zero < 0 & per_df_non_zero != 0])
    per_avg <- mean(per_df[[3]])
    per_sd <- sd(per_df[[3]])
    f_measure <- (2 * (precision * sensitivity)) / (precision + sensitivity)
    #enter data into csv

    system(paste("echo ", "sensitivity", ",", sensitivity, " >", " ",
    organ_base_dir, pipe, "/overall/", pipe, "_calculated.csv", sep = ""))

    system(paste("echo ", "precision", ",", precision, " >>", " ",
    organ_base_dir, pipe, "/overall/", pipe, "_calculated.csv", sep = ""))

    system(paste("echo ", "overall_error", ",", per_avg, " >>", " ",
    organ_base_dir, pipe, "/overall/", pipe, "_calculated.csv", sep = ""))

    system(paste("echo ", "overall_sd", ",", per_sd, " >>", " ",
    organ_base_dir, pipe, "/overall/", pipe, "_calculated.csv", sep = ""))

    system(paste("echo ", "pos_error", ",", per_avg0_pos, " >>", " ",
    organ_base_dir, pipe, "/overall/", pipe, "_calculated.csv", sep = ""))

    system(paste("echo ", "pos_sd", ",", per_sd0_pos, " >>", " ",
    organ_base_dir, pipe, "/overall/", pipe, "_calculated.csv", sep = ""))

    system(paste("echo ", "neg_error", ",", per_avg0_neg, " >>", " ",
    organ_base_dir, pipe, "/overall/", pipe, "_calculated.csv", sep = ""))

    system(paste("echo ", "neg_sd", ",", per_sd0_neg, " >>", " ",
    organ_base_dir, pipe, "/overall/", pipe, "_calculated.csv", sep = ""))

    system(paste("echo ", "f_measure", ",", f_measure, " >>", " ",
    organ_base_dir, pipe, "/overall/", pipe, "_calculated.csv", sep = ""))

    #enter into the big csv

    system(paste("echo ", "sensitivity", ",", sensitivity, 
    ",", pipe, ",overall", " >>", " ",
    organ_base_dir, "biotype_calculated.csv", sep = ""))

    system(paste("echo ", "precision", ",", precision, 
    ",", pipe, ",overall", " >>", " ",
    organ_base_dir, "biotype_calculated.csv", sep = ""))

    system(paste("echo ", "overall_error", ",", per_avg,
    ",", pipe, ",overall", " >>", " ",
    organ_base_dir, "biotype_calculated.csv", sep = ""))

    system(paste("echo ", "overall_sd", ",", per_sd,
    ",", pipe, ",overall", " >>", " ",
    organ_base_dir, "biotype_calculated.csv", sep = ""))

    system(paste("echo ", "pos_error", ",", per_avg0_pos, 
    ",", pipe, ",overall", " >>", " ",
    organ_base_dir, "biotype_calculated.csv", sep = ""))

    system(paste("echo ", "pos_sd", ",", per_sd0_pos,
    ",", pipe, ",overall", " >>", " ",
    organ_base_dir, "biotype_calculated.csv", sep = ""))

    system(paste("echo ", "neg_error", ",", per_avg0_neg,
    ",", pipe, ",overall", " >>", " ",
    organ_base_dir, "biotype_calculated.csv", sep = ""))

    system(paste("echo ", "neg_sd", ",", per_sd0_neg,
    ",", pipe, ",overall", " >>", " ",
    organ_base_dir, "biotype_calculated.csv", sep = ""))

    system(paste("echo ", "f_measure", ",", f_measure,
    ",", pipe, ",overall", " >>", " ",
    organ_base_dir, "biotype_calculated.csv", sep = ""))

#make biotype calculations
    for (biotype in bio_types){
        biotype_truth_df <- melt(separated[[pipe]][[biotype]][2:(samples + 2)])
        biotype_exp_df <- melt(separated[[pipe]][[biotype]]
        [(samples + 8):(samples + 7 + samples)])

#add names
        biotype_df <- cbind(biotype_truth_df, biotype_exp_df)
        colnames(biotype_df) <- c("geneid", "tru_sample",
        "tru_value", "exp_sample", "exp_value")

        false_neg <- colSums(biotype_df["exp_value"] == 0)
        true_pos <- colSums(biotype_df["exp_value"] != 0)

        sensitivity <- ((true_pos) / (true_pos + false_neg))

        print(biotype)
        print(true_pos/10)
        #precision

        fp_df <- melt(fp_separated[[pipe]][[biotype]]
        [(samples + 8):(samples + 7 + samples)])
        false_pos <- colSums(fp_df["value"] != 0)

        precision <- ((true_pos) / (false_pos + true_pos))

        #percent error
        per_df <- melt(per_separated[[pipe]][[biotype]][2:(samples + 2)])
        per_df_non_zero <- per_df[3][per_df[3] != 0]
        per_avg0_pos <- mean(per_df_non_zero[per_df_non_zero > 0])
        per_sd0_pos <- sd(per_df_non_zero[per_df_non_zero > 0])
        per_avg0_neg <- mean(per_df_non_zero[per_df_non_zero < 0 & per_df_non_zero != -1])
        per_sd0_neg <- sd(per_df_non_zero[per_df_non_zero < 0 & per_df_non_zero != -1])
        per_avg <- mean(per_df[[3]])
        per_sd <- sd(per_df[[3]])
        f_measure <- (2 * (precision * sensitivity)) / (precision + sensitivity)

        system(paste("echo ", "sensitivity", ",", sensitivity, " >", " ",
        organ_base_dir, pipe, "/", "biotypes", "/", biotype,
        "/", pipe, "_", biotype, "_", "_calculated.csv", sep = ""))

        system(paste("echo ", "precision", ",", precision, " >>", " ",
        organ_base_dir, pipe, "/", "biotypes", "/", biotype,
        "/", pipe, "_", biotype, "_", "_calculated.csv", sep = ""))

        system(paste("echo ", "overall_error", ",", per_avg, " >>", " ",
        organ_base_dir, pipe, "/", "biotypes", "/", biotype,
        "/", pipe, "_", biotype, "_", "_calculated.csv", sep = ""))

        system(paste("echo ", "overall_sd", ",", per_sd, " >>", " ",
        organ_base_dir, pipe, "/", "biotypes", "/", biotype,
        "/", pipe, "_", biotype, "_", "_calculated.csv", sep = ""))

        system(paste("echo ", "pos_error", ",", per_avg0_pos, " >>", " ",
        organ_base_dir, pipe, "/", "biotypes", "/", biotype,
        "/", pipe, "_", biotype, "_", "_calculated.csv", sep = ""))

        system(paste("echo ", "pos_sd", ",", per_sd0_pos, " >>", " ",
        organ_base_dir, pipe, "/", "biotypes", "/", biotype,
        "/", pipe, "_", biotype, "_", "_calculated.csv", sep = ""))

        system(paste("echo ", "neg_error", ",", per_avg0_neg, " >>", " ",
        organ_base_dir, pipe, "/", "biotypes", "/", biotype,
        "/", pipe, "_", biotype, "_", "_calculated.csv", sep = ""))

        system(paste("echo ", "neg_sd", ",", per_sd0_neg, " >>", " ",
        organ_base_dir, pipe, "/", "biotypes", "/", biotype,
        "/", pipe, "_", biotype, "_", "_calculated.csv", sep = ""))

        system(paste("echo ", "f_measure", ",", f_measure, " >>", " ",
        organ_base_dir, pipe, "/", "biotypes", "/", biotype,
        "/", pipe, "_", biotype, "_", "_calculated.csv", sep = ""))

#write to a csv file for all pipes for easy analysis

        system(paste("echo ", "sensitivity", ",", 
        sensitivity, ",", pipe, ",", biotype, " >>", " ",
        organ_base_dir, "biotype_calculated.csv", sep = ""))

        system(paste("echo ", "precision", ",", 
        precision, ",", pipe, ",", biotype, " >>", " ",
        organ_base_dir, "biotype_calculated.csv", sep = ""))

        system(paste("echo ", "overall_error", ",", per_avg, 
        ",", pipe, ",", biotype, " >>", " ",
        organ_base_dir, "biotype_calculated.csv", sep = ""))

        system(paste("echo ", "overall_sd", ",", 
        per_sd, ",", pipe, ",", biotype, " >>", " ",
        organ_base_dir, "biotype_calculated.csv", sep = ""))

        system(paste("echo ", "pos_error", ",", 
        per_avg0_pos, ",", pipe, ",", biotype, " >>", " ",
        organ_base_dir, "biotype_calculated.csv", sep = ""))

        system(paste("echo ", "pos_sd", ",", 
        per_sd0_pos, ",", pipe, ",", biotype, " >>", " ",
        organ_base_dir, "biotype_calculated.csv", sep = ""))

        system(paste("echo ", "neg_error", ",", 
        per_avg0_neg, ",", pipe, ",", biotype, " >>", " ",
        organ_base_dir, "biotype_calculated.csv", sep = ""))

        system(paste("echo ", "neg_sd", ",", 
        per_sd0_neg, ",", pipe, ",", biotype, " >>", " ",
        organ_base_dir, "biotype_calculated.csv", sep = ""))

        system(paste("echo ", "f_measure", ",", 
        f_measure, ",", pipe, ",", biotype, " >>", " ",
        organ_base_dir, "biotype_calculated.csv", sep = ""))
    }
}

#make graphs that show biotype and pipeline specific metrics
#import biotype calculated csv

metrics_df <- read.csv(paste(organ_base_dir, 
"biotype_calculated.csv", sep = ""))

measurement_list <- unique(metrics_df[1])
measurement_list <- list(measurement_list)
measurement_list <- unlist(measurement_list)

#we will have each be made into bar graph, with %error, and sd
n <- 1
for (metric in measurement_list) {

    if (grepl("error", metric)) {
        sd <- measurement_list[n + 1]

        filtered_dataset <- metrics_df %>% filter(str_detect(measurement, metric))
        sd_val <- metrics_df %>% filter(str_detect(measurement, sd))

        colnames(sd_val) <- c("measurement", "sd", "pipeline", "biotype")

        combined <- cbind(filtered_dataset, sd_val["sd"])

        metrics_plot <- ggplot(combined, aes(fill = biotype, 
        x = measurement, y = value)) +
        geom_bar(position = "dodge", stat = "identity") +
        ggtitle(paste("Pipeline ", metric, " for Different Biotypes", sep = "")) +
        facet_wrap(~pipeline) +
        ylab(metric) +
        scale_fill_jama()

        ggsave(paste(organ_base_dir, "overall_",
        metric, ".pdf", sep = ""), plot = metrics_plot)

    

    } else if (grepl("sd", metric)) {
        
        next

    } else {
    filtered_dataset <- metrics_df %>% filter(str_detect(measurement, metric))

    metrics_plot <- ggplot(filtered_dataset, aes(fill = biotype, 
    x = measurement, y = value)) +
    geom_bar(position = "dodge", stat = "identity") +
    ggtitle(paste("Pipeline ", metric, " for Different Biotypes", sep = "")) +
    facet_wrap(~pipeline) +
    ylab(metric) +
    scale_fill_jama()

    ggsave(paste(organ_base_dir, "overall_", 
    metric, ".pdf", sep = ""), plot = metrics_plot)
    }

    n <- n + 1
}

#============#
