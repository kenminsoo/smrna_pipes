# smRNA Pipes

This script is to run default parameters and set up pipelines related to processing small RNA sequencing data given a dir of processed fasta files. 

Pipelines/Workflows:
Pre-processing => AASRA alignment pipeline/method => featuecounts
Standard RNA-seq workflow (FeatureCounts style output compatible with comparisons.R)
COMPSRA (Non-feature count style of output, separated)

# Data Workflow
We will run many trials, increasing the complexity of datasets to see how different pipelines will respond. 

# Data Format: {trial/experiment/prefix}-{pipeline}-{data_type}.csv

## Data Types
FP = False Positives
PER_ERROR = Percentage error for each read that was supposed to be in the datasets
RAW_COUNTS = Raw counts to generate regressions to evaluate performance
RAW_ERROR = Raw difference between experimental and base truth

# Scripts for the processing of files to allow for pipeline compatability

None for now. See issues for goals. Completed scripts will be placed here with descriptions for use. 

# Current List of Compounding Complexities:

## Base-Perfect Sequencing and Aside 1 & 2

1) Base-Perfect sequencing data with no adapters or errors. No small RNA sequences that occur multiple times in the genome (According to ITAS annotation). Includes only piRNAs, miRNAs, tRNAs, and rRNAs.

Prefix: base

2) Base-Perfect sequencing data with no adapters or errors. We include small RNAs that occur multiple times in the genome (addition of tRNA fragments) to see how pipelines adapt reads that occur in multiple places in the genome. 

Prefix: 

3) Base-Perfect sequencing data with no adapters or errors. We include other small RNAs such as snoRNAs, snRNAs, and siRNAs, adding each type to see how each effects the ability of pipelines to find unique reads. 

Prefix: 

Aside Inquiry 1) Identify what transcripts were omitted from the ITAS annotation file that we use. Identify what snoRNAs, snRNAs, etc. were removed if any from our combined annotation file. Run literature review (basic search) on each to see how relevant they are and if real data is being missed out on. 

Aside Inquiry 2) Identify embedded small RNAs (i.e. piRNA derived snoRNAs) and see how pipes respond when reads from both are generated. Try to identify if there is a good way to uniquely identify or predict the presence of one or the other. 

Aside Inquiry 3) How well does the ratio-based normalization method work? 

## Addition of random error to the model 

4) 

## Addition of experiment-type specific error to model
plasma, tissue, poop(?)

## Investigation in pre-processing errors 

)

## Investigation of addition of other species to sequencing data (Via literature review of commonly found non-human species in sequencing data whether through contamination or true signal)


## Validation with Real Datasets via RT-qPCR from lab or existing datasets 

)