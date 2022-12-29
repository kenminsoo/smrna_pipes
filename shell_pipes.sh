#!/bin/bash

#this script will run a series of pipes 

dir="/home/kenmn/nfs_scratch/simulated_data_1/*.fasta"

export PATH=$PATH:/home/kenmn/nfs_scratch/temp_pipes/AASRA
#meant to be used with MANA HPC, comment this code out and activate following tools manually

module load lang/Anaconda3/2022.05
eval "$(conda shell.bash hook)"
conda activate aasra
module load bio/Bowtie2/2.3.4.2-intel-2018.5.274

#with the given reference

AASRA-index -i human_smallRNA_reference.fa -l CCCCCCCCCC -r GGGGGGGGGG -s index.saf

for file in *.fasta
do
    base=${file%.*}
    seqtk seq -F "#" $file > $base".fastq"
done

for file in *.fastq
do
    AASRA -p 1 -i $file -l CCC -r GGG -b "/home/kenmn/nfs_scratch/temp_pipes/AASRA/anchored_human_smallRNA_reference.fa"
done

#with a generated reference