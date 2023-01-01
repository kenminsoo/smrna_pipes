#!/bin/bash

#this script will run a series of pipes 
#use with pre processed data if not simulated

###=======================###
#CONFIG START

#Data location
dir="/home/kenmn/nfs_scratch/simulated_data_1"
#Trial/Experiment/Prefix
pre="base"
#results dir name
res_dir="/home/kenmn/nfs_scratch/sim_results"

#genome path (fasta) for indexing when neeeded
genome_path=/home/kenmn/bioinfo_tools/hg38/hg38.fasta
genome_name=hg38
annotation_dir=/home/kenmn/bioinfo_tools/smrnaseq_annotation/

#threads
threads=2

#the data folder ($dir) will have references for the following tools:
#bowtie

module load lang/Anaconda3/2022.05
eval "$(conda shell.bash hook)"

#CONFIG END
###=======================###
#AASRA START

export PATH=$PATH:/home/kenmn/nfs_scratch/temp_pipes/AASRA
#meant to be used with MANA HPC, comment this code out and activate following tools manually
#2 cpu
#64G

#activate conda environment with subread & seqtk
#conda environments will be exported and put into the repo when complete
conda activate aasra
module load bio/Bowtie2/2.3.4.2-intel-2018.5.274

cd $dir

#convert to fastq
for file in *.fasta
do
    base=${file%.*}
    seqtk seq -F "I" $file > $base".fastq"
done

#with the given reference 
#TRY TO FIX-featureCounts runs into an issue with this particular index
#I think we need to shorten the saf file. 

###AASRA-index -i human_smallRNA_reference.fa -l CCCCCCCCCC -r GGGGGGGGGG -s index.saf

###for file in *.fastq
###do
###    AASRA -p 1 -i $file -l CCC -r GGG -b "/home/kenmn/nfs_scratch/temp_pipes/AASRA/anchored_human_smallRNA_reference.fa"
###done

#with a generated reference
#This generated reference is from the ITAS annotation paper, including piRNA, miRNA, rRNA, and tRNA
#We will introduce snoRNAs, snRNAs, (from ensembl) when we begin to increase complexity of data
AASRA-index -i unique_transcripts.fa -l CCCCCCCCCC -r GGGGGGGGGG -s index_uni.saf

for file in *.fastq
do
    AASRA -p 1 -i $file -l CCC -r GGG -b "/home/kenmn/nfs_scratch/temp_pipes/AASRA/anchored_unique_transcripts.fa"
done

#now count the transcripts and compare to the original
###featureCounts -T 2 -a "/home/kenmn/nfs_scratch/temp_pipes/AASRA/index_base.saf" -F SAF -o counts_base.tsv *.sam

featureCounts -T 2 -a "/home/kenmn/nfs_scratch/temp_pipes/AASRA/index_uni.saf" -F SAF -o counts_uni.tsv *.sam


#organize data 
#pwd -> is $dir

rm anchored*
mv counts_uni* ..
cd ..
mkdir $res_dir
mv counts_uni* $res_dir
mkdir $res_dir"/AASRA"
mkdir $res_dir"/AASRA/"$pre
mv $res_dir/counts_uni* $res_dir"/AASRA/"$pre

cd $dir

#AASRA END
###=======================###

module purge
conda deactivate

###=======================###
#MANUAL PIPELINE START

conda activate smRNAseq

bowtie-build --threads $threads $genome_path $genome_name

for f in *.fastq; do
    file="$f"
    sampname="${file%%.*}"
    bowtie -x $genome_name $f -v 1 -S "${sampname}".sam -p $threads --reorder
done

featureCounts -T $threads -a $annotation_dir"human_allRNA.gtf" -F 'GTF' -g 'transcript_id' -o 'base-manual-counts.tsv' *.sam -O -M --fraction

#organize data

mv base-manual-counts.* $res_dir
mkdir $res_dir"/manual"
mkdir $res_dir"/manual/"$pre
mv $res_dir/base-manual-counts.* $res_dir"/manual/"$pre
rm *.sam

#MANUAL PIPELINE END
###=======================###

module purge
conda deactivate

###=======================###
#COMPSRA START

compsra_dir=/home/kenmn/nfs_scratch/temp_pipes/COMPSRA

module load lang/Java/16.0.1

cp $compsra_dir ../temp_pipes -r

cd ../temp_pipes/COMPSRA

java -jar COMPSRA.jar -tk -dr -ck miRNA_hg38,piRNA_hg38,tRNA_hg38,snoRNA_hg38,snRNA_hg38,circRNA_hg38

java -jar COMPSRA.jar -tk -dr -ck star_hg38

for f in $dir/*.fastq; do
    echo $f >> sample.list
done

mkdir exp_out

#after running for the first time, do not run again with the alignment mbi option

java -jar COMPSRA.jar -ref hg38 -qc -ra TGGAATTCTCGGGTGCCAAGG -rb 4 -rh 20 -rt 20 -rr 20 -rlh 8,17 -aln -mt star -mbi -ann -ac 1,2,3,4,5,6 -inf sample.list -out ./exp_out/

#uncomment if first time 
##java -jar COMPSRA.jar -ref hg38 -qc -rb 4 -rh 20 -rt 20 -rr 20 -rlh 8,17 -aln -mt star -mbi -ann -ac 1,2,3,4,5,6 -inf sample.list -out ./exp_out/

java -jar COMPSRA.jar -ref hg38 -qc -rb 4 -rh 20 -rt 20 -rr 20 -rlh 8,17 -aln -mt star -ann -ac 1,2,3,4,5,6 -inf sample.list -out ./exp_out/

java -jar COMPSRA.jar -ref hg38 -fun -fd -fdclass 1,2,3,4,5,6 -fdcase 1-5 -fdctrl 6-10 -fdnorm cpm -fdtest mwu -fdann -pro COMPSRA_DEG -inf sample.list -out ./exp_out/

java -jar COMPSRA.jar -ref hg38 -fun -fm -fms 1-10 -fdclass 1,2,3,4,5,6 -fdann -pro COMPSRA_MERGE -inf sample.list -out ./exp_out/

#We obtain some count txt files from the COMPSRA pipeline...
#However, to determine whether the false positives detected are truly FP 
#We will use featurecounts on the mapped bam files to see if our annotation file
#Is not covering everything, thereby missing annotations

conda activate smRNAseq

cd $compsra_dir
mkdir featurecount
for f in exp_out/sample*; do
    for bam in $f"/sample*_Aligned.out.bam"; do
        cp $bam featurecount
    done
done

cd featurecount

featureCounts -T $threads -a $annotation_dir"human_allRNA.gtf" -F 'GTF' -g 'transcript_id' -o 'base-compsra-counts.tsv' *.bam -O -M --fraction

mv base-compsra-counts.* $res_dir
mkdir $res_dir"/compsra"
mkdir $res_dir"/compsra/"$pre
mv $res_dir/base-compsra-counts.* $res_dir"/compsra/"$pre

#COMPSRA END
###=======================###
