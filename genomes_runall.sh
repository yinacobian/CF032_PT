#!/bin/bash

#bash genomes_runall.sh IDS.txt /home/acobian/CF032_2020 40

#To run do: thisscrit.sh [IDS.txt] [path to main folder] [number of threads to use in the system]

#1) Create a working directory: mkdir CF032_2020
#2) Save this file in the folder CF032_2020
#3) Create a folder for the raw reads, inside your main project folder: mkdir CF032_2020/P00_raw
#4) Put all reads files, including R1 and R1 inside the folder CF032_2020/P00_raw
#5) Create a list with the sample names and save it as IDS.txt, place it in the folder CF032_2020
## For example use this command line: ls P00_raw/ | cut -d '_' -f 1,2 | uniq > IDS.txt
#6) Open a screes session: screen -DR CF032 
#7) Run the command: bash genomes_runall.sh IDS.txt /home/acobian/CF032_2020 40

#1.- Quality filtering pair end : prinseq++
mkdir $2/P01_prinseq_output
cat $1 | xargs -t -I{fileID} sh -c "prinseq++ -fastq $2/P00_raw/{fileID}_R1_001.fastq -fastq2 $2/P00_raw/{fileID}_R2_001.fastq -lc_entropy=0.5 -trim_qual_right=15 -trim_qual_left=15 -trim_qual_type mean -trim_qual_rule lt -trim_qual_window 2 -min_len 30 -min_qual_mean 20  -rm_header -out_name $2/P01_prinseq_output/{fileID} -threads $3 -out_format 1"
#prinseq++ -fastq /home/acobian/CF032_2020/P00_raw/AC5204_S63_R1_001.fastq -fastq2 /home/acobian/CF032_2020/P00_raw/AC5204_S63_R2_001.fastq -lc_entropy=0.5 -trim_qual_right=15 -trim_qual_left=15 -trim_qual_type mean -trim_qual_rule lt -trim_qual_window 2 -min_len 30 -min_qual_mean 20  -rm_header -out_name /home/acobian/CF032_2020/P01_prinseq_output/AC5204_S63 -threads 40 -out_format 1

#2.- Denovo assembly and comparison to NT
mkdir $2/P02_denovo
cat $1 | xargs -I{fileID} sh -c "spades.py -1 $2/P01_prinseq_output/{fileID}_good_out_R1.fasta -2 $2/P01_prinseq_output/{fileID}_good_out_R2.fasta -t $3 --only-assembler -o $2/P02_denovo/spades_{fileID}"
#make sure you have the script removesmalls.pl in your home
cat $1 | xargs -I{fileID} sh -c "perl /home/acobian/bin/removesmalls.pl 900 $2/P02_denovo/spades_{fileID}/contigs.fasta > $2/P02_denovo/more900_contigs_{fileID}.fasta"

cat $1 | xargs -I{fileID} sh -c "blastn -query $2/P02_denovo/more900_contigs_{fileID}.fasta -db /home/DATABASES/blast/nt/nt -out $2/P02_denovo/vs_NT_more900_contigs_{fileID}.blastn -evalue 0.1 -num_threads $3 -max_target_seqs 1 -outfmt '6 qseqid sseqid pident length mismatchgapopen qstart qend sstart send evalue bitscore sskingdoms sscinames'"
cat $1 | xargs -I{fileID} sh -c "perl /home/acobian/bin/MYSCRIPTS/besthitblast.pl $2/P02_denovo/vs_NT_more900_contigs_{fileID}.blastn > $2/P02_denovo/besthit_vs_NT_more900_contigs_{fileID}.blastn"

#cat $1 | xargs -I{fileID} sh -c "cut -f 1 $2/P07_denovo/besthit_vs_NT_more900_contigs_{fileID}.blastn > $2/P07_denovo/toremove_{fileID}.txt"
#cat $1 | xargs -I{fileID} sh -c "perl /home/acobian/bin/remove_fasta_id.pl $2/P07_denovo/toremove_{fileID}.txt $2/P07_denovo/more900_contigs_{fileID}.fasta > $2/P07_denovo/nohits_more900_contigs_{fileID}.fasta"
#cat $1 | xargs -I{fileID} sh -c "perl /home/acobian/bin/removesmalls.pl 2000 $2/P07_denovo/nohits_more900_contigs_{fileID}.fasta > $2/P07_denovo/nohits_more2000_contigs_{fileID}.fasta"

### rename the contigs to include the sample ID. 
### submitt to RAST and PATRIC for annotation 
### get the CDS 


#3.- Denovo assembly with a small amount of reads, n=100,000
mkdir $2/P03_subsample_denovo
#subsample 50,000 reads
cat $1 | xargs -I{fileID} sh -c "perl /home/acobian/bin/MYSCRIPTS/random-sample-fasta.pl -i $2/P01_prinseq_output/{fileID}_good_out_R1.fasta -o $2/P03_subsample_denovo/{fileID}_good_out_R1_subsample_50000.fasta -r -n 50000"
#subsample 100,000 reads
cat $1 | xargs -I{fileID} sh -c "perl /home/acobian/bin/MYSCRIPTS/random-sample-fasta.pl -i $2/P01_prinseq_output/{fileID}_good_out_R1.fasta -o $2/P03_subsample_denovo/{fileID}_good_out_R1_subsample_100000.fasta -r -n 100000"

#denovo assemble
cat $1 | xargs -I{fileID} sh -c "spades.py -s $2/P03_subsample_denovo/{fileID}_good_out_R1_subsample_50000.fasta --only-assembler  -t $3 -o $2/P03_subsample_denovo/spades_50000_{fileID}"
cat $1 | xargs -I{fileID} sh -c "spades.py -s $2/P03_subsample_denovo/{fileID}_good_out_R1_subsample_100000.fasta --only-assembler  -t $3 -o $2/P03_subsample_denovo/spades_100000_{fileID}"

#remove smalls
cat $1 | xargs -I{fileID} sh -c "perl /home/acobian/bin/removesmalls.pl 900 $2/P03_subsample_denovo/spades_50000_{fileID}contigs.fasta > $2/P03_subsample_denovo/spades50000_more900_contigs_{fileID}.fasta"
cat $1 | xargs -I{fileID} sh -c "perl /home/acobian/bin/removesmalls.pl 900 $2/P03_subsample_denovo/spades_100000_{fileID}contigs.fasta > $2/P03_subsample_denovo/spades100000_more900_contigs_{fileID}.fasta"

#blastN vs NT 
cat $1 | xargs -I{fileID} sh -c "blastn -query $2/P03_subsample_denovo/spades50000_more900_contigs_{fileID}.fasta -db /home/DATABASES/blast/nt/nt -out $2/P03_subsample_denovo/vs_NT_spades50000_more900_contigs_{fileID}.blastn -evalue 0.1 -num_threads $3 -max_target_seqs 1 -outfmt '6 qseqid sseqid pident length mismatchgapopen qstart qend sstart send evalue bitscore sskingdoms sscinames'"
cat $1 | xargs -I{fileID} sh -c "perl /home/acobian/bin/MYSCRIPTS/besthitblast.pl $2/P03_subsample_denovo/vs_NT_spades50000_more900_contigs_{fileID}.blastn > $2/P03_subsample_denovo/besthit_vs_NT_spades50000_more900_contigs_{fileID}.blastn"
cat $1 | xargs -I{fileID} sh -c "blastn -query $2/P03_subsample_denovo/spades100000_more900_contigs_{fileID}.fasta -db /home/DATABASES/blast/nt/nt -out $2/P03_subsample_denovo/vs_NT_spades100000_more900_contigs_{fileID}.blastn -evalue 0.1 -num_threads $3 -max_target_seqs 1 -outfmt '6 qseqid sseqid pident length mismatchgapopen qstart qend sstart send evalue bitscore sskingdoms sscinames'"
cat $1 | xargs -I{fileID} sh -c "perl /home/acobian/bin/MYSCRIPTS/besthitblast.pl $2/P03_subsample_denovo/vs_NT_spades100000_more900_contigs_{fileID}.blastn > $2/P03_subsample_denovo/besthit_vs_NT_spades100000_more900_contigs_{fileID}.blastn"

