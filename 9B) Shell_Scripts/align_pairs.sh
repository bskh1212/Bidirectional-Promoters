##Align Pairs

##Run with current environment (-V) and in the current directory (-cwd)
#$ -V -cwd
##Request time
#$ -l h_rt=00:15:00
##Request some memory per core
#$ -l h_vmem=1G
##Get email at start and end of job
#$ -m be

##import anaconda module
module add anaconda

##initialise
conda init bash

##create environment
##conda create --name blast_env

##load environment
source activate blast_env

##use anaconda to import blast
conda install -c bioconda blast

##align pairs
blastn -query bip_forward_sequences.fa -subject bip_reverse_sequences.fa -out pair_alignment_results.txt -outfmt "6 qseqid sseqid pident qcovs evalue" -max_hsps 1
