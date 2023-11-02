# Bidirectional-Promoters
The supplementary information to "Identifying and Understanding the Importance and Role of Bidirectional Promoters in Wheat" a Masters thesis produced by Katie Hawkins while studying at the University of Leeds in the academic year 2022/2023.

File Structure Overview

1A. bidirectionally_arranged_genes.R
R script that generates a list of bidirectionally arranged genes for wheat, rice, thale cress and other genomes of interest.
1B. Protein_Coding
Folder containing input files for the script in 1A.
1C. BiPs
Folder containing lists of bidirectional genes outputted from the Bidirectionally_Arranged_Genes.R script as csv files.
1D. Proximal
Folder containing the script Controls.R, adapted from supplementary 1A, and output file proximal_triticum_aestivum.csv which contains a list of proximally arranged genes where proximal genes are defined as protein coding genes with transcription start sites within a thousand base pairs of one another and located on the same strand.

2A. filter_if_conserved.R
R script that outputs two excel files containing orthologue information of bidirectionally arranged wheat genes
2B. all_orthologues
Folder containing wheat orthologues for bidirectionally arranged thale cress and rice genes as two text files. 
2C. BIP_Orthologue_Information.csv
Csv file that contains a list of all bidirectional wheat genes with corresponding orthologue information. The column AT_Freq indicates the number of thale cress orthologues for that gene, and AT_Confidence_Sum is a sum of the confidence indicator. The confidence indicator can only be 1 or 0. Equally OS_Freq and OS_Confidence_Sum provide the same information for rice.

3A. BIP_GO_terms.csv
A list of bidirectional genes and their associated GO terms as downloaded from EnsemblPlants using the biomart tool.
3B. BIP_GO_enrichment_biological_processes.csv
GO enrichment analysis for bidirectional genes for the domain biological processes
3C. Proximal_GO_terms.csv
A list of proximal genes and their associated GO terms as downloaded from EnsemblPlants using the biomart tool.

4A. Get_Promoter_Regions
Folder containing a python script get_promoter_regions.py that creates two fasta files, one with the promoter regions for the bidirectional pairs and one with the promoter regions for the proximal pairs. Also contains all required input text files.
4B. BIP_promoter_sequences_fowards.fa
Output fasta file containing the promoter sequences of the bidirectional pairs from the forward strand.
4C. Proximal_promoter_sequences_forwards.fa
Output fasta file containing the promoter sequences of the proximal pairs from the forward strand.

5A. coexpression.R
R script that imports all tpm data for proximal and bidirectional genes, calculates the co-expression coefficient and outputs the data as a csv file.
5B. Coexpression_Data.csv
A csv file containing the pair ID, gene names, description, and co-expression coefficient for proximal and bidirectional genes.

6A. Temperature_GO_Terms.csv
A csv file containing the GO term and description of bidirectional genes associated with one or more GO terms whose description contains one or more of the following key terms: temperature, cold, and heat.
6B. Pull_Temp_BIP_GO_Terms.R
An R script that outputs a list of bidirectionally arranged genes with temperature related GO terms as a csv file.
6C. BIP_GO_temp.csv
A list of bidirectionally arranged genes with temperature related GO terms as a csv file.

7A. Cold_differential_gene_translation.csv
Contains the mapping from the old to new genome translation for all differentially expressed genes in the cold study.
7B. Drought_heat_differential_gene_translation
Contains the mapping from the old to new genome translation for all differentially expressed genes in the drought heat study.
7C. Microspore_differential_gene_translation
Contains the mapping from the old to new genome translation for all differentially expressed genes in the microspore study.
7D. BIP_temperature_info.csv
As csv file containing, pair ID, co-expression coefficient, a note of what temperature studies the gene was differentially expressed in, and how many temperature studies, and a count value from 0-25 that indicates its likelihood to be related to temperature response for all bidirectional genes.

8A. Streme_output
A folder containing the streme output of motifs enriched in bidirectional and proximal promoters in html format.
8B. Get_Promoter_Regions
A folder containing three python scripts “get_promoter_regions_BIP.py”, “get_promoter_regions_proximal_part1.py”, and “get_promoter_regions_proximal_part2.py” that produce the promoter sequences of bidirectional and proximal genes.
8C. 1000bp_upsteam
A folder containing “BIP_1000bp_upstream_sequence_forward.txt” and “BIP_1000bp_upstream_sequence_reverse.txt”, which contain the sequences of the 1000 bp upstream of forward and reverse bidirectional promoter respectively  
8D. Bidirectional_Promoter_Sequences
A zip folder containing BIP_promoter_sequences_forwards.fa and BIP_promoter_sequences_reverse.fa
8E. proximal_promoter_sequences_head.fa
A fasta file containing the intermediate region between the tail end of one proximal gene and the TSS of the next.

9A. Gene_Sequences
A folder containing four fasta files containing bidirectional and proximal genes sequences split into forwards and reverse genes, in the case of bidirectional genes, and split into gene 1 and gene 2 in the case of proximal genes.
9B. Shell_scripts
A folder containing “align_pairs.sh” and “align_proximal_pairs.sh”, shell scripts that align the gene sequence of bidirectional and proximal genes with their pair, respectively.
9C. Shell_script_output
A folder containing two text files with the output of the align_pairs.sh and align_proximal_pairs.sh scripts, “pair_alignment_results.txt” and “proximal_pair_alignment_results.txt”. Additionally contains the full output of the bidirectional shell script, which includes the visual alignment between the two sequences “pair_alignment_full_results.txt”. 
9D. Proximal_triticum_aestivum_orientation.csv
Csv file containing a list of proximal genes along with their head/tail classification.
