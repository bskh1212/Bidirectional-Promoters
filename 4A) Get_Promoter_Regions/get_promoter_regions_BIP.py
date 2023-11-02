##create a file with promoters for all BIPs

##import csv module
import csv

##import ensembl plant data with sequence of 1000bp upstream flanking regions
##for all BIP genes 
##save data to dictionary called upstream

upstream_f = {}
gene_id = "placeholder"
gene_bases = "placeholder"

with open("BIP_1000bp_upstream_sequence_forward.txt", "r", newline='') as fasta_file:
    for line in fasta_file:
        if ">" in line:

            ##append current sequence info to the upstream_f dictionary
            upstream_f[gene_id] = gene_bases

            ##set the new gene id
            gene_id = line.strip("\n")
            gene_id = gene_id.strip(">")

            ##clear the gene bases string
            gene_bases = ""
            
        elif ">" not in line:
            ##update the gene bases string with the next line of the file
            ##as long as the file is not an id line
            gene_bases = gene_bases + line.strip("\n")
    upstream_f[gene_id] = gene_bases
        
upstream_f.pop("placeholder")

##reverse
upstream_r = {}
gene_id = "placeholder"
gene_bases = "placeholder"

with open("BIP_1000bp_upstream_sequence_reverse.txt", "r", newline='') as fasta_file:
    for line in fasta_file:
        if ">" in line:

            ##append current sequence info to the upstream_f dictionary
            upstream_r[gene_id] = gene_bases

            ##set the new gene id
            gene_id = line.strip("\n")
            gene_id = gene_id.strip(">")

            ##clear the gene bases string
            gene_bases = ""
            
        elif ">" not in line:
            ##update the gene bases string with the next line of the file
            ##as long as the file is not an id line
            gene_bases = gene_bases + line.strip("\n")
    upstream_r[gene_id] = gene_bases

upstream_r.pop("placeholder")

##import gene pairs as dictionaries
##forward containing keys=gene names and values=pair_ID on 1 strand
##reverse containing keys=gene names and values=pair_ID on -1 strand
##distance_BIP containing keys=pair_ID and values=distance between TSS for BIP
##distance_prox containing keys=pair_ID for proximal TSS
##pair containing keys=Pair_ID and values=list of two genes

reverse = {}
forward = {}
distance_BIP = {}
pair = {}

##file headers: PairID	Gene_1	1_Description	Gene_2	2_Description
##Distance	Coexpression_Coefficient	P_Value	type	Strand
##BH

with open("../../6. Co-expression/Coexpression_Data.csv") as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=",")
    for line in csv_reader:
        if line[8] == "BIP":
            forward[line[1]] = line[0] ##gene1 = pairID
            reverse[line[3]] = line[0] ##gene2 = pairID
            pair[line[0]] = [line[1], line[3]] ##pairID = gene1, gene2
            distance_BIP[line[0]] = line[5] ##PairID = distance

##at this stage the following infomation has been imported:
##1000bp upstream sequence of all BIP genes
##pair IDs and key values
##distances between TSS for all pairs
##next cut the 1000bp lengths using the distance parameter so you only
##have the sequence for the promoter

##forward
##define function to return the promoter sequence for a given bip gene
def get_promoter(gene):
    ##first obtain the pair ID
    if gene in forward:
           pair_ID = forward[gene]
    elif gene in reverse:
        pair_ID = reverse[gene]
    
    ##use the pair ID and promoter type to obtain the distance
    d = int(distance_BIP[pair_ID])

    ##get the sequence and cut it

    ##for forwards genes
    if gene in forward:
        sequence = upstream_f[gene]
        sequence = sequence[1000-d:1000]

    ##for reverse genes
    elif gene in reverse:
        sequence = upstream_r[gene]
        sequence = sequence[::-1]
        sequence = sequence[0:d]

    ##the outout of the function is the promoter sequence
    return sequence

##for all genes perform get_promoter function
##and add the output to one of three dictionaries(key=gene ID, value=promoter
##sequence) one for forwards, reverse, and proximal (forwards)

promoters_f = {}
for gene in forward:
    v = get_promoter(gene)
    pair_ID = forward[gene]
    promoters_f[pair_ID] = v
    
promoters_r = {}
for gene in reverse:
    v = get_promoter(gene)
    pair_ID = reverse[gene]
    promoters_r[pair_ID] = v


##two sequences for each promoter the forward and reverse complement
##saved to dictionaries promoters_f and promoters_r respectively

##spot check: are promoter values for forwards and reverse complimentary?
##(which they should be because the orientation of the reverse was
##reversed earlier)
print(promoters_f["Pair_201"])
print(promoters_r["Pair_201"])
print(promoters_f["Pair_654"])
print(promoters_r["Pair_654"])


##function that adds new lines every 60 characters
def insert_newlines(string, every=60):
    return '\n'.join(string[i:i+every] for i in range(0, len(string), every))

####write the forward genes to fasta file
##with open("BIP_promoter_sequences_fowards.fa", "w") as fasta_file:
##    for gene in promoters_f:
##        fasta_file.write(">"+gene+" "+pair[gene][0]+" "+pair[gene][1]+"\n")
##        fasta_file.write(insert_newlines(promoters_f[gene])+"\n")
##

