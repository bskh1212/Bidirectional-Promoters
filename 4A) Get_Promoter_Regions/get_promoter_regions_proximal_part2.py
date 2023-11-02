##get proximal promoter regions
import fasta
import pandas as pd
import math

##part 2: cut the correct sequence out of the 1000bp upstream flanking region

##import fasta file with the sequences
sequences = fasta.read_fasta("1000bp_upstream_proximal_head_genes.txt")

##import the gene info using pandas 
df = pd.read_csv("proximal_triticum_aestivum_orientation.csv")

##subset to get heads/tails rows in different data sets
tails = df[df["Orientation"] == "tail"]
heads = df[df["Orientation"] == "head"]

##work out gene length for tail genes
gene_length_dict = {}

def populate_dict(row):
    key = row["Gene_ID"]
    gene_length = abs(row["Gene_End"] - row["Gene_Start"])
    gene_length_dict[key] = gene_length

df.apply(lambda row: populate_dict(row), axis=1)

##pair ID distance dict
distance_dict = {}

def populate_dict(row):
    if math.isnan(row["Distance"]):
        pass
    else:
        row["Distance"] != "NaN"
        key = row["PairID"]
        value = row["Distance"]

        distance_dict[key] = value

df.apply(lambda row: populate_dict(row), axis=1)

##calculate promoter length for all tail genes
promoter_length_dict = {}

def get_promoter_length(row):
    gene_ID = row["Gene_ID"]
    pairID = row["PairID"]
    ##get distance and length from respective dicts
    distance = distance_dict[pairID]
    gene_length = gene_length_dict[gene_ID]
    ##calculate promoter length
    promoter_length = distance - gene_length

    ##populate dict
    promoter_length_dict[pairID] = promoter_length

tails.apply(lambda row: get_promoter_length(row), axis=1)

##create a head gene:pair ID dict
head_to_pair = {}

def populate_dict(row):
    key = row["Gene_ID"]
    value = row["PairID"]

    head_to_pair[key] = value

heads.apply(lambda row: populate_dict(row), axis=1)

##create a new sequence dict pairID:sequence
sequences1 = {}

for k, v in sequences.items():
    pairID = head_to_pair[k]
    sequences1[pairID] = v

##for every pairID cut sequence (1000-promoter length):1000
promoter_sequence_dict = {}
    
for pairID, sequence in sequences1.items():
    promoter_length = int(promoter_length_dict[pairID])
    promoter_seq = sequence[1000-promoter_length:1000]
    promoter_sequence_dict[pairID] = promoter_seq


##write the promoter sequences to a fasta file
def insert_newlines(string, every=60):
    return '\n'.join(string[i:i+every] for i in range(0, len(string), every))

with open("proximal_promoter_sequences_head.fa", "w") as fasta_file:
    for pairID, seq in promoter_sequence_dict.items():        
        if len(seq) >= 100:
            fasta_file.write(">"+pairID+"\n")
            fasta_file.write(insert_newlines(seq+"\n"))
        else: pass

##count number of promoter regions >= 100bp
count = 0
for seq in promoter_sequence_dict.values():
    if len(seq) >= 100:
        count = count + 1
    else: pass
print(count)


