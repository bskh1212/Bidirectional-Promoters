##get the proximal promoter regions
##part 1: determine orientation relative to other gene in the pair

##head = TSS faces promoter, tail = the end of the gene faces the promoter

##note: some genes are present in two proximal pairs, in which case they are
##both in a head AND tail orientation depending on which pair you're condidering

##import modules
import csv

##create dictionaries to store proximal gene info
pair_to_genes = {}
gene_to_pair = {}

strand = {}
gene_start = {}
gene_end = {}
TSS = {}
distance = {}

##read info from file

##headers
##[0]Gene_ID
##[1]Strand
##[2]Chromosome
##[3]Gene_Description
##[4]Gene_Start
##[5]Gene_End
##[6]TSS
##[7]Distance
##[8]PairID

rows = []
with open("../../6. Co-expression/proximal_triticum_aestivum.csv") as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=",")
    for line in csv_reader:
        gene = line[0]
        strand[gene] = line[1]
        gene_start[gene] = line[4]
        gene_end[gene] = line[5]
        TSS[gene] = line[6]
        distance[gene] = [line[7]]

        pair_ID = line[8]
        ##gene to pair
        gene_to_pair[gene] = pair_ID
        ##pair to genes
        ##does the key already exist?
        existing = pair_to_genes.get(pair_ID)

        if existing:
            pair_to_genes[pair_ID].append(gene)
        else: pair_to_genes[pair_ID] = [gene]

        rows.append([line[0], line[1], line[2], line[3], line[4],
                    line[5], line[6], line[7], line[8]])
   
##determine orientation
##key must be gene AND pair ID
orientation = {}

for pair, gene in pair_to_genes.items():
    if pair == "PairID":
        pass
    else:
        first = int(TSS[gene[0]]) * int(strand[gene[0]])
        second = int(TSS[gene[1]]) * int(strand[gene[1]])

        ##define keys
        key0 = gene[0] + "_" + pair
        key1 = gene[1] + "_" + pair

            
        if second > first:
            orientation[key0] = "tail"
            orientation[key1] = "head"
        elif second < first:
            orientation[key0] = "head"
            orientation[key1] = "tail"

##write file with added orientation information
with open("../../6. Co-expression/proximal_triticum_aestivum_orientation.csv", 'w', newline="") as csv_file:
    writer = csv.writer(csv_file)
    writer.writerow(["Gene_ID",
                     "Strand",
                     "Chromosome",
                     "Gene_Description",
                     "Gene_Start",
                     "Gene_End",
                     "TSS",
                     "Distance",
                     "PairID",
                     "Orientation"])
    for n in rows:
        if n[0] != "Gene_ID":

            key = n[0] + "_" + n[8]
            orientation_value = orientation[key]
            n.append(orientation_value)

            writer.writerow(n)

csv_file.close()
