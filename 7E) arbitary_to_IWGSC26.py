##Translate the arbitrary codes from the differential expression analysis
##from the supplementary data of the microspore paper from the wheat
##expression browser

##import required modules for handling of csv files and regex expressions
import csv
import re

##create a dictionary where the key is the arbitrary code used in the paper
##and the value is its gene name in the IWGSC26 genome assembly
##this information is pulled from the document
##"transcript_gene_name_conversion.csv"
translate = {}
with open("transcript_gene_name_conversion.csv") as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=",")
    for row in csv_reader:
        value = re.sub("\.[1-9]*", "", row[1])
        translate[row[0]] = value

##remove instances of the dictionary where the value is blank, i.e. there is no
##equivalent gene in the document 
translate_copy = translate.copy()
for k,v in translate_copy.items():
    if v == '':
       del translate[k]

##populate lists with the arbitrary reference for each differential regulation
##condition T1 and T2 are defined in the paper, up and down refers to the
##direction of the differential expression for that gene
T1up = []
T1down = []
T2up = []
T2down = []
ref_list = [T1up, T1down, T2up, T2down]

with open("regulation.csv") as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=",")
    for row in csv_reader:
        T1up.append(row[0])
        T1down.append(row[1])
        T2up.append(row[2])
        T2down.append(row[3])

##remove blanks from the lists
def remove_blank(str_list):
    while '' in str_list:
        str_list.remove('')
    return str_list

T1up = remove_blank(T1up)
T1down = remove_blank(T1down)
T2up = remove_blank(T2up)
T2down = remove_blank(T2down)

##populate lists with gene names by using the previously defined translate
##dictionary to convert the arbitrary reference lists
T1upGenes = []
T1downGenes = []
T2upGenes = []
T2downGenes = []
gene_list = [T1upGenes, T1downGenes, T2upGenes, T2downGenes]

notInKeyList = []

def ref_to_gene(str_list, new_list):
    for arbitrary_ref in str_list[1:]:
        if arbitrary_ref in translate:
            new_list.append(translate[arbitrary_ref])
        else:
            notInKeyList.append(arbitrary_ref)
    return new_list

for n in range(0,4):
    ref_to_gene(ref_list[n], gene_list[n])

##pad lists with white space so they're the same length
list_lengths = [(len(T1upGenes)),
                (len(T1downGenes)),
                (len(T2upGenes)),
                (len(T2downGenes))]

longest_list = sorted(list_lengths, reverse=True)[0]

def padd(str_list, til=longest_list):
    while len(str_list) < til:
        str_list.append("")
    return str_list

for n in gene_list:
    n = padd(n)

##write gene names to file "differential_expression_IGWSC.26.csv"
with open('differential_expression_IGWSC.26.csv', 'w', newline='') as file:
    writer = csv.writer(file)
    writer.writerow(["T1up", "T1down", "T2up", "T2down"])
    
    for n in range(0, longest_list):
        row = [T1upGenes[n], T1downGenes[n], T2upGenes[n], T2downGenes[n]]
        writer.writerow(row)

##to access other output
##T[1/2][up/down] contains the arbitrary reference
##T[1/2][up/down]Gene contains the translated genes
##notInKeyList contains the arbitrary references unable to be translated



