##Bidirectionally Arranged Genes

##Use information from biomart (ensembl plants) to find head-to-head genes on
##opposite strands within 1000bp of one another (i.e. potential bidirectional 
##promoters)

##set working directory (change to location you've saved the script)
setwd("~/Leeds/Project/2. Head-to-Head Genes")

##define function
##input_file: should be a csv file with headers containing data frame with the 
##following columns in order: Gene_ID, Strand, Chromosome, Gene_Description,
##Gene_Start, Gene_End, TSS
##The order is important!
##return: function will return a dataframe of BIP genes with all the original
##columns plus Distance and Pair_ID. Distance is between the gene and the
##following gene and will be NA if the following gene is not paired
##bidirectionally with that following gene

locate_BIP <- function(input_file){
  ##import biomart output
  gene_data <- read.delim(input_file,
                          sep = ",",
                          header=T)
  ##rename columns
  colnames(gene_data) <- c("Gene_ID", "Strand", "Chromosome", "Gene_Description",
                           "Gene_Start", "Gene_End", "TSS")
  
  ##remove genes on an unknown (Un) chromosome from the list
  gene_data <- gene_data[!grepl("Un", gene_data$Chromosome),]
    
  ##remove genes from Mt and Pt chromosome from the list
  gene_data <- gene_data[gene_data$Chromosome != "Mt",]
  gene_data <- gene_data[gene_data$Chromosome != "Pt",]
  
  ##order the data first by TSS and then by chromosome
  gene_data <- gene_data[order(gene_data$TSS),]
  gene_data <- gene_data[order(gene_data$Chromosome),]
  
  ##order the data first by TSS and then by chromosome
  gene_data <- gene_data[order(gene_data$TSS),]
  gene_data <- gene_data[order(gene_data$Chromosome),]
  
  ##Add columns for distance and pair IDs
  gene_data$Distance <- NA
  gene_data$PairID <- NA
  
  ##assign input data frame length to len and count to 0
  len <- dim(gene_data)[1]
  count <- 0
  
  ##initialise empty dataframe in which to append putative gene pairs
  putative_bidirectional <- data.frame()
  
  ##loop through gene_data and compare the TSS to the TSS of the next gene in
  ##the data frame.
  for (row in 1:(len-1)){
    
    ##assign variables
    A <- gene_data$TSS[row]
    B <- gene_data$TSS[row+1]
    strandA <- gene_data$Strand[row]
    strandB <- gene_data$Strand[row+1]
    chromeA <- gene_data$Chromosome[row]
    chromeB <- gene_data$Chromosome[row+1]
    dist <- abs(A-B)
    
    ##update distance column
    gene_data$Distance[row] <- dist
    
    ##criteria:
    ##less than 1000bp apart
    ##on opposite strands 
    ##further than 100bp apart
    ##on the same chromosome
    ##the first in the pair is on the -1 strand 
    ##(and thus excludes overlapping genes)
    if (dist < 1000 & strandA != strandB & dist > 100 & chromeA == chromeB &
        strandA == -1){
      
      ##increase count
      count <- count + 1
      
      ##append bidirectional pair to dataframe
      putative_bidirectional <- rbind(putative_bidirectional, gene_data[row,])
      putative_bidirectional <- rbind(putative_bidirectional, gene_data[row+1,])
      
      ##assign the index of the last entry in the bidirectional gene list, which
      ##was just added
      index <- as.integer(dim(putative_bidirectional)[1])
      
      ##use that index to update the pair ID
      padded_count <- formatC(count, width = 4, format = "d", flag = "0")
      putative_bidirectional$PairID[index-1] <- paste("Pair", as.character(padded_count), sep="_")
      putative_bidirectional$PairID[index] <- paste("Pair", as.character(padded_count), sep="_")
    }
  }
  putative_bidirectional$Distance <- putative_bidirectional$Distance - 1
  return(putative_bidirectional)
}

##Triticum aestivum IWGSC (2100 genes)
putative_bidirectional_TA <- locate_BIP("Protein_Coding/triticum_aestivum_protein_coding_info.txt")
write.csv(putative_bidirectional_TA,
          "BIPs/putative_BIP_triticum_aestivum.csv",
          row.names = FALSE)

##Arabidopsis thanliana (4632 genes)
putative_bidirectional_AT <- locate_BIP("Protein_Coding/arabidopsis_thaliana_protein_coding_info.txt")
write.csv(putative_bidirectional_AT,
          "BIPs/putative_BIP_arabidopsis_thaliana.csv",
          row.names = FALSE)

##Oryza sativa japonica (2008 genes)
putative_bidirectional_OS <- locate_BIP("Protein_Coding/oryza_sativa_protein_coding_info.txt")
write.csv(putative_bidirectional_OS,
          "BIPs/putative_BIP_oryza_sativa.csv",
          row.names = FALSE)


##Oat (1592 genes) 
putative_bidirectional_oat <- locate_BIP("Protein_Coding/oat_protein_coding_info.txt")
write.csv(putative_bidirectional_oat,
          "BIPs/putative_BIP_oat.csv",
          row.names = FALSE)

##Maize (1362 genes) 
putative_bidirectional_maize <- locate_BIP("Protein_Coding/maize_protein_coding_info.txt")
write.csv(putative_bidirectional_maize,
          "BIPs/putative_BIP_maize.csv",
          row.names = FALSE)
table(putative_bidirectional_maize$Chromosome)

##Cherry Tomato (1106 genes) 
putative_bidirectional_tomato <- locate_BIP("Protein_Coding/cherry_tomato_protein_coding_info.txt")
write.csv(putative_bidirectional_tomato,
          "BIPs/putative_BIP_cherry_tomato.csv",
          row.names = FALSE)
table(putative_bidirectional_tomato$Chromosome)

##Black Cottonwood (1468 genes) 
putative_bidirectional_cottonwood <- locate_BIP("Protein_Coding/black_cottonwood_protein_coding_info.txt")
write.csv(putative_bidirectional_cottonwood,
          "BIPs/putative_BIP_black_cottonwood.csv",
          row.names = FALSE)
table(putative_bidirectional_cottonwood$Chromosome)

##Field Mustard (3412 genes) 
putative_bidirectional_mustard <- locate_BIP("Protein_Coding/field_mustard_protein_coding_info.txt")
write.csv(putative_bidirectional_mustard,
          "BIPs/putative_BIP_field_mustard.csv",
          row.names = FALSE)
table(putative_bidirectional_mustard$Chromosome)

##For a different organism:
custom <- locate_BIP("path/input_file_name.csv")
write.csv(putative_bidirectional_mustard,
          "name_of_output_file.csv",
          row.names = FALSE)
table(custom$Chromosome)

