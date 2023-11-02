##In Proximity but Not Bidirectionally Arranged
##Adapted from "Bidirectionally Arranged Genes" code

##Use information from biomart (ensembl plants) to find proximal genes on
##opposite strands within 1000bp of one another (i.e. potential bidirectional 
##promoters)

##set working directory
setwd("~/Leeds/Project/6. Co-expression/Proximal")

##define function
##input_file: should be a csv file with headers containing data frame with the 
##following columns in order: Gene_ID, Strand, Chromosome, Gene_Description,
##Gene_Start, Gene_End, TSS
##return: function will return a dataframe of proximal genes with all the 
##original columns plus Distance. Distance is between the gene and the following
##gene and will be NA if the following gene is not within the 100-1000bp range

locate_proximal <- function(input_file){
  ##import biomart output
  gene_data <- read.delim(input_file,
                          sep = ",",
                          header=T)
  ##rename columns
  colnames(gene_data) <- c("Gene_ID", "Strand", "Chromosome", "Gene_Description",
                           "Gene_Start", "Gene_End", "TSS")
  
  ##remove genes on an unknown (Un) chromosome from the list
  gene_data <- gene_data[gene_data$Chromosome != "Un",] 
  ##remove genes from Mt and Pt chromosome from the list
  gene_data <- gene_data[gene_data$Chromosome != "Mt",]
  gene_data <- gene_data[gene_data$Chromosome != "Pt",]
  
  ##order the data first by TSS and then by chromosome
  gene_data <- gene_data[order(gene_data$TSS),]
  gene_data <- gene_data[order(gene_data$Chromosome),]
  
  ##order the data first by TSS and then by chromosome
  gene_data <- gene_data[order(gene_data$TSS),]
  gene_data <- gene_data[order(gene_data$Chromosome),]
  
  ##Add column for distance and PairID
  gene_data$Distance <- NA
  gene_data$PairID <- NA
  
  ##assign input data frame length to len
  len <- dim(gene_data)[1]
  
  ##initialise count
  count <- 0

  ##initialise empty dataframe in which to append proximal gene pairs
  proximal <- data.frame()
  
  
  ##loop through gene_data and compare the TSS to the TSS of the next gene in
  ##the data frame
  for (row in 1:(len-1)){
    
    ##assign variables
    A <- gene_data$TSS[row]
    B <- gene_data$TSS[row+1]
    nameA <- gene_data$Gene_ID[row]
    nameB <- gene_data$Gene_ID[row+1]
    strandA <- gene_data$Strand[row]
    strandB <- gene_data$Strand[row+1]
    chromeA <- gene_data$Chromosome[row]
    chromeB <- gene_data$Chromosome[row+1]
    dist <- abs(A-B)
    
    ##update distance column
    gene_data$Distance[row] <- dist
    
    ##criteria:
    ##less than 1000bp apart
    ##on the same strand
    ##further than 100bp apart
    ##on the same chromosome
    ##on the same strand 
    ##there's multiple results for some genes so discount them as pairs if they're
    ##the same gene
    if (dist < 1000 & strandA == strandB & dist > 100 & chromeA == chromeB 
        & nameA != nameB){
      
      ##increase count
      count <- count + 1
      
      ##append to dataframe
      proximal <- rbind(proximal, gene_data[row,])
      proximal <- rbind(proximal, gene_data[row+1,])
      
      ##assign the index of the last entry in the bidirectional gene list, which
      ##was just added
      index <- as.integer(dim(proximal)[1])
      
      ##use that index to update the pair ID
      proximal$PairID[index-1] <- paste("Pair", as.character(count), sep="_")
      proximal$PairID[index] <- paste("Pair", as.character(count), sep="_")

    }
  }
  proximal$Distance <- proximal$Distance - 1
  return(proximal)
}

##Triticum aestivum IWGSC (3060 pairs)
proximal_TA <- locate_proximal("../../2. Head-to-Head Genes/Protein_Coding/triticum_aestivum_protein_coding_info.txt")
##write.csv(proximal_TA,
##          "proximal_triticum_aestivum.csv",
##          row.names = FALSE)

##
unique_proximal_TA <- unique(proximal_TA$Gene_ID)
##write.csv(unique_proximal_TA,
##          "unique_proximal_TA_list.csv",
##          row.names = FALSE)

##Arabidopsis thanliana ()
##putative_bidirectional_AT <- locate_BIP("Protein_Coding/protein_coding_arabidopsis_thaliana.txt")
##write.csv(putative_bidirectional_AT,
##          "BIPs/putative_BIP_arabidopsis_thaliana.csv",
##          row.names = FALSE)

##Oryza sativa japonica ()
##putative_bidirectional_OS <- locate_BIP("Protein_Coding/protein_coding_oryza_sativa.txt")
##write.csv(putative_bidirectional_OS,
##          "BIPs/putative_BIP_oryza_sativa.csv",
##          row.names = FALSE)
