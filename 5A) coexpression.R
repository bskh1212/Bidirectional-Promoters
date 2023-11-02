##Coexpression Data

##note the BIP list has been refined since the data was downloaded, so there are
##superfluous entries and they are not in order nor necessarily in pairs

##import required libraries
library(Hmisc) #rcorr()
library(rlist) #t()

##set working directory
setwd("~/Leeds/Project/6. Co-expression")

##read in BIP genes
BIP <- read.delim("../3. Orthology/BIP_Orthologue_Information.csv",
                  sep=",",
                  header=T)
##keep genes IDs, pair IDs, description, strand, and distance
BIP <- data.frame(BIP$Gene_ID, BIP$PairID, BIP$Gene_Description, BIP$Strand,
                  BIP$Distance)

##rename columns
colnames(BIP) <- c("Gene_ID", "Pair_ID", "Description", "Strand", "Distance")

##read in proximal genes (control)
proximal <- read.delim("Proximal/proximal_triticum_aestivum.csv",
                       sep=",",
                       header=T)

####get proximal gene tpm data

##define empty data frame
tpm_data_proximal <- data.frame()
##get list of proximal tpm files
files <- dir(path="proximal_tpm")

##import a file to get variable information, i.e. the first 12 columns, and
##append it to tpm_data_proximal
variables <- read.delim("proximal_tpm/proximal_tpm_181-240.tsv", sep="\t",
                        header=T)
tpm_data_proximal <- rbind(tpm_data_proximal, variables[,1:12])

##loop through the files and append the tpm data only to tpm_data_proximal
for (n in files){
  file_data <- read.delim(paste("proximal_tpm/", n, sep=""), header=T, sep="\t")
  file_data <- file_data[,13:dim(file_data)[2]]
  tpm_data_proximal <- cbind(tpm_data_proximal, file_data)
}

##change column names
colnames(tpm_data_proximal)[2:12] <- c("Study", "High_Level_Tissue", "Tissue", 
                              "High_Level_Age", "Age",
                              "High_Level_Stress_Disease", "Stress_Disease",
                              "High_Level_Variety", "Variety", "Intermediate",
                              "Intermediate_Stress")

##create a copy of proximal data with only pairs where at least one has a value
##greater than 0.5

proximal$check <- TRUE

for (n in unique(proximal$pairs)){
  A <- proximal[proximal$Pair_ID==n, "Gene_ID"][1]
  B <- proximal[proximal$Pair_ID==n, "Gene_ID"][2]

  max_A <- max(tpm_data_proximal[A])
  max_B <- max(tpm_data_proximal[B])

  cutoff <- max(max_A, max_B)
  
  if(cutoff<0.5){
    proximal[proximal$Pair_ID==n, "check"] <- FALSE
  }
  
}

keep <- proximal[proximal$check==T, "Gene_ID"]
tpm_data_proximal_copy <- tpm_data_proximal[,colnames(tpm_data_proximal) %in% keep]

##order tpm data copy alphabetically
tpm_data_proximal_copy <- tpm_data_proximal_copy[, order(colnames(tpm_data_proximal_copy))]

##coexpression matrix of tpm data
##mat_proximal <- round(cor(tpm_data_proximal_copy),2)
##calculate the coexpression coefficients of all bidirectional genes against
##all other bidirectional genes
mat_proximal_all <- rcorr(as.matrix(tpm_data_proximal_copy),type="pearson")
mat_proximal <- mat_proximal_all$r
mat_proximal_p <- mat_proximal_all$P

##arrange proximal data so that one gene is in one column, the other gene
##in the second column, their respective descriptions, pair ID, correlation
##coefficient, and unadjusted p value in the final column

##pull those with odd row names
##get number of rows
rows <- nrow(proximal)
##extract odd rows 
odd_rows <- seq_len(rows) %% 2
# getting data from odd data frame
odd <- proximal[odd_rows == 1, c("Gene_ID", "Gene_Description", "Pair_ID",
                                 "Distance","Gene_End")]
colnames(odd)[1:2] <- c("Gene_2", "2_Description")
colnames(odd)[5] <- "Gene_End_1"

##pull those with even row names
##even = gene 2s
even <- proximal[odd_rows == 0, c("Gene_ID", "Gene_Description", "Pair_ID",
                                  "Strand", "Gene_Start")]
colnames(even)[1:2] <- c("Gene_1", "1_Description")
colnames(even)[5] <- "Gene_Start_2" 

##merge forward and reverse genes column ways
proximal_co <- merge(even, odd, by="Pair_ID") 

##define function that subsets dataframe "mat_proximal" (containing the 
##coefficients) by variables from proximal_co, which contains the proximal pairs
##on each row (and also mat_proximal_p, to pull the p values out)
subset_mat <- function(n){
  A <- proximal_co$Gene_1[n]
  B <- proximal_co$Gene_2[n]
  
  if (A %in% colnames(mat_proximal) && B %in% rownames(mat_proximal)){
    co <- mat_proximal[A,B]
    p_value <- mat_proximal_p[A,B]
    return(c(co, p_value))
  }
  else {return(NA)}
}

##use function with sapply to return list 
co_list <- sapply(1:dim(proximal_co)[1],subset_mat)

##bind list to proximal_co
co_list <- t(co_list)
proximal_co <- cbind(proximal_co, co_list)

##rename co list column
colnames(proximal_co)[10:11] <- c("Coexpression_Coefficient", "P_Value")

##########

####read tpm data files (BIP TA)

##define empty data frame
tpm_data <- data.frame()
##get list of tpm files
files <- dir(path="tpm")

##import a file to get variable information, i.e. the first 12 columns, and
##append it to tpm_data
variables <- read.delim("tpm/tpm_raw_1-60.tsv", sep="\t", header=T)
tpm_data <- rbind(tpm_data, variables[,1:12])

##loop through the files and append the tpm data only to tpm_data
for (n in files){
  file_data <- read.delim(paste("tpm/", n, sep=""), header=T, sep="\t")
  file_data <- file_data[,13:dim(file_data)[2]]
  tpm_data <- cbind(tpm_data, file_data)
}

##change column names
colnames(tpm_data)[2:12] <- c("Study", "High_Level_Tissue", "Tissue", 
                              "High_Level_Age", "Age",
                              "High_Level_Stress_Disease", "Stress_Disease",
                              "High_Level_Variety", "Variety", "Intermediate",
                              "Intermediate_Stress")

####filter the tpm data so only the BIP genes are taken fowards

##get tpm column names that are gene names but NOT in the BIP gene list
##start with all col names as a dataframe
tpm_colnames <- as.data.frame(colnames(tpm_data)[13:dim(tpm_data)[2]])
##rename the colnames of the dataframe
colnames(tpm_colnames) <- "Gene_ID"
##create a column with logical values indicating whether the gene is in the BIP
##data frame or not
tpm_colnames$logical <- tpm_colnames$Gene_ID %in% BIP$Gene_ID
##keep only the false ones
NOT_tpm_colnames <- tpm_colnames[tpm_colnames$logical==FALSE,]

##remove the genes in the NOT_tpm_columns dataframe (i.e. not BIP genes) from
##the tpm dataframe
for (n in NOT_tpm_colnames$Gene_ID){
  tpm_data[n] <- NULL
}

##create a copy of BIP data with only pairs where at least one has a value
##greater than 0.5

BIP$check <- TRUE

for (n in unique(BIP$pairs)){
  A <- BIP[BIP$Pair_ID==n, "Gene_ID"][1]
  B <- BIP[BIP$Pair_ID==n, "Gene_ID"][2]
  
  max_A <- max(tpm_data[A])
  max_B <- max(tpm_data[B])
  
  cutoff <- max(max_A, max_B)
  
  if(cutoff<0.5){
    BIP[BIP$Pair_ID==n, "check"] <- FALSE
  }
  
}

keep <- BIP[BIP$check==T, "Gene_ID"]
tpm_data_copy <- tpm_data[,colnames(tpm_data) %in% keep]

##order tpm data copy alphabetically
tpm_data_copy <- tpm_data_copy[, order(colnames(tpm_data_copy))]

##coexpression matrix of tpm data
mat_all <- rcorr(as.matrix(tpm_data_copy),type="pearson")
mat <- mat_all$r
mat_p <- mat_all$P

##arrange BIP data so that the forward gene is in one column, the reverse gene
##in the second column, their respective descriptions, pair ID, and correlation
##coefficient in the final column

##pull forward genes
forward <- BIP[BIP$Strand==1, c("Gene_ID", "Description", "Pair_ID")]
colnames(forward)[1:2] <- c("Forward", "F_Description") 

##pull reverse genes
reverse <- BIP[BIP$Strand==-1, c("Gene_ID", "Description","Pair_ID", "Distance")]
colnames(reverse)[1:2] <- c("Reverse", "R_Description")

##merge forward and reverse genes column ways
BIP_co <- merge(forward, reverse, by="Pair_ID") 

##define function that subsets dataframe "mat" (containing the coefficients) by
##variables from BIP_co, which contains the bidirectional pairs on each row
subset_mat <- function(n){
  A <- BIP_co$Forward[n]
  B <- BIP_co$Reverse[n]
  
  if (A %in% colnames(mat) && B %in% rownames(mat)){
    co <- mat[A,B]
    p_value <- mat_p[A,B]
    return(c(co, p_value))
  }
  else {return(NA)}
}

##use function with sapply to return list 
co_list <- sapply(1:dim(BIP_co)[1],subset_mat)

##bind list to BIP_co
co_list <- t(list.cbind(co_list))
BIP_co <- cbind(BIP_co, co_list)

##rename co list column
colnames(BIP_co)[7:8] <- c("Coexpression_Coefficient", "P_Value")

##Put all coexpression data in one dataframe
BIP_co_copy <- BIP_co
colnames(BIP_co_copy) <- c("Pair_ID", "Gene_1", "1_Description", "Gene_2",
                           "2_Description", "Distance",
                           "Coexpression_Coefficient", "P_Value")

BIP_co_copy$type <- "BIP"
proximal_co$type <- "proximal"

proximal_co_copy <- proximal_co
proximal_co_copy$Strand <- NULL
proximal_co_copy$Gene_Start_2 <- NULL
proximal_co_copy$Gene_End_1 <- NULL

all_data <- rbind(BIP_co_copy, proximal_co_copy)

##multiple test corrections
all_data$BH <- p.adjust(all_data$P_Value, method = "BH")


##all_data dimensions: 1519 - 10

##dataframe key
##all_data = coefficient data in pairs for BIP and proximal genes with corrected
##p values BIP_co = coefficient data in pairs, only BIP, uncorrected p values
##proximal_co = the same as above but for proximal instead of BIP
##tpm_data = data from wheat expression browser
##tpm_data_copy = data from wheat expression browser but with variables and zero
##value columns removed (unless their pair has expression)

##write to file
#write.csv(all_data, "Coexpression_Data.csv", row.names = F)
#write.csv(proximal_co, "Coexpression_Data_proximal_only.csv", row.names=F)
#write.csv(tpm_data, "tpm_data.csv", row.names = F)
#write.csv(tpm_data_proximal, "proximal_tpm_data.csv", row.names = F)
#write.csv(mat, "r_value_matrix.csv", row.names = F)
