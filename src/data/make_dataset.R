# install.packages(c("RSQLite", "tidyr", "dplyr", "Hmisc", "data.table"))

library(RSQLite)
library(data.table)
library(tidyr)
library(dplyr)


dtbinom <- function(fraction, trials, refFraction) {
  dt <- data.table(fraction , trials, refFraction)
  dt[, P := binom.test(x=fraction, n=trials, p = refFraction)$p.value, seq_len(nrow(dt))]$P
}

appendToDB = F

minReads = 1000
minMutReads = 200
filterReads = F
filterMutReads = F
minCountEvents = 2
removeSNVs = T
removeHDR_1MM = T
transformTDs = F
maxOutcomes = 20
writeOutcomes = TRUE
writeTornado = F
useGeneSubset = F
useRandomSubset = F
randomSubsetSize = 500
maxDistanceToCutSite = 2

if(useGeneSubset == TRUE){
  con = dbConnect(RSQLite::SQLite(), dbname = "data/processed/MBCrisprMBAgain_1.2Subset_cseale.db")
} else {
  con = dbConnect(RSQLite::SQLite(), dbname = "data/processed/MBCrisprMBAgain_1.2_cseale.db")
}

# list of files
listFile = c("MB01", "MB02", "MB03", "MB04", "MB05","MB06","MB07","MB08")
# just process one file for noe
# listFile = c("MB01")

# select DDR subset
GeneSubset = fread("src/data/GeneSubset.txt", header=T, sep = "\t", stringsAsFactors = F)

##actually some barcodes were not correct, they are correct in this file
correctBarcodeFile = "src/data/yusa_corrected.txt"
correctedBarcode = read.csv(correctBarcodeFile,sep = "\t", stringsAsFactors = F)
correctedBarcode = correctedBarcode %>% select(ID_2, ID, correctBarcode, correctGene) %>% rename(Barcode = ID)

##load the non-targetting controls
nontargettingFile = "src/data/nontarget.txt"
nontargetting = read.csv(nontargettingFile,sep = "\t", stringsAsFactors = F)
nontargetting = nontargetting %>% rename(ID_2 = id, Gene_NonT = Gene, Barcode_NonT = Barcode)

##overwrite the barcodes and genes that are used as non-targetting controls
correctedBarcode = merge(correctedBarcode, nontargetting, all.x = T, by="ID_2")
correctedBarcode = correctedBarcode %>% mutate(correctBarcode = ifelse(is.na(Barcode_NonT),correctBarcode, Barcode_NonT))
correctedBarcode = correctedBarcode %>% mutate(correctGene = ifelse(is.na(Gene_NonT),correctGene, Gene_NonT))
correctedBarcode = correctedBarcode %>% select(-Gene_NonT, -Barcode_NonT, -ID_2)
correctedBarcode = rbind(correctedBarcode,c("Empty-1","Empty-1","Empty"))
correctedBarcode = rbind(correctedBarcode,c("control","Empty-1","Empty"))


interm_non_targetting_gene_names_file = "data/interim/non-targeting_names.txt"
write.table(correctedBarcode %>% filter(correctGene == "Non-targeting"),file = interm_non_targetting_gene_names_file, row.names = F)

controlsFixed = fread("data/raw/MBcontrols_fixed_reducedCols.txt", sep = "\t", header=T, stringsAsFactors = FALSE)

randomGenes = NULL

# raw data directory
raw_data_file_dir = "data/raw/"

# process subset for debugging
row_limit = Inf

# process all data
# row_limit = 1000

for (f in listFile) {
  print(f)
  if(useGeneSubset == FALSE){
    file = paste0(raw_data_file_dir,f,"_reducedCols.txt")
    dt = fread(file, sep = "\t", header=T, stringsAsFactors = FALSE, nrows=row_limit)
  } else {
    file = paste0(raw_data_file_dir,f,"_subset_reducedCols.txt")
    dt = fread(file, sep = "\t", header=T, stringsAsFactors = FALSE, nrows=row_limit)
  }

  dt = rbind(dt, controlsFixed %>% filter(Alias == dt$Alias[1]))
  rows = nrow(dt)
  dt = merge(dt, correctedBarcode, by="Barcode")
  dt$Barcode = dt$correctBarcode
  dt$Gene = dt$correctGene
  dt = dt %>% select(-correctBarcode, -correctGene)
  
  ## Subset for genelist
  if(useGeneSubset == TRUE){
    if(useRandomSubset == TRUE){
      if(is.null(randomGenes)){
        randomGenes = dt %>% select(Gene) %>% distinct() %>% sample_n(size = randomSubsetSize)
      }
      genes = c(randomGenes$Gene,GeneSubset$Gene)
      dt = dt %>% filter(Gene %in% genes)
    } else{
      dt = dt %>% filter(Gene %in% GeneSubset$Gene)
    }
  }

  #remove non-overlapping events based on distance to relative cut site
  dt = dt %>% filter(delRelativeStartRight <= maxDistanceToCutSite & delRelativeEnd >= -maxDistanceToCutSite+1)
 
  if(removeSNVs){
    dt = dt %>% filter(Type != "SNV")
  }

  if(removeHDR_1MM){
    dt = dt %>% filter(Type != "HDR1MM")
  }

  if(transformTDs == T){
    ##remove TDs and put them as DELINS here
    dt = dt %>% mutate(Type = ifelse(Type == "TANDEMDUPLICATION" | Type == "TANDEMDUPLICATION_COMPOUND", "DELINS", Type))
  }

  ### remove sgrnas with too little events
  sgRNACount = dt %>% group_by(Barcode) %>% count(wt=countEvents)
  sgRNAMutCount = dt %>% filter(Type != "WT") %>% group_by(Barcode) %>% count(wt=countEvents)
  
  sgRNAsKeep = sgRNACount %>% filter(n>=minReads)
  sgRNAsKeep2 = sgRNAMutCount %>% filter(n>=minMutReads)
  
  if(filterReads == T & filterMutReads == F) {
    dt = dt %>% filter(Barcode %in% sgRNAsKeep$Barcode)
  } else if(filterReads == F & filterMutReads == T) {
    dt = dt %>% filter(Barcode %in% sgRNAsKeep2$Barcode)
  } else if(filterReads == T & filterMutReads == T){
    dt = dt %>% filter(Barcode %in% sgRNAsKeep$Barcode & Barcode %in% sgRNAsKeep2$Barcode)
  }

  ##add a countTable to the DB
  countTable = dt %>% group_by(Alias, Barcode, Gene) %>% summarise(allEvents = sum(countEvents))
  countTableNonWT = dt %>% filter(Type != "WT") %>% group_by(Barcode) %>% summarise(mutEvents = sum(countEvents))
  countTable = merge(countTable, countTableNonWT, by = "Barcode")
  
  dbWriteTable(conn = con, name = "countTable", countTable, append = appendToDB, overwrite = !appendToDB)  
  
  dt$SubType = dt$Type

  dt$outcome <- paste(dt$Type,dt$delRelativeStart,dt$delRelativeEnd,dt$insert,paste0(dt$homologyLength,"bp"),sep = "|")
  currentRows = nrow(dt)

  print(paste("Initial rows", rows," reduced to ",currentRows, "fraction:",currentRows/rows))


  if (writeOutcomes) {
    print("Preparing outcomes table")

    dtOutcomes = dt %>% group_by(Alias, Barcode) %>% mutate(fraction_per_barcode = fraction/sum(fraction))
    outcomes = dtOutcomes %>% select(Alias, Gene, Barcode, outcome, fraction_per_barcode, Type, countEvents, delRelativeStart, delRelativeEnd, delSize, insSize, homologyLength, insertion)
    test = merge(outcomes, countTable %>% select(-Alias, -Gene), by = "Barcode")

    # We will ignore these p-value tests for now, may come back to these later
    # test = merge(outcomes, countTableGene, by = "Gene")
    # test = test %>% mutate(trials = ifelse(outcomeTop == "WT|0|0||-1bp", countEvents, mutEvents))
    # test = test %>% mutate(counts = floor(fraction*trials))
    # test = test %>% group_by(outcomeTop) %>% mutate(mean = mean(fraction))
    # test$Pvalue = -log10(dtbinom(test$counts, test$trials, test$mean))
    # test = test %>% mutate(Pvalue = ifelse(Pvalue == Inf, 300, Pvalue))
    
    print(paste("Writing outcomes gene",nrow(outcomes)))
    dbWriteTable(conn = con, name = "outcomes", test, append = appendToDB, overwrite = !appendToDB)
    print("Writing done genes")
  }
  appendToDB = T
}


dbExecute(con, 'create index outcomes_full_idx on outcomes (Alias, Gene, Barcode, outcome, fraction_per_barcode, countEvents)')


dbDisconnect(conn = con)

## clean up
dt = NULL
dtBarcode = NULL
dtGene = NULL
