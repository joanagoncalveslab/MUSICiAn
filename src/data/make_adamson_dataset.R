# install.packages(c("RSQLite", "tidyr", "dplyr", "Hmisc", "data.table"))

library(RSQLite)
library(data.table)
library(tidyr)
library(dplyr)

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

con = dbConnect(RSQLite::SQLite(), dbname = "data/processed/Adamson.db")


# list of files
listFile = c("target-1","target-2","target-3","target-4")
# currentDir = "E:\\Joost_Repair_Seq\\outSlurm\\"

# raw data directory
raw_data_file_dir = "data/raw/"


file = paste0(raw_data_file_dir,"20221205_hdr_run_rep_1_2_reducedColsR(1).txt")

dt = fread(file, sep = "\t", header=T, stringsAsFactors = FALSE)
rows = nrow(dt)

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



dbExecute(con, 'create index outcomes_full_idx on outcomes (Alias, Barcode, outcome, fraction_per_barcode, countEvents)')


dbDisconnect(conn = con)

## clean up
# dt = NULL
# dtBarcode = NULL
# dtGene = NULL
