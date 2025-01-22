# install.packages(c("RSQLite", "tidyr", "dplyr", "Hmisc", "data.table"))

library(RSQLite)
library(data.table)
library(tidyr)
library(dplyr)
library(ggplot2)
library(Hmisc)
library(ggbeeswarm)


dtbinom <- function(fraction, trials, refFraction) {
  dt <- data.table(fraction , trials, refFraction)
  dt[, P := binom.test(x=fraction, n=trials, p = refFraction)$p.value, seq_len(nrow(dt))]$P
}

appendToDB = F

minReads = 1000
minMutReads = 200 ##changed to 1000 for the subscreen
filterReads = F
filterMutReads = T
minCountEvents = 2
removeSNVs = T
removeHDR_1MM = T
transformTDs = T
maxOutcomes = 100 ##just a test, change later
writeOutcomes = T ##should be TRUE
writeTornado = T
maxDistanceToCutSite = 2
Pvalue_inf = 300
writePValues = T


# list of files
# listFile = c("MBsub1_reseq_full.txt", "MBsub2_reseq_full.txt", "MBsub3_reseq_full.txt") #correct ones
listFile = c("LV581_rep1.txt", "LV581_rep2.txt", "LV581_rep3.txt",
    "LV582_rep1.txt", "LV582_rep2.txt", "LV582_rep3.txt",
    "LV583_rep1.txt", "LV583_rep2.txt", "LV583_rep3.txt") #correct ones

processedDir = "~/repos/MUSICian/data/processed/"
rawDir = "~/repos/MUSICian/data/raw/"
# db_name = "MBCrisprMBSubscreeen_Full_2.db"
db_name = "MBCrisprMBSubscreeen_TUD.db"

db_path = file.path(processedDir,db_name)

con = dbConnect(RSQLite::SQLite(), dbname = db_path)


for (f in listFile) {
  print(f)
  file = file.path(rawDir, f)
  dt = fread(file, sep = "\t", header=T, stringsAsFactors = FALSE, fill = T)
  
  ##the number of initial rows is use later
  rows = nrow(dt)
  
  ##remove some columns are they are not needed, saves memory
  #more columns can be added
  dt = dt %>% select(-Name, -Split, -Dir, -File, -Raw, -getSubjectComments, -getIDPart, -possibleDouble, -ClassName, -InZone, -jumpedLeft, -jumpedRight, -entireQueryUsed)

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
