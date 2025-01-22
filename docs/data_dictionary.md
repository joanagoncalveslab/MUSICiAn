# MUSIC Data Dictionary

## Storage

Data stored locally and on OneDrive at:

```/Users/colm/OneDrive - Delft University of Technology/My Drafts/Predicting Repair Pathways by Mutational Signatures/transfer_1687439_files_dacda79d/```


Or on the U: drive under:

```/tudelft.net/staff-umbrella/protonddr/MUSICian/data```

## Samples

 There are 8 samples, numbered 01 through 08. The below table shows the sample details:

| Sample   |      Target Site      |  Genotype             |
|----------|-----------------------|-----------------------|
| 01, 02   | 01                    | WT                    |
| 03, 04   | 02                    | WT                    |
| 05, 06   | 03                    | WT                    |
| 07       | 01                    | NHEJ-deficient        |
| 08       | 01                    | POLQ-deficient        |

## Raw data

Data was generated using the SIQ tool [1]. The user guide for the tool can be found at https://github.com/RobinVanSchendel/SIQ/blob/master/SIQ%20and%20SIQPlotteR%20User%20guide.docx. 

Raw data stored as tab seperated files `MBXX_reductedCols.txt`, where `XX` is the number of the sample. The columns are described below:


| Column        |      Description      | 
|---------------|-----------------------|
| countEvents   | Number of reads assigned to this outcome |
| fraction      | Fraction of total reads assigned to this outcome |
| Barcode       | ??? |
| Gene          | ??? |
| Alias         | Sample Name (e.g MB01) |
| Subject       | ??? |
| del           | Deleted nucleotide sequence |
| insertion     | Inserted nucleotide sequence |
| delStart      | Starting location of the deletion relative to the start of the target sequence|
| delEnd        | Ending location of the deletion relative to the start of the target sequence |
| delRelativeStart | Starting location of the deletion relative to the cut site |
| delRelativeEnd | Ending location of the deletion relative to the cut site |
| delRelativeStartRight | Same as delRelativeStart |
| delRelativeEndRight | Same as delRelativeEnd |
| delRelativeStartTD | Used by SIQPlotter |
| delRelativeEndTD | Used by SIQPlotter |
| getHomologyColor | Used by SIQPlotter |
| homology | Homologous sequence flanking the cut site |
| homologyLength | Size of homologous sequence |
| delSize | Size of deletion |
| insSize | Size of insertion |
| SNVMutation | Describes the type of SNV (e.g GC->TA) |
| Type | See mutation type table below |


## Mutation Types

The mutation types are described below:
| Type        |      Description      | 
|-------------|-----------------------|
| DELETION    | Simple deletion  |
| INSERTION   | Simple insertion |
| TINS        | Templated Insertion (Deletion with insertion copied from flanks) |
| HDR | Homology directed repair |
| HDR1MM      | HDR with 1 mismatch, usually caused by sequencing errors |
| TANDEMDUPLICATION | Tandem Duplication, *double check definition with Robin* |
| TANDEMDUPLICATION_COMPOUND | Tandem Duplication with additional insertion |
| SNV | Single nucleotide variation |
| WT | Wild-type |

Note: HDR and HDR_1MM events only called with HDR sequence is provided.



## SQLite Database
You can connect to the SQLite DB easily using `TablePlus`. See https://tableplus.com/. Some data has already been preprocessed and is stored under `MBCrisprMBAgain_1.2.db` and `MBCrisprMBAgain_12Subset.db`. *Need to check what has been done here.*

There are 8 tables: 
[('countTable',), ('geneAlt',), ('barcodeAlt',), ('barcodecount',), ('outcomes',), ('topoutcomes',), ('outcomesGene',), ('topoutcomesGene',)]

The processed data is further subdivided into multiple mutation types:
| outcomeTop | Description |
| -----------|-------------|
| DELETION | A deletion |
| LARGE_DELETION | A deletion > than ? nucleotides |
| HDR | HDR event |
| INSERTION_1bp | Insertion of 1bp |
| PQ_DELETION | Specific deletion associated with POLQ knockouts conducted in previous experiments. Can be treated as a regular deletion |
| TINS | Templated insertion |
| WT | Wild-type |
| DELINS | Deletion with insertion |


## Other Datasets

### GeneSubset2
Manually selected genes from the literature/experience and from previous screens which were focused on the Fanconi anemia pathway. Not all these genes are positive controls, for example, Brca2 did not seem to have a strong effect.

### GeneSubsetSD30 / 2023-01-06_CustomLib

Method of analysis:
- Default QC setting for reads
- Split data per barcode, per alias (3 target sites/2 replicates per site = 6 aliases)
- At least 2 reads of barcode + alias required to support an outcome
 
At least 200 mutagenic (non-WT) reads per barcode + alias to maintain the barcode for that alias > assume the profile is “reliable”
- Outcomes are aggregated into types > “INSERTION” “INSERTION_1BP” “TINS” “WT” “LARGE_DELETION” (>14BP) “PQ_DELETION” (manual attributed) “HDR” “DELINS” “DELETION”
- Per Alias, the mutation profiles/types from all barcodes targeting the same gene are averaged to generate a gene profile > each barcode has an equal weight irrespective of #reads
- For each Alias, the Types’ “average” fraction over the full dataset is compared to the fraction found per gene, to calculate a log-fold-change (LFC)
- For each Alias + gene, the #reads is used to calculate the probability of finding the observed fraction, to come to a P-value
- Per Type, the mean LFC over all aliases is calculated and ranked > each alias has an equal weight. The ‘logic’ behind this is that “noise” will have positive and negative LFC which cancel each other out.
- Per Type, the SD of the LFC is determined
 
Gene hit selection:
- At least 3 barcodes per gene need to survive the filters for the gene to be included
- WT fraction can be no more than 80%, as higher WT% seems to correlate with toxicity
- A genes’ meanLFC needs to exceed at least 2.5x SD of the Type, AND;
- A genes + type’ meanLFC needs to exceed at least 0.25 (2LOG) from the full data set, AND;
- A genes + type’ meanPVal needs to exceed at least 7.5, to be called a hit
- Per type, the top20 min-LFC and max-LFC genes are always included, regardless of the above settings;
- Per type, at most 150 genes are included, ranked per meanPvalue
 
I made a table “GeneHitListFull”, where you can find a list of genes (first column). Then in the other columns several data sets:
- Yusa = Yusa-library (original genome wide)
- GeneHits = gene hits from MUSIC, based on the above settings
- AdamsongLong = The first, ‘full’ library in Adamson et al.
- AdamsonShort = The second, subselection, library in Adamson et al.
- AdamsonAct = A third, ‘100 most active barcodes’, library in Adamson et al.
- Kosicki = The genes that were tested in Kosicki et al.


Data Description:
| Name       | Description |
| -----------|-------------|
| Gene | - |
| Barcode | - |
| Sequence | - |
| Oligo | - |
| In_GeneHits | | 
| In_Yusa | Group that made the original "genome wide" library, but few genes were missing, including Rad54 |
| In_Adamson | Genes common with the Adamson dataset |
| In_Kosicki | Another paper (ask Marco for reference) |
| In_JMList | Genes a postdoc in Marcos group was interested in |
| In_RBList | Genes Robin was interested in |
| BcHits_fromYusa | - |
| BcRest_fromYusa | - |
| Bclast_fromYusa | - |
| BcNT_fromGecko | Non-targetting library from Gecko |



### GeneLookupTable
File with synonyms for different gene names. Not perfect, manually curated by Marco. Can also look at genecards for aliases. Can also check the genomic location of the genes.


# References

[1] van Schendel, R., Schimmel, J., & Tijsterman, M. (2022). SIQ: easy quantitative measurement of mutation profiles in sequencing data. bioRxiv.