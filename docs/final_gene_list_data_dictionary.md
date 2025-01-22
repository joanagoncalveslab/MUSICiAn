# Data Dictionary

The excel file final_gene_list.xlsx contains 4 sheets:
- final_gene_list
- final_mutational_profile_T1
- final_mutational_profile_T2
- final_mutational_profile_T3
- config (this sheet can be ignored, just included as a check)

T1, T2, and T3 represent the 3 target sites: LV581, LV582, and LV583 respectively. 

| Sample            |      LUMC notation      |      Delft notation      | 
|-------------------|-----------------------|-----------------------|
| MB01  | LV581 | T1 |
| MB02  | LV581 | T1 |
| MB03  | LV582 | T2 |
| MB04  | LV582 | T2 |
| MB05  | LV583 | T3 |
| MB06  | LV583 | T3 |

Below is a description of the sheets and each of the columns contained therein.

---

## final_gene_list

This file contains the gene list. The columns are divided into four groups; T1, T2, T3, and Overall. Groups T1-T3 each has a single column, "Score". This is the Mahalanobis distance value, a measure of how far away a mutational profile is from the center of the distribution of all mutational profiles. In more practical terms, it is a measure of "outlierness".

The group "Overall" has the following columns:

| Column            |      Description      | 
|-------------------|-----------------------|
| Combined score    | Unified score for "outlierness". |
| Rank              | Outlier ranking based on the combined score, 1 being the top outlying gene. |
| Pseudo-control | Was this gene used as a pseudo-control for visualisation? True / False |

**Note:** Pseudo-controls were chosen as the 100 most centrally located mutational profiles relative to the overall distribution of all mutational profiles.

**Note:** You can use the overall ranking provided to select the number of genes of interest, such as the top 100.

---



## final_mutational_profile_<target_site>

This file contains the mutational profiles used for outlier detection, per target site. The first column contains the gene name, and every column after that describes the frequency of that mutational outcome. The categories are:

| Mutational Outcome             |      Note             | 
|--------------------------------|-----------------------|
| Any Insertion                  |                       |
| Deletion 0bp microhomology     |                       |
| Deletion 1bp microhomology     |                       |    
| Deletion 2bp microhomology     |                       |
| Deletion 3+bp microhomology    |                       |
| Deletion with insertion        | Includes templated insertions, tandem duplications, and compound tandem duplications.  |
| Homology Directed Repair       |                       |

**Reminder:** Very briefly, after filtering barcodes, we take an average mutational profile across all remaining barcodes, per gene, per replicate. We then take the average mutational profile for a gene between two replicates of the same target site as our final mutational profile, resulting in three mutational profiles per gene (one per target site). 