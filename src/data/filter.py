import pandas as pd
import numpy as np
from tqdm import tqdm
from src.config import get_interim_dir, FILTER_COUNT, CORR_THRESHOLD, BAD_TARGET_THRESHOLD, NUM_GRNAS_THRESHOLD
from src.data.load_dataset import get_target_site


READ_COUNT_THRESHOLD = FILTER_COUNT

# accepts a list of gRNAs with multi-index {Alias, Gene, Barocde} and repair outcomes read counts
# filters high quality gRNAs
def filter_barcodes(profiles):

    print("Start filtering...")

    print("Before: {} gRNAs, {} genes".format(profiles.shape[0], profiles.index.unique(level="Gene").shape[0]))

    gene_results = pd.DataFrame(0, index=profiles.index, columns=["Filtered"])
    gene_results["Counts"] = profiles.sum(axis=1)

    ######################################
    ## Filter by read counts
    ######################################
    read_counts_less = profiles.sum(axis=1) < READ_COUNT_THRESHOLD
    gene_results.loc[read_counts_less[read_counts_less].index, "Filtered"] = 1
    profiles = profiles.loc[gene_results[gene_results["Filtered"] == 0].index, :]
    print("After {} read count cutoff: {} gRNAs, {} genes".format(READ_COUNT_THRESHOLD, gene_results[gene_results["Filtered"] == 0].shape[0], gene_results[gene_results["Filtered"] == 0].index.unique(level="Gene").shape[0]))


    ######################################
    ## Filter by gRNA counts per gene
    ######################################
    guide_counts_per_gene = profiles.groupby(["Target", "Gene"]).size()
    guide_counts_per_gene = guide_counts_per_gene[guide_counts_per_gene < NUM_GRNAS_THRESHOLD]
    for idx in guide_counts_per_gene.index:
        islice = pd.IndexSlice
        gene_results.loc[islice[idx[0],:,idx[1]], "Filtered"] = 2
    profiles = profiles.loc[gene_results[gene_results["Filtered"] == 0].index, :]
    print("After {} plus gRNA cutoff: {} gRNAs, {} genes".format(NUM_GRNAS_THRESHOLD, gene_results[gene_results["Filtered"] == 0].shape[0], gene_results[gene_results["Filtered"] == 0].index.unique(level="Gene").shape[0]))


    ######################################
    ## Filter by median correlation gRNAs in the same replicate
    ######################################
    # idx = pd.IndexSlice
    # corr_results = []
    # for i in range(1, 7):
    #     sample = "MB0" + str(i)
    #     target = get_target_site(sample)
    #     genes = profiles.loc[(target, sample)].index.unique(level="Gene").to_series()
    #     for gene in tqdm(genes):
    #         # barcodes = profiles.loc[(sample, gene)].index.to_list()
    #         corr_ = profiles.loc[idx[target, sample, gene, :]].T.corr()
    #         np.fill_diagonal(corr_.values, np.nan)
    #         corr_ = corr_.median()
    #         r = corr_.reset_index()
    #         r["Target"] = target
    #         r["Alias"] = sample
    #         r["Gene"] = gene
    #         corr_results.append(r)
    # corr_results = pd.concat(corr_results)
    # corr_results = corr_results.set_index(["Target", "Alias", "Gene", "Barcode"])
    # corr_results.columns = ["Corr"]
    # gene_results = gene_results.join(corr_results)
    # gene_results.loc[corr_results[corr_results["Corr"].notnull() & (corr_results["Corr"] < CORR_THRESHOLD)].index, "Filtered"] = 3
    # profiles = profiles.loc[gene_results[gene_results["Filtered"] == 0].index, :]


    ######################################
    ## Filter by median correlation to other gRNAs for the same replicate
    ######################################
    idx = pd.IndexSlice
    corr_results_within = []
    corr_results_between = []
    for target in profiles.index.unique(level="Target"):
        # target = "T{}".format(i)
        genes = profiles.loc[target].index.unique(level="Gene").to_series()
        for gene in tqdm(genes):
            # barcodes = profiles.loc[(sample, gene)].index.to_list()
            corr_ = profiles.loc[idx[target, :, gene, :]].T.corr()
            np.fill_diagonal(corr_.values, np.nan)
            aliases = corr_.index.unique("Alias")
            for alias in aliases:
                    acorr_ = corr_.loc[alias, alias]
                    acorr_ = acorr_.median()
                    r  = acorr_.reset_index()
                    r["Target"] = target
                    r["Alias"] = alias
                    r["Gene"] = gene
                    corr_results_within.append(r)

    corr_results_within = pd.concat(corr_results_within)
    corr_results_within = corr_results_within.set_index(["Target", "Alias", "Gene", "Barcode"])
    corr_results_within.columns = ["Corr_Within"]
    gene_results = gene_results.join(corr_results_within)
    gene_results.loc[gene_results[gene_results["Corr_Within"].notnull() & (gene_results["Corr_Within"] < CORR_THRESHOLD)].index, "Filtered"] = 3
    profiles = profiles.loc[gene_results[gene_results["Filtered"] == 0].index, :]

    ######################################
    ## Filter by median correlation to other gRNAs for the other replicate
    ######################################
    idx = pd.IndexSlice
    corr_results_within = []
    corr_results_between = []
    for target in profiles.index.unique(level="Target"):
        genes = profiles.loc[target].index.unique(level="Gene").to_series()
        for gene in tqdm(genes):
            # barcodes = profiles.loc[(sample, gene)].index.to_list()
            corr_ = profiles.loc[idx[target, :, gene, :]].T.corr()
            np.fill_diagonal(corr_.values, np.nan)
            aliases = corr_.index.unique("Alias")
            if aliases.shape[0] == 2:
                a, b = aliases
                acorr_ = corr_.loc[b, a]
                acorr_ = acorr_.median()
                r  = acorr_.reset_index()
                r["Target"] = target
                r["Gene"] = gene
                r["Alias"] = a
                corr_results_between.append(r)

                acorr_ = corr_.loc[a, b]
                acorr_ = acorr_.median()
                r  = acorr_.reset_index()
                r["Target"] = target
                r["Gene"] = gene
                r["Alias"] = b
                corr_results_between.append(r)

    if len(corr_results_between) > 0:
        corr_results_between = pd.concat(corr_results_between)
        corr_results_between = corr_results_between.set_index(["Target", "Alias", "Gene", "Barcode"])
        corr_results_between.columns = ["Corr_Between"]
        gene_results = gene_results.join(corr_results_between)
        gene_results.loc[gene_results[gene_results["Corr_Between"].notnull() & (gene_results["Corr_Within"] >= CORR_THRESHOLD) & (gene_results["Corr_Between"] < CORR_THRESHOLD)].index, "Filtered"] = 4
        profiles = profiles.loc[gene_results[gene_results["Filtered"] == 0].index, :]


    # less_corr_threshold_gRNAs = gene_results[gene_results["Filtered"] == 3].groupby("Barcode").size() > 1
    # less_corr_threshold_gRNAs = less_corr_threshold_gRNAs[less_corr_threshold_gRNAs].index.to_numpy()
    # profiles = profiles.loc[~profiles.index.get_level_values("Barcode").isin(less_corr_threshold_gRNAs), :]
    # gene_results.loc[gene_results.index.get_level_values("Barcode").isin(less_corr_threshold_gRNAs) & (gene_results["Filtered"] == 3), "Filtered"] = 3.1
    # gene_results.loc[gene_results.index.get_level_values("Barcode").isin(less_corr_threshold_gRNAs) & (gene_results["Filtered"] == 0), "Filtered"] = 3.2
    # gene_results.loc[gene_results["Filtered"] == 3, "Filtered"] = 0    
    print("After {} median correlation cutoff: {} gRNAs, {} genes".format(CORR_THRESHOLD, gene_results[gene_results["Filtered"] == 0].shape[0], gene_results[gene_results["Filtered"] == 0].index.unique(level="Gene").shape[0]))


    # profiles = profiles.div(profiles.sum(axis=1), axis=0).groupby(["Target", "Gene"]).mean()
    # print("Pooling gRNAs for different replicates together")

    # find barcodes which have failed in all replicates for a target site
    barcode_counts = gene_results.groupby(["Barcode", "Target"]).size()
    failed_counts = gene_results[gene_results["Filtered"] != 0].groupby(["Barcode", "Target"]).size()
    barcode_counts = barcode_counts.to_frame(name="total").join(failed_counts.to_frame(name="failed")).fillna(0)
    barcode_counts = barcode_counts[barcode_counts["total"] == barcode_counts["failed"]].groupby("Barcode").size()
    # filter out any where all replicates have failed in two target sites
    failed_barcodes = barcode_counts[barcode_counts >= BAD_TARGET_THRESHOLD]
    failed_barcodes = failed_barcodes.index.to_numpy()
    profiles = profiles.loc[~profiles.index.get_level_values("Barcode").isin(failed_barcodes)]
    gene_results.loc[gene_results.index.get_level_values("Barcode").isin(failed_barcodes) & (gene_results["Filtered"] == 0), "Filtered"] = 5

    print("After filtering barcodes which were bad in {} targets: {} gRNAs, {} genes".format(BAD_TARGET_THRESHOLD, gene_results[gene_results["Filtered"] == 0].shape[0], gene_results[gene_results["Filtered"] == 0].index.unique(level="Gene").shape[0]))

    return gene_results, profiles


if __name__ == "__main__":
    p = pd.read_pickle(get_interim_dir() + "test.pkl")
    # p = p[p.index.get_level_values("Alias") != "MB01"]
    q, r = filter_barcodes(p)
    print("Done.")
