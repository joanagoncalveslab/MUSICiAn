import pandas as pd
import numpy as np
from tqdm import tqdm

from src.data.load_dataset import load_lumc_outcomes
from src.data.categories import aggregate_outcomes
from src.config import get_common_barcodes

def threshold_read_counts(profiles):
    profiles = profiles.iloc[:,:7].copy()
    column_name = "Read Count"
    read_counts_less = profiles.sum(axis=1) 
    profiles[column_name] = read_counts_less
    return profiles[column_name]

def filter_by_median_correlation_to_other_gRNAs(profiles):
    profiles = profiles.iloc[:,:7].copy()

    idx = pd.IndexSlice
    corr_results_within = []
    column_name = "corr w/ other gRNAs"

    for target in profiles.index.unique(level="Target"):
        genes = profiles.loc[target].index.unique(level="Gene").to_series()
        for gene in tqdm(genes):
            corr_ = profiles.loc[idx[target, :, gene, :]].iloc[:,:7].T.corr()
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
    corr_results_within.columns = [column_name]
    profiles = profiles.join(corr_results_within)
    return profiles[column_name]

def filter_by_median_correlation_to_other_gRNAs_in_paired_replicate(profiles):
    profiles = profiles.iloc[:,:7].copy()

    idx = pd.IndexSlice
    corr_results_between = []
    column_name = "corr w/ other gRNAs in paired replicate"

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
        corr_results_between.columns = [column_name]
        profiles = profiles.join(corr_results_between)
    return profiles[column_name]


if __name__ == "__main__":

    ALIASES = ["MB0{}".format(i) for i in range(1, 7)]
    TARGETS = ["T{}".format(i) for i in range(1, 4)]
    LIMIT = None

    # takes about 6 minutes to load the data into memory
    p = {}
    for alias in ALIASES:
        df = load_lumc_outcomes(limit=LIMIT, sample_name=alias)
        # convert to profiles
        df = df.reset_index().pivot(index=["Target", "Alias", "Gene", "Barcode"], columns=["lumc_category", "del_len", "start"], values="countEvents").fillna(0)
        p[alias] = df
    # p = pd.concat(p, axis=0).fillna(0)

    # aggregate categories
    idx = pd.IndexSlice
    Xs = []
    for alias in tqdm(ALIASES):
        # X = p.loc[p.index.get_level_values("Alias") == alias, :]
        X = p[alias]
        X = aggregate_outcomes(X, alias)
        Xs.append(X)
    profiles = pd.concat(Xs).fillna(0)

    passed_read_counts = threshold_read_counts(profiles)
    passed_read_counts

    passed_within_corr = filter_by_median_correlation_to_other_gRNAs(profiles)
    passed_within_corr

    passed_between_corr = filter_by_median_correlation_to_other_gRNAs_in_paired_replicate(profiles)
    passed_between_corr

    pd.concat((passed_read_counts, passed_within_corr, passed_between_corr), axis=1).to_csv("candidates_sgRNA_analysis.txt", sep="\t")