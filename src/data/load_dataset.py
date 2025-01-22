from random import sample
import sqlite3
from src.config import get_db_file, get_adamson_db_file, get_validation_db_file
import pandas as pd

MINIMUM_MUTATED_EVENTS = None
MINIMUM_FREQUENCY = None

def load_lumc_outcomes(includeWT=False, sample_name=None, limit=None, validation=False, downsample=None):
    query_str = "select Target, Alias, Gene, Barcode, fraction_per_barcode, lumc_category, del_len, start, countEvents from {}"
    con = sqlite3.connect(get_db_file() if not validation else get_validation_db_file())
    
    table_name = "lumc_outcomes"

    if validation:
        table_name = table_name + "_" + downsample
        
    query_str = query_str.format(table_name)    

    if sample_name is not None:
        query_str += " where Alias is '{}'".format(sample_name)
    if limit is not None:
        query_str += " limit {}".format(limit)
    df = pd.read_sql_query(query_str, con)
    df = df.set_index(["Target", "Alias", "Gene", "Barcode"])
    if not includeWT:
        df = df[df["lumc_category"] != "Wild-Type"]
        df = normalise_dataset_per_barcode(df)
    return df

def load_adamson_lumc_outcomes(includeWT=False, sample_name=None, limit=None):
    query_str = "select Target, Alias, Gene, Transcript, Barcode, fraction_per_barcode, lumc_category, del_len, start, countEvents from lumc_outcomes"
    con = sqlite3.connect(get_adamson_db_file())
    if sample_name is not None:
        query_str += " where Alias is '{}'".format(sample_name)
    if limit is not None:
        query_str += " limit {}".format(limit)
    df = pd.read_sql_query(query_str, con)
    df = df.set_index(["Target", "Alias", "Gene", "Transcript", "Barcode"])
    if not includeWT:
        df = df[df["lumc_category"] != "Wild-Type"]
        df = normalise_dataset_per_barcode(df)
    return df

def load_samples_01_to_06(raw = False, min_mut_events = MINIMUM_MUTATED_EVENTS, min_frequency = MINIMUM_FREQUENCY):
    query_str = 'select Barcode, Gene, Alias, outcome, \
        fraction_per_barcode, mutEvents, countEvents \
        from outcomes where Alias in ("MB01", "MB02", "MB03", "MB04", "MB05", "MB06") and \
        Gene is not "Empty" and Gene is not "Non-targeting"'

    if min_mut_events is not None:
        query_str = query_str + " and mutEvents >= {}".format(min_mut_events)

    if min_frequency is not None:
        query_str = query_str + " and fraction_per_barcode >= {}".format(min_frequency)

    con = sqlite3.connect(get_db_file())
    df = pd.read_sql_query(query_str, con)
    df_idx = df.set_index(["Alias", "Gene", "Barcode", "outcome"])

    print("loaded outcomes")   
    
    if raw:
        return df_idx

    df_idx_norm = normalise_dataset_per_barcode(df_idx)

    return df_idx_norm

def load_sample(sample_name, raw = False, min_mut_events = MINIMUM_MUTATED_EVENTS, min_frequency = MINIMUM_FREQUENCY, limit=None, set_index=True, con=None):
    query_str = 'select Barcode, Gene, Alias, outcome, \
        fraction_per_barcode, mutEvents, countEvents \
        from outcomes where Alias is "{}" and \
        Gene is not "Empty" and Gene is not "Non-targeting"'.format(sample_name)

    if min_mut_events is not None:
        query_str = query_str + " and mutEvents >= {}".format(min_mut_events)

    if min_frequency is not None:
        query_str = query_str + " and fraction_per_barcode >= {}".format(min_frequency)

    if limit is not None:
        query_str = query_str + " limit {}".format(limit)

    if con is None:
        con = sqlite3.connect(get_db_file())

    print("running query:")
    print(query_str)
    
    df = pd.read_sql_query(query_str, con)

    if set_index:
        df = df.set_index(["Alias", "Gene", "Barcode", "outcome"])

    print("loaded outcomes")   
    
    if raw:
        return df

    df_norm = normalise_dataset_per_barcode(df)

    return df_norm

def load_adamson_sample(sample_name, min_mut_events = MINIMUM_MUTATED_EVENTS, min_frequency = MINIMUM_FREQUENCY, limit=None, set_index=True):
    query_str = 'select Barcode, Alias, Gene, outcome, \
        fraction_per_barcode, mutEvents, countEvents, insSize \
        from outcomes where Alias is "{}"'.format(sample_name)
    
    if min_mut_events is not None:
        query_str = query_str + " and mutEvents >= {}".format(min_mut_events)

    if min_frequency is not None:
        query_str = query_str + " and fraction_per_barcode >= {}".format(min_frequency)

    if limit is not None:
        query_str = query_str + " limit {}".format(limit)

    con = sqlite3.connect(get_adamson_db_file())
    df = pd.read_sql_query(query_str, con)

    if set_index:
        df = df.set_index(["Alias", "Gene", "Barcode", "outcome"])

    print("loaded outcomes")   
    
    return df

def normalise_dataset_per_barcode(df_idx):
    if "Transcript" in df_idx.index.names:
        totals = df_idx[["fraction_per_barcode"]].groupby(["Target", "Alias", "Gene", "Transcript", "Barcode"]).sum()
    else:
        totals = df_idx[["fraction_per_barcode"]].groupby(["Target", "Alias", "Gene", "Barcode"]).sum()
    df_idx[["fraction_per_barcode"]] = df_idx[["fraction_per_barcode"]].div(totals)
    print("Normalised per Gene") 
    return df_idx


def normalise_dataset_per_gene(df_idx):
    df_idx_norm_mean = df_idx.groupby(level=["Alias", "Gene", "outcome"]).mean()
    totals = df_idx_norm_mean[["fraction_per_barcode"]].sum(level="Gene")
    df_idx_norm_mean = df_idx_norm_mean[["fraction_per_barcode"]].div(totals, level="Gene").droplevel(level="Alias")
    print("Normalised per Gene") 
    return df_idx_norm_mean


def remove_rare_outcomes(df_idx_norm_mean, min_frequency=200):
    num_genes = df_idx_norm_mean.index.get_level_values("Gene").unique().shape[0]
    frequent_outcomes = (df_idx_norm_mean.groupby("outcome").sum().div(num_genes)) > min_frequency
    frequent_outcomes = frequent_outcomes[frequent_outcomes["fraction_per_barcode"]]
    frequent_outcomes = frequent_outcomes.index

    df_idx_norm_mean_reduced = df_idx_norm_mean.loc[df_idx_norm_mean.index.get_level_values("outcome").isin(frequent_outcomes)]
    totals = df_idx_norm_mean_reduced[["fraction_per_barcode"]].sum(level="Gene")
    df_idx_norm_mean_reduced = df_idx_norm_mean_reduced[["fraction_per_barcode"]].div(totals, level="Gene")
    print("Normalised frequent outcomes > {}".format(min_frequency))
    return df_idx_norm_mean_reduced

def pivot_repair_outcome_profiles(df_idx_norm_mean_reduced):
    if df_idx_norm_mean_reduced["outcome"].unique().shape[0] > 100:
        raise Exception("There many be too many outcomes, pivot could take a long time")

    df_outcome_profiles = df_idx_norm_mean_reduced.reset_index().pivot(index=["Gene"], columns="outcome", values="fraction_per_barcode").fillna(0)
    df_outcome_profiles.index = df_outcome_profiles.index.str.lower()
    print("Repair outcome profiles created")
    return df_outcome_profiles


def get_target_site(sample_name):
    if sample_name in ["MB01", "MB02", "MB07", "MB08", "MBsub1", "MBsub2", "MBsub3"]:
        return "T1"
    if sample_name in ["MB03", "MB04"]:
        return "T2"
    if sample_name in ["MB05", "MB06"]:
        return "T3"

def get_target_site_TUD(sample_name):
    if "LV581" in sample_name:
        return "T1"
    elif "LV582" in sample_name:
        return "T2"
    elif "LV583" in sample_name:
        return "T3"

if __name__ == "__main__":
    df = load_lumc_outcomes(limit=1000000)
    print(df.shape)