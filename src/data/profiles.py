import pandas as pd
from scipy import stats
from tqdm import tqdm
from src.data.load_dataset import load_lumc_outcomes
from src.config import get_common_barcodes, get_adamson_barcodes, get_validation_barcodes


def load_data_from_database(limit=None, sample_name=None):
    # takes about 6 minutes to load the data into memory
    df = load_lumc_outcomes(limit=limit, sample_name=sample_name)

    X = df.reset_index()

    # combine repair outcome category and deletion length
    # X["Feature"] = X["lumc_category"] + ", " + X["del_len"].astype(str)

    # pivot to create profiles
    X = X.pivot(index=["Target", "Alias", "Gene", "Barcode"], columns=["lumc_category", "del_len", "start"], values="countEvents").fillna(0)

    return X


def generate_genewise_repair_outcome_profiles(data, level="Alias", source="LUMC", filter_barcodes=True, validation=False, downsample=None):
    idx = pd.IndexSlice

    # targets = ["T1", "T2", "T3"]
    # samples = ["MB01", "MB02", "MB03", "MB04", "MB05", "MB06"]

    iterable = data.index.unique(level=level).to_list()
    profiles = {}

    for i in tqdm(iterable):
        X = data.loc[data.index.get_level_values(level) == i, :]
    
        # read barcodes, and choose the valid barcodes to keep 
        if source == "LUMC":
            filtered_barcodes = pd.read_csv(get_common_barcodes() if not validation else get_validation_barcodes(downsample=downsample), sep="\t")
        elif source == "Adamson":
            filtered_barcodes = pd.read_csv(get_adamson_barcodes(), sep="\t")
        else:
            filtered_barcodes = pd.read_csv(get_common_barcodes() if not validation else get_validation_barcodes(downsample=downsample), sep="\t")

        if filter_barcodes:
            filtered_barcodes = filtered_barcodes[filtered_barcodes["Filtered"].isin([0,5])][["Target", "Alias", "Gene", "Barcode"]]
        filtered_barcodes = filtered_barcodes[filtered_barcodes[level] == i]
        kept_indexes = pd.MultiIndex.from_frame(filtered_barcodes)
        # remove bad gRNAs
        profiles[i] = X.loc[kept_indexes] 
    return profiles