import pandas as pd
from musician.musician.util import encode_indel_properties, get_raw_file_for_sample
from tqdm import tqdm

def preprocess_raw_MUSIC_sample(sample_number):
    print("Loading raw data, takes about 30 seconds")
    raw_file = get_raw_file_for_sample(sample_number)
    df = pd.read_csv(raw_file, sep="\t", nrows=1000)

    print("Adding unique encoding per outcome")
    tqdm.pandas(desc="Progress")
    df["outcome"] = df.progress_apply(encode_indel_properties, axis=1)

    print("Selecting columns of interest")
    count_columns = ["Alias", "Gene", "Barcode", "outcome", "countEvents", "fraction"]
    count_df = df[count_columns]

    print("Saving processed file")
    file_to_save = raw_file.replace("/data/", "/data_processed/").replace(".txt", ".processed.txt")
    count_df.to_csv(file_to_save, index=False, sep="\t")

    print("Simple sanity check to see that each outcome and gene is unique")
    dups = count_df[["Gene", , "outcome"]].duplicated()
    if dups.any():
        print(count_df[dups].iloc[0,])
        raise Exception("Some gene/outcome pairs rows have been duplicated. Recheck function logic.")

    print("Finished", sample_number)

if __name__ == "__main__":
    preprocess_raw_MUSIC_sample("01")
    

