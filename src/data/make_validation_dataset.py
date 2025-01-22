import sqlite3
import math
import pandas as pd
import numpy as np
from tqdm import tqdm
from src.config import get_validation_db_file, PROCESSED_DATA_PATH
from src.data.load_dataset import load_sample, get_target_site, get_target_site_TUD
from src.data.categories import convert_to_lumc_categories


# Reduce dataframe by categorising outcomes as LUMC categories
experiment = "TUD"

def create_LUMC_outcomes_table(downsample=None, limit=None):
    conn = sqlite3.connect(get_validation_db_file(experiment=experiment))
    replace = True
    table_name = "lumc_outcomes"

    # if (downsample is not None) and downsample not in ["half", "quarter", "hundreth"]:
    #     raise Exception("Invalid option")
    # elif downsample is None:
    print(table_name)
    # else:
    #     table_name = table_name + "_" + downsample
    #     print(table_name)

    if experiment == "TUD":
        all_samples = ['105973-001-001_LV581_rep1', '105973-001-002_LV582_rep1', '105973-001-003_LV583_rep1', \
            '105973-001-004_LV581_rep2', '105973-001-005_LV582_rep2', '105973-001-006_LV583_rep2', \
            '105973-001-007_LV581_rep3', '105973-001-008_LV582_rep3', '105973-001-009_LV583_rep3']
    else:
        all_samples = ["MBsub1", "MBsub2", "MBsub3"] 


    for sample in all_samples:
        print("Working on", sample)
        df = load_sample(sample, raw = True, set_index=False, limit=limit, con=conn)

        if downsample is not None:
            downsampled_dfs = []
            df = df.set_index("Barcode")
            barcodes = df.index.unique()

            for bc in tqdm(barcodes):
                barcode_df = df.loc[bc]
                barcode_df = barcode_df.loc[~barcode_df["outcome"].str.contains("WT")]
                bc_gene = barcode_df["Gene"].iloc[0]
                bc_alias = barcode_df["Alias"].iloc[0]

                bc_total_count = barcode_df["countEvents"].sum()
                bc_outcomes = barcode_df["outcome"]
                bc_probs = barcode_df["countEvents"]/bc_total_count
                if downsample == "half":
                    n = math.floor(bc_total_count * 0.5)
                elif downsample == "hundreth":
                    n = math.floor(bc_total_count * 0.01)
                else:
                    n = math.floor(bc_total_count * 0.25)
                draw = np.random.choice(bc_outcomes, n, p=bc_probs)
                tmp = pd.Series(draw).value_counts()
                tmp.name = "count"
                tmp = tmp.to_frame().reset_index().rename(columns={"index":"outcome", "count": "countEvents"})
                tmp["Barcode"] = bc
                tmp["Gene"] = bc_gene
                tmp["Alias"] = bc_alias
                tmp["fraction_per_barcode"] = tmp["countEvents"]/n
                tmp["mutEvents"] = n
                downsampled_dfs.append(tmp)
                tmp=None
            
            df = pd.concat(downsampled_dfs, axis=0)

        df["lumc_category"] = df["outcome"].apply(convert_to_lumc_categories)
        if experiment is None:
            df["Target"] = df["Alias"].apply(get_target_site)
        elif experiment == "TUD":
            df["Target"] = df["Alias"].apply(get_target_site_TUD)
        

        df_reduced = df[["Target", "Alias", "Gene", "Barcode", "lumc_category", "fraction_per_barcode", "countEvents"]].groupby(["Target", "Alias", "Gene", "Barcode", "lumc_category"]).sum().reset_index()
        df_reduced[["lumc_category", "start", "del_len", "ins_len"]] = pd.DataFrame(df_reduced["lumc_category"].to_list(), index=df_reduced.index, columns=["lumc_category", "start", "del_len", "ins_len"])
        df_reduced = df_reduced.groupby(["Target", "Alias", "Gene", "Barcode", "lumc_category", "del_len", "start"]).sum().reset_index()
        df_reduced.to_sql(name=table_name, con=conn, if_exists="replace" if replace else "append", index=False)

        print("Reduced rows in {} from {} to {} ({}%)".format(sample, df.shape[0], df_reduced.shape[0], df_reduced.shape[0]/df.shape[0]))

        replace = False

    cursor = conn.cursor()
    cursor.execute("create index {0}_idx on {0} \
                   (Target, Alias, Gene, Barcode, lumc_category, del_len, start, fraction_per_barcode, countEvents)".format(table_name))
    conn.close()


if __name__ == "__main__":
    create_LUMC_outcomes_table(downsample=None, limit=None)
