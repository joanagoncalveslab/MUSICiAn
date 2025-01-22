import sqlite3
import pandas as pd
from tqdm import tqdm
from src.config import get_adamson_db_file, PROCESSED_DATA_PATH
from src.data.load_dataset import load_adamson_sample, get_target_site
from src.data.categories import convert_to_lumc_categories


# Reduce dataframe by categorising outcomes as LUMC categories

def create_Adamson_outcomes_table():
    conn = sqlite3.connect(get_adamson_db_file())
    replace = True

    for s in tqdm(range(1, 3)):
        sample = "K562_SpCas9_target-1_HDR_oBA701_AX227_rep_{}".format(s)
        df = load_adamson_sample(sample, set_index=False, limit=None)

        df["lumc_category"] = df.apply(lambda x: convert_to_lumc_categories(x["outcome"], x["insSize"]), axis=1)
        df["Target"] = "Adamson_T1"
        tmp = df["Gene"].str.split("-", n=1, expand=True)
        df["Gene"] = tmp[0].str.capitalize()
        df["Transcript"] = tmp[1]

        df_reduced = df[["Target", "Alias", "Gene", "Barcode", "Transcript", "lumc_category", "fraction_per_barcode", "countEvents"]].groupby(["Target", "Alias", "Gene", "Transcript", "Barcode", "lumc_category"]).sum().reset_index()
        df_reduced[["lumc_category", "start", "del_len", "ins_len"]] = pd.DataFrame(df_reduced["lumc_category"].to_list(), index=df_reduced.index, columns=["lumc_category", "start", "del_len", "ins_len"])
        df_reduced = df_reduced.groupby(["Target", "Alias", "Gene", "Transcript", "Barcode", "lumc_category", "del_len", "start"]).sum().reset_index()
        df_reduced.to_sql(name="lumc_outcomes", con=conn, if_exists="replace" if replace else "append", index=False)

        print("Reduced rows in {} from {} to {} ({}%)".format(sample, df.shape[0], df_reduced.shape[0], df_reduced.shape[0]/df.shape[0]))

        replace = False

    cursor = conn.cursor()
    cursor.execute("create index lumc_outcomes_idx on lumc_outcomes (Target, Alias, Gene, Transcript, Barcode, lumc_category, del_len, start, fraction_per_barcode, countEvents)")
    conn.close()


if __name__ == "__main__":
    create_Adamson_outcomes_table()
