import sys
import pandas as pd
from tqdm import tqdm
from src.config import get_adamson_barcodes, get_interim_dir
from src.data.make_adamson_dataset import create_Adamson_outcomes_table
from src.data.load_dataset import load_adamson_lumc_outcomes
from src.data.filter import filter_barcodes
from src.data.categories import aggregate_outcomes
from src.data.profiles import generate_genewise_repair_outcome_profiles



if __name__ == "__main__":
    CMDS = ["all", "filter", "finalise"]
    LIMIT = None
    ALIASES = ["K562_SpCas9_target-1_HDR_oBA701_AX227_rep_{}".format(i) for i in range(1, 3)]
    TARGETS = ["Adamson_T{}".format(i) for i in range(1, 2)]

    test = (len(sys.argv) > 2) and (sys.argv[2] == "test")
    if test:
        LIMIT = 10000

    cmd = sys.argv[1] if len(sys.argv) > 1 else "filter"

    if cmd not in CMDS:
        raise Exception("{} is not a valid command, valid commands are {}".format(cmd, CMDS))

    if cmd == "all":
        # create reduced database table
        create_Adamson_outcomes_table()


    if test:
        dfs = []
        for sample in ALIASES:
            dfs.append(load_adamson_lumc_outcomes(limit=LIMIT, sample_name=sample))
        df = pd.concat(dfs)
    else:
        df = load_adamson_lumc_outcomes(limit=None)

    df = df.reset_index()
    df["Gene Symbol"] = df["Gene"]
    df["Gene"] = df["Gene"] + "-" + df["Transcript"]
    # convert to profiles
    p = df.pivot(index=["Target", "Alias", "Gene", "Barcode"], columns=["lumc_category", "del_len", "start"], values="countEvents").fillna(0)

    # aggregate categories
    idx = pd.IndexSlice
    Xs = []
    for alias in tqdm(ALIASES):
        X = p.loc[p.index.get_level_values("Alias") == alias, :]
        X = aggregate_outcomes(X, alias)
        Xs.append(X)
    profiles = pd.concat(Xs).fillna(0)

    if cmd in ["filter", "finalise", "all"]:
        # filter barcodes and save filtered barcodes file
        gene_results, _ = filter_barcodes(profiles)
        gene_results.reset_index().to_csv(get_adamson_barcodes(), sep="\t", index=False)

    if cmd in ["filter", "finalise", "aggregate", "all"]:
        # generated aggregated profiles per replicate
        samplewisemutspectra = generate_genewise_repair_outcome_profiles(data=profiles, level="Alias", source="Adamson")
        for k in samplewisemutspectra:
            samplewisemutspectra[k] = samplewisemutspectra[k].reset_index()
            tmp = samplewisemutspectra[k]["Gene"].str.split("-", n=1, expand=True)
            samplewisemutspectra[k]["Gene"] = tmp[0].str.capitalize()
            samplewisemutspectra[k]["Transcript"] = tmp[1]
            samplewisemutspectra[k] = samplewisemutspectra[k].set_index(["Target", "Alias", "Gene", "Transcript", "Barcode"])
            samplewisemutspectra[k].to_pickle(get_interim_dir() + "{}_gRNAwise_reduced_repair_outcome_profiles.pkl".format(k))
        # generate_genewise_repair_outcome_profiles(data=profiles, level="Target", source="Adamson")

        targetwisemutspectra = {}
        rep1 = samplewisemutspectra["K562_SpCas9_target-1_HDR_oBA701_AX227_rep_1"]
        rep2 = samplewisemutspectra["K562_SpCas9_target-1_HDR_oBA701_AX227_rep_2"]
        targetwisemutspectra["Adamson_T1"] = pd.concat([
            rep1,
            rep2,
        ])

        for k in targetwisemutspectra:
            print(k)
            targetwisemutspectra[k].to_pickle(get_interim_dir() + "{}_gRNAwise_reduced_repair_outcome_profiles.pkl".format(k))
    
    print("Done.")
