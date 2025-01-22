import sys
import pandas as pd

from tqdm import tqdm

from src.config import get_common_barcodes, get_interim_dir, get_validation_barcodes, FILTER_COUNT
from src.data.make_dataset import create_LUMC_outcomes_table
from src.data.make_validation_dataset import create_LUMC_outcomes_table as create_validation_LUMC_outcomes_table
from src.data.load_dataset import load_lumc_outcomes
from src.data.filter import filter_barcodes
from src.data.categories import aggregate_outcomes
from src.data.profiles import generate_genewise_repair_outcome_profiles

if __name__ == "__main__":
    CMDS = ["all", "filter", "finalise"]
    VALIDATION = True 
    DOWNSAMPLE = "hundreth" if VALIDATION else None #[None, "half", "quarter", "hundreth"]
    LIMIT = None
    ALIASES = ["MB0{}".format(i) for i in range(1, 7)] if not VALIDATION else ["MBsub{}".format(i) for i in range(1, 4)]
    TARGETS = ["T{}".format(i) for i in range(1, 4)] if not VALIDATION else ["T1"]
    FILTER_BARCODES = True

    test = (len(sys.argv) > 2) and (sys.argv[2] == "test")
    if test:
        LIMIT = 1000

    cmd = sys.argv[1] if len(sys.argv) > 1 else "filter"

    if cmd not in CMDS:
        raise Exception("{} is not a valid command, valid commands are {}".format(cmd, CMDS))

    if cmd == "all":
        # create reduced database table
        if not VALIDATION:
            create_LUMC_outcomes_table()
        else:
            create_validation_LUMC_outcomes_table(downsample=DOWNSAMPLE)


    # takes about 6 minutes to load the data into memory
    p = {}
    print("loading data into memory")
    for alias in tqdm(ALIASES):
        df = load_lumc_outcomes(limit=LIMIT, sample_name=alias, validation=VALIDATION, downsample=None if not VALIDATION else DOWNSAMPLE)
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

    if cmd in ["filter", "all"]:
        # filter barcodes and save filtered barcodes file
        gene_results, _ = filter_barcodes(profiles)
        gene_results.reset_index().to_csv(get_common_barcodes() if not VALIDATION else get_validation_barcodes(downsample=DOWNSAMPLE), sep="\t", index=False)

    if cmd in ["filter", "finalise", "all"]:
        # generated aggregated profiles per replicate
        samplewisemutspectra = generate_genewise_repair_outcome_profiles(data=profiles, level="Alias", source="LUMC", \
                                                                         filter_barcodes=FILTER_BARCODES, validation=VALIDATION, downsample=DOWNSAMPLE)
        for k in samplewisemutspectra:
            if not VALIDATION:
                samplewisemutspectra[k].to_pickle(get_interim_dir() + "{}_gRNAwise_{}_repair_outcome_profiles.{}.pkl"\
                                                .format(k, "reduced" if FILTER_BARCODES else "all", FILTER_COUNT))
            else:
                samplewisemutspectra[k].to_pickle(get_interim_dir() + "{}_gRNAwise_{}_repair_outcome_profiles.{}.{}.pkl"\
                                                .format(k, "reduced" if FILTER_BARCODES else "all", FILTER_COUNT, \
                                                          "full" if DOWNSAMPLE is None else DOWNSAMPLE))

        targetwisemutspectra = {}
        if not VALIDATION:
            mb01df = samplewisemutspectra["MB01"]
            mb02df = samplewisemutspectra["MB02"]
            targetwisemutspectra["T1"] = pd.concat([
                mb01df,
                mb02df,
            ])
            mb03df = samplewisemutspectra["MB03"]
            mb04df = samplewisemutspectra["MB04"]
            targetwisemutspectra["T2"] = pd.concat([
                mb03df,
                mb04df,
            ])

            mb05df = samplewisemutspectra["MB05"]
            mb06df = samplewisemutspectra["MB06"]
            targetwisemutspectra["T3"] = pd.concat([
                mb05df,
                mb06df,
            ])
        else:
            mbsub1df = samplewisemutspectra["MBsub1"]
            mbsub2df = samplewisemutspectra["MBsub2"]
            mbsub3df = samplewisemutspectra["MBsub3"]
            targetwisemutspectra["T1"] = pd.concat([
                mbsub1df,
                mbsub2df,
                mbsub3df
            ])

        if not VALIDATION:
            for k in targetwisemutspectra:
                targetwisemutspectra[k].to_pickle(get_interim_dir() + "{}_gRNAwise_{}_repair_outcome_profiles.{}.pkl".format(k, "reduced" if FILTER_BARCODES else "all", FILTER_COUNT))
        else:
            for k in targetwisemutspectra:
                targetwisemutspectra[k].to_pickle(get_interim_dir() + "{}_validation_gRNAwise_{}_repair_outcome_profiles.{}.{}.pkl"\
                                                  .format(k, "reduced" if FILTER_BARCODES else "all", FILTER_COUNT, \
                                                          "full" if DOWNSAMPLE is None else DOWNSAMPLE))
    
    print("Done.")
