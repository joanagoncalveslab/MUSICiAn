import pandas as pd

from src.config import get_interim_dir

targetwisemutspectra = {}
mb01df = pd.read_pickle(get_interim_dir() + "MB01_gRNAwise_reduced_repair_outcome_profiles.pkl")
mb02df = pd.read_pickle(get_interim_dir() + "MB02_gRNAwise_reduced_repair_outcome_profiles.pkl")
targetwisemutspectra["T1"] = pd.concat([
    mb01df,
    mb02df,
])

mb03df = pd.read_pickle(get_interim_dir() + "MB03_gRNAwise_reduced_repair_outcome_profiles.pkl")
mb04df = pd.read_pickle(get_interim_dir() + "MB04_gRNAwise_reduced_repair_outcome_profiles.pkl")
targetwisemutspectra["T2"] = pd.concat([
    mb03df,
    mb04df,
])

mb05df = pd.read_pickle(get_interim_dir() + "MB05_gRNAwise_reduced_repair_outcome_profiles.pkl")
mb06df = pd.read_pickle(get_interim_dir() + "MB06_gRNAwise_reduced_repair_outcome_profiles.pkl")
targetwisemutspectra["T3"] = pd.concat([
    mb05df,
    mb06df,
])

for k in targetwisemutspectra:
    targetwisemutspectra[k].to_pickle(get_interim_dir() + "{}_gRNAwise_reduced_repair_outcome_profiles.pkl".format(k))
