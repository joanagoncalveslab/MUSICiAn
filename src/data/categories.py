import numpy as np

def split_outcome_to_parts(outcome):
    type_, delRelativeStart, delRelativeEnd, inserted_seq, mh_len = outcome.split("|")
    mh_len = int(mh_len[:-2])
    return type_, delRelativeStart, delRelativeEnd, inserted_seq, mh_len

def convert_to_lumc_categories(outcome, ins_len = None):
    type_, delRelativeStart, delRelativeEnd, inserted_seq, mh_len = split_outcome_to_parts(outcome)
    del_len = int(delRelativeEnd) - int(delRelativeStart)
    ins_len = len(inserted_seq) if ins_len is None else ins_len
    category = ""

    if type_ == "WT":
        category = "Wild-Type"
        return (category, 0, 0, 0)

    if type_ == "DELINS":
        category = "Deletion with insertion"
        return (category, 0, 0, 0)

    if type_ == "TINS":
        category = "Deletion with templated insertion"
        return (category, 0, 0, 0)

    if type_ == "INSERTION":
        if ins_len == 1:
            category = "1bp insertion - {}".format('?' if inserted_seq == '' else inserted_seq)
            return (category, 0, 0, ins_len)
        elif ins_len == 2:
            category = "2bp insertion - {}".format('??' if inserted_seq == '' else inserted_seq)
            return (category, 0, 0, ins_len)
        else:
            category = "3+bp insertion"
            return (category, 0, 0, 3)

    if type_ == "TANDEMDUPLICATION":
        category = "Tandem Duplication"
        return (category, 0, 0, 0)

    if type_ == "TANDEMDUPLICATION_COMPOUND":
        category = "Tandem Duplication plus"
        return (category, 0, 0, 0)

    if type_ == "HDR":
        category = "Homology Directed Repair"
        return (category, 0, 0, del_len)

    if type_ == "DELETION":
        if mh_len <5:
            category =  "Deletion {}bp microhomology".format(mh_len)
        elif mh_len <= 15:
            category =  "Deletion 5-15bp microhomology"
        elif mh_len > 15:
            category =  "Deletion >15bp microhomology"

        return (category, delRelativeStart if del_len < 30 else 30, del_len if del_len < 30 else 30, ins_len)

    raise Exception("No category was found for", type_)



def aggregate_outcomes(X, alias):
    counts = X.sum(axis=1)
    X = X.groupby("lumc_category", axis=1).sum()
    X = X.loc[:, (X != 0).any(axis=0)]

    insertioncols = list(X.columns.values[X.columns.get_level_values(0).to_series().str.contains("1bp insertion")]) +\
        list(X.columns.values[X.columns.get_level_values(0).to_series().str.contains("2bp insertion")]) + ["3+bp insertion"]
    X["Any Insertion"] = X.loc[:, insertioncols].sum(axis=1)
    # insertionothercols = list(X.columns.values[X.columns.get_level_values(0).to_series().str.contains("2bp insertion")]) + ["3+bp insertion"]
    # X["Other Insertion"] = X.loc[:, insertionothercols].sum(axis=1)

    largemhcols = ["Deletion 3bp microhomology", "Deletion 4bp microhomology", "Deletion 5-15bp microhomology"]
    X["Deletion 3+bp microhomology"] = X.loc[:, np.intersect1d(X.columns, largemhcols)].sum(axis=1)

    tandemduplications = X.columns.values[X.columns.get_level_values(0).to_series().str.contains("Tandem")]
    X["Deletion with insertion"] = X[["Deletion with templated insertion"]].sum(axis=1) + X.loc[:, tandemduplications].sum(axis=1) + X[["Deletion with insertion"]].sum(axis=1)

    X = X.drop(insertioncols, axis=1)
    X = X.drop(tandemduplications, axis=1)
    X = X.drop("Deletion with templated insertion", axis=1)
    X = X.drop(np.intersect1d(X.columns, largemhcols), axis=1)

    if not X.sum(axis=1).eq(counts).all():
        raise Exception("Reads are missing, recheck categorisation")
    return X