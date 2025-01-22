import pandas as pd
import numpy as np
from scipy import stats

def calculate_difference_to_controls(X, controls, samples, return_log2_fold_change=False):
    if X.columns.nlevels != 1 or X.index.nlevels != 2:
        error_str = """
        X should be the following format:
        Index should be two levels, first level is the replicate/targetsite. Second level is the gene.
        Columns should be a single level, with the mutation type
        """
        raise Exception(error_str)

    X = X.loc[samples]
    X = X.dropna()
    print(X.shape)

    Xpseudo = X.loc[X.index.get_level_values(level=1).isin(controls)]
    print(Xpseudo.shape)

    Xpseudo_gmean = pd.DataFrame(Xpseudo.groupby(level=0, axis=0).apply(stats.gmean).values.tolist(), index=X.index.unique(0), columns=X.columns)
    Xpseudo_gmean = Xpseudo_gmean.div(Xpseudo_gmean.sum(axis=1), axis=0)
    Xpseudo_gmean

    X = X.stack().reset_index().rename(columns={"level_0": "Sample"}).pivot(index="Gene", columns=["Sample", "lumc_category"], values=0).dropna()
    print(X.shape)

    Xpseudo_gmean = Xpseudo_gmean.stack()
    print(Xpseudo_gmean.shape)

    if return_log2_fold_change:
        Xdiff = np.log2(X/Xpseudo_gmean)
    else:
        Xdiff = X - Xpseudo_gmean
        print(Xdiff.shape)
    return Xdiff
