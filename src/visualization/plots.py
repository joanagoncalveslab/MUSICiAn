import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from scipy.cluster import hierarchy
from sklearn.metrics.pairwise import pairwise_distances

FANCONI_ANEMIA_GENE_SET = ['Atrip', 'Atr', 'Fancm', 'Faap24', 'Cenps', 'Cenpx', 'Telo2', 'Hes1'\
 , 'Faap100', 'Fanca', 'Fancb', 'Fancc', 'Fance', 'Fancf', 'Fancg', 'Fancl'\
    , 'Wdr48', 'Usp1', 'Ube2t', 'Fanci', 'Fancd2', 'Brca2', 'Palb2', 'Rad51c'\
        , 'Rad51', 'Brca1', 'Brip1', 'Fan1', 'Mlh1', 'Pms2', 'Rev1', 'Rev3l', 'Polh'\
            , 'Poli', 'Polk', 'Poln', 'Rmi1', 'Rmi2', 'Top3b', 'Top3a', 'Blm', 'Rpa1', 'Rpa2'\
                , 'Rpa3', 'Mus81', 'Eme1', 'Eme2', 'Ercc4', 'Ercc1', 'Slx1b', 'Slx4']

HR_GENE_SET = ['Ssbp1', 'Rad50', 'Mre11a', 'Nbn', 'Atm', 'Brca1', 'Bard1', 'Rbbp8', 'Brip1', 'Topbp1'\
      , 'Abraxas1', 'Uimc1', 'Babam1', 'Babam2', 'Brcc3', 'Brcc3dc', 'Palb2', 'Brca2', 'Sem1'\
        , 'Sycp3', 'Rpa1', 'Rpa2', 'Rpa3', 'Rad51', 'Rad52', 'Rad51b', 'Rad51c', 'Rad51d', 'Xrcc2', 'Xrcc3', 'Rad54l'\
      , 'Rad54b', 'Pold1', 'Pold2', 'Pold3', 'Pold4', 'Blm', 'Top3b', 'Top3a', 'Mus81', 'Eme1']

NHEJ_GENE_SET = ['Xrcc6', 'Xrcc5', 'Dclre1c', 'Prkdc', 'Poll', 'Polm', 'Dntt', 'Lig4', 'Xrcc4', 'Nhej1', 'Rad50', 'Mre11a', 'Fen1']

POLYMERASE_GENE_SET = ['Polr3a','Pold4','Poln','Polr1c','Pole','Polr3g','Polr3gl','Pole2','Polr2m','Polr2h','Polr2g','Polr1has',
               'Polr2e','Polr3k','Polr2i','Polr3f','Polr1h','Polr1f','Polr1b','Polr2c','Polr2b','Polg2','Pold2','Polb','Pola1',
               'Polr1g','Polr3d','Polq','Polr1a','Polr2j','Pole4','Polr2l','Pold1','Polm','Polr2k','Poll','Polr2a','Polk','Polr2f',
               'Polr1e','Pold3','Polr1d','Polr3b','Polr2d','Polh', 'Polg', 'Poldip3', 'Poldip2', 'Polr3c', 'Pole3','Polrmt','Poli','Pola2','Polr3e','Polr3h']

# plotting functions
def biplot(score, coeff, labels=None, components=None, plot_arrows=True, groups=None, group_order=None, palette="Set1"):
    xs = score[:,0]
    ys = score[:,1]
    n = coeff.shape[0]
    scalex = 1.0/(xs.max() - xs.min())
    scaley = 1.0/(ys.max() - ys.min())

    fig, ax = plt.subplots()

    if groups is None:
        ax.plot(xs * scalex,ys * scaley, alpha=1, marker='o', linestyle='')
    else:
        uniq_groups = groups.unique()
        
        if group_order is None:
            group_order = uniq_groups

        palette = sns.color_palette(palette, 12)
        colormapping = {group_order[i]: palette[i] for i in range(0, len(group_order))}
        colormapping

        for name in group_order:
            idx = (groups == name)
            ax.plot(xs[idx] * scalex, ys[idx] * scaley, alpha=0.5, marker='o', label=name, linestyle='', c=colormapping[name])
        ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))

    if plot_arrows:
        for i in range(n):
            ax.arrow(0, 0, coeff[i,0], coeff[i,1], color = 'black',alpha = 0.5, zorder=10)
            if labels is not None:
                ax.text(coeff[i,0]* 1.15, coeff[i,1] * 1.15, labels[i], color = 'black', ha = 'center', va = 'center')
    ax.set_xlim(-1,1)
    ax.set_ylim(-1,1)
    ax.set_xlabel(components[0])
    ax.set_ylabel(components[1])
    plt.grid()
    return scalex, scaley

def screeplot(explained_variance_ratio, ax, alias):
    PC_values = np.arange(len(explained_variance_ratio)) + 1
    ax.plot(PC_values, explained_variance_ratio/explained_variance_ratio.sum(), 'o-', linewidth=2, color='royalblue')
    ax.set_title('{} Scree Plot'.format(alias))
    ax.set_xlabel('Principal Component')
    ax.set_ylabel('Variance Explained')
    plt.ylim(bottom=0)


def loadingsplot(loadings, features, cluster=False, eigenvalues=None, figsize=None):
    num_pc = loadings.shape[1]
    pc_list = ["PC"+str(i) for i in list(range(1, num_pc+1))]
    loadings_df = pd.DataFrame(loadings, columns=pc_list)
    loadings_df['variable'] = features
    loadings_df = loadings_df.set_index('variable')
    z = None
    if not cluster:
        plt.figure(figsize=figsize)
        sns.heatmap(loadings_df.round(4), annot=loadings_df.shape[0] <= 15, cmap='vlag', annot_kws={"fontsize":7})
    if cluster:
        if eigenvalues is not None:
            dist = pairwise_distances((loadings * (np.max(eigenvalues) - eigenvalues)), metric = "euclidean")
            z = hierarchy.linkage(dist, method="ward")
            sns.clustermap(loadings_df.round(4), annot=loadings_df.shape[0] <= 15, \
                cmap='vlag', annot_kws={"fontsize":7}, col_cluster=False, row_linkage=z, figsize=figsize)
        else:
            sns.clustermap(loadings_df.round(4), annot=loadings_df.shape[0] <= 15, \
                cmap='vlag', annot_kws={"fontsize":7}, method="ward", col_cluster=False, figsize=figsize)
    return loadings_df, z


def psuedocontrolcomparisonplot(X_plot, pseudo_controls=[], compare_genes = [], background_genes = None, num_outcomes = 30, include_heatmap=True):
    columns = X_plot.columns
    X_plot = X_plot.div(X_plot.sum(axis=1), axis=0)
    X_plot = X_plot.groupby(["Gene"]).mean()
    X_plot = pd.DataFrame(X_plot.values.tolist(), index=X_plot.index, columns=columns)
    X_plot.columns = X_plot.columns.to_list()


    linepdata = X_plot.loc[background_genes if background_genes is not None else pseudo_controls, :].T
    linepdata.index.name = "Feature"
    outcome_order = list(linepdata.T.mean().sort_values(ascending=False).index)[:num_outcomes]

    if include_heatmap:
        fig, (ax, ax2) = plt.subplots(1, 2, sharey=True, figsize=(num_outcomes/5, num_outcomes/8), gridspec_kw={'width_ratios': [1.5, 1]})
    else:
        fig, ax = plt.subplots(figsize=(5.5, num_outcomes/8))

    for col in linepdata.columns:
        p = linepdata.loc[outcome_order[::-1], col].reset_index()
        y = p.index
        x = p[col]
        labels = p["Feature"]
        ax.plot(x, y, color="grey", alpha=0.1)

    if pseudo_controls is not None:
        psudeo_control_profile = X_plot.loc[pseudo_controls, :]
        psudeo_control_avg_profile = psudeo_control_profile.mean()
        logfolddf = np.log2(X_plot/psudeo_control_avg_profile)
        logfolddf.head()
        ax.plot(linepdata.loc[outcome_order[::-1],:].mean(axis=1), y, color="black", label="Pseudo-Control")

    for gene in compare_genes:
        p = X_plot.T.loc[outcome_order[::-1], gene].reset_index()
        y = p.index
        x = p[gene]
        labels = p["Feature"]
        ax.plot(x, y, label=gene)

    ax.legend()
    ax.set_yticks(y, labels)
    ax.set_xlabel("Frequency")
    ax.set_ylabel("Outcome")
    ax.grid()

    if include_heatmap:
        if pseudo_controls is None:
            raise Exception("If including heatmap, pseudo-controls must be provided.")
        limit = logfolddf.loc[compare_genes].abs().max().max()
        implotdata = logfolddf.loc[compare_genes, outcome_order].T.iloc[::-1]
        im = ax2.imshow(implotdata, vmax=limit, vmin=-limit, cmap='bwr', aspect='auto')
        plt.colorbar(im)
        ax2.set_xticklabels(ax2.get_xticks(), rotation = 45)
        ax2.set_xticks(np.arange(len(compare_genes)), compare_genes)

    return None

def gene_sample_frequency_change_heatmap(change_df, gene_to_vis, counts=None, figsize=(5, 5)):
    idx = pd.IndexSlice
    heatmap_df = change_df.loc[idx[:,gene_to_vis], np.sign(change_df.loc[idx[:,gene_to_vis], :].replace({-np.inf: np.nan})).sum(axis=0).abs().sort_values(ascending=False).index]
    # gene = "Rrm1"
    heatmap_df = heatmap_df.T.sort_index().T
    fig, ax2 = plt.subplots(1, 1, sharey=True, figsize=figsize)
    limit = heatmap_df.replace({-np.inf: np.nan}).abs().max().max()
    implotdata = heatmap_df.T
    im = ax2.imshow(implotdata, vmax=limit, vmin=-limit, cmap='bwr', aspect='auto')
    plt.colorbar(im)
    ax2.set_xticklabels(ax2.get_xticks(), rotation = 45)
    if counts is not None:
        ax2.set_xticks(np.arange(heatmap_df.shape[0]), ["{}\n({:.0f})".format(a, counts[a].xs(gene_to_vis, level="Gene").mean() if gene_to_vis in counts[a].index.unique("Gene") else 0) for a in heatmap_df.index.get_level_values("Sample")])
    else:
        ax2.set_xticks(np.arange(heatmap_df.shape[0]), [a for a in heatmap_df.index.get_level_values("Sample")]) 
    ax2.set_yticks(np.arange(heatmap_df.shape[1]), heatmap_df.columns.to_list())
    plt.title(gene_to_vis)


def mutational_spectra_clustermap(X_imputed, outliers, psuedo_controls, outlier_df=None, log_fold_change=False, vmin=None, vmax=None, method="ward", metric="correlation", figsize=(30, 12)):
    Xall = pd.concat(X_imputed, axis=0).dropna()
    print("All:", Xall.shape)

    Xout = Xall.loc[Xall.index.get_level_values(level=1).isin(outliers)]
    Xout.index.set_names("Sample", level=0, inplace=True)
    print("Outliers:", Xout.shape)

    Xpseudo = Xall.loc[Xall.index.get_level_values(level=1).isin(psuedo_controls)]
    print("Psuedo controls:", Xpseudo.shape)

    Xpseudo_gmean = pd.DataFrame(Xpseudo.groupby(level=0, axis=0).apply(stats.gmean).values.tolist(), index=list(X_imputed.keys()), columns=Xall.columns)
    Xpseudo_gmean = Xpseudo_gmean.div(Xpseudo_gmean.sum(axis=1), axis=0)
    Xpseudo_gmean

    Xout = Xout.stack().reset_index().pivot(index="Gene", columns=["Sample", "lumc_category"], values=0).dropna()
    print("Outliers, post pivot:", Xout.shape)

    Xpseudo_gmean = Xpseudo_gmean.stack()
    print("Mean psuedo control:", Xpseudo_gmean.shape)

    if log_fold_change:
        Xdiff = np.log2(Xout / Xpseudo_gmean)
    else:
        Xdiff = Xout - Xpseudo_gmean
    print("Log2 fold change:",Xdiff.shape)

    # Xdiff.columns = ['{}_{}'.format(i, j) for i, j in Xdiff.columns]
    Xdiff.index = Xdiff.index.to_list()

    if outlier_df is not None:
        allsampleoutliers = outlier_df.loc[Xdiff.index, :]
        nonrepairoutliers = allsampleoutliers[~allsampleoutliers[("Global", "isGORepair")]].index.to_series()
        boolrepairoutliers = Xdiff.index.to_series().isin(nonrepairoutliers)
        nonrepairoutliers

        pal = sns.color_palette('Dark2', 2)
        lut = dict(zip([True, False], pal))
        row_colors = pd.Series(boolrepairoutliers).map(lut)

    pal = sns.color_palette('Dark2', Xdiff.columns.get_level_values("lumc_category").unique().shape[0])
    lut = dict(zip(Xdiff.columns.get_level_values("lumc_category").unique(), pal))
    col_colors = Xdiff.columns.get_level_values("lumc_category").to_series().map(lut)

    ### TODO: Allow cluster genes by jensen shannon and cluster rows by correlation matrices
    # z_genes = hierarchy.linkage(Xdiff, metric="jensenshannon", method="average")
    # z_outcomes = hierarchy.linkage(Xdiff.T, metric="correlation", method="average")


    if vmin is not None or vmax is not None:
        cg = sns.clustermap(Xdiff.T, metric=metric, method=method, figsize=figsize, center=0, dendrogram_ratio=(.05, .1), cmap="RdBu", cbar_pos=(0, .9, .001, .01), \
            vmin=None, vmax=None)
    else:
        cg = sns.clustermap(Xdiff.T, metric=metric, method=method, figsize=figsize, center=0, dendrogram_ratio=(.05, .1), cmap="RdBu", cbar_pos=(0, .9, .001, .01), \
            robust=True)
    cg.ax_heatmap.yaxis.tick_left()
    cg.ax_heatmap.xaxis.tick_top()
    
    for i, tick_label in enumerate(cg.ax_heatmap.axes.get_xticklabels()):
        if outlier_df is not None:
            tick_text = tick_label.get_text()
            tick_label.set_color(row_colors[tick_text])
        tick_label.set_rotation(90)
    
    for i, tick_label in enumerate(cg.ax_heatmap.axes.get_yticklabels()):
        tick_text = tick_label.get_text().split("-")[1]
        tick_label.set_color(lut[tick_text])
        tick_label.set_rotation(0)
    plt.tight_layout()

    # return hierarchy.linkage(Xdiff.T, metric=metric, method=method)

def get_labels(Xdiff, labels_to_plot):
    labels = pd.Series(np.repeat("", Xdiff.shape[0]), index=Xdiff.index)
    inter_labels = np.intersect1d(Xdiff.index, labels_to_plot)
    # labels[inter_labels] = inter_labels
    return pd.Series(inter_labels)

def offdiag_f(x, y, **kwargs):
    ax = sns.scatterplot(x=x, y=y, palette=kwargs["palette"], hue=kwargs["hue"], alpha=kwargs["alpha"])

    if "labels" in kwargs:
        z = kwargs["labels"]
        del kwargs["labels"]
            
        if ("label" in kwargs) and (kwargs["label"] != "Other"):
            for i in range(len(x)):
                ax.annotate(str(z.values[i]), xy=(x.values[i], y.values[i]), fontsize=8,
                            textcoords="offset points",
                        color=kwargs.get("color","k"),
                        bbox=dict(pad=0.0, alpha=0.0, color='black'),
                        va='center', ha='center')

def pathway_pairplot(pp_df, pathway_name, pathway_gene_set, show_only_pathway_genes=False, 
                     hue=None, palette=None):
    if pathway_name is None:
        g = sns.PairGrid(data=pp_df)
        g.map_offdiag(offdiag_f, alpha=0.5, palette=palette, hue=hue)
        g.map_diag(sns.kdeplot, common_norm=False)
        plt.tight_layout()
        g.add_legend()
        return g

    pp_df[pathway_name] = pp_df.index.to_series()\
        .isin(pathway_gene_set).apply(lambda x: pathway_name if x else "Other")
    if show_only_pathway_genes:
        g =  sns.PairGrid(data=pp_df.loc[pp_df[pathway_name] == pathway_name]\
                          , hue=pathway_name\
                            , hue_order=[pathway_name])
    else:
        g = sns.PairGrid(data=pp_df, hue=pathway_name, hue_order=[pathway_name, "Other"])

    g.map_offdiag(offdiag_f, labels=get_labels(pp_df, pathway_gene_set), alpha=0.5, palette=palette, hue=hue)
    g.map_diag(sns.kdeplot, common_norm=False)
    g.fig.suptitle(pathway_name + " (# Genes: {})".format(pp_df.shape[0])) 
    plt.tight_layout()
    g.add_legend()
    return g
