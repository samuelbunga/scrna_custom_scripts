import os
import sys
import gseapy
import numpy as np
import pandas as pd

# Setting global variables
gmt_database = "/Users/sbunga/PycharmProjects/INDRA/scrna/test-project/annotation_script/gmt_database/"


def set_wd(outdir):
    """Set working directory to the provided path"""
    try:
        os.mkdir(outdir)
        os.chdir(outdir)
        os.mkdir('rnk')
        os.mkdir('label_data')
        print("Working directory set to: "+os.getcwd())

    except FileExistsError:
        pass


def to_df(seurat_data_raw, delim):
    """
    :param seurat_data_raw:
    Convert the raw seurat tab/comma separated file
    into a pandas dataframe
    :return:
    Pandas dataframe
    """
    df = pd.read_csv(seurat_data_raw, sep=delim, header=0)
    return df


def to_csv(genes, pvals, cluster):
    """
    :param genes:
    List of genes
    :param pvals:
    List of pvalues
    :param cluster:
    Cluster number
    :return:
    Create a rank file and save it to ./rnk directory
    """
    data = {"Gene": genes,
            "p_vals": pvals}
    rank_df = pd.DataFrame(data)
    rank_df.to_csv("rnk/cluster_" + str(cluster) + ".rnk", sep="\t",
                   index=False, header=False)


def make_rank(seurat_df):
    """ This function accepts a seurat cluster gene
    dataframe and creates a rank file for each cluster """
    full_cluster = {0}
    genes, pvals = [], []
    count = 0
    for i in range(0, len(seurat_df)):
        cluster = df.iloc[i]['cluster']
        if cluster not in full_cluster:
            full_cluster.add(cluster)
            to_csv(genes, pvals, count)
            count = cluster
            # reset the genes and pvals list to empty
            genes, pvals = [], []
        genes.append(seurat_df.iloc[i]['gene'])
        pvals.append(seurat_df.iloc[i]['avg_logFC'])

    # writing the data from the last loop
    to_csv(genes, pvals, count)


def annotate_clusters(nclusters, species):
    """
    :param nclusters: total number of clusters
    :param species: Species name
    :return: creates label data for each cluster and saves them into
    ./label_data directory
    Also, returns a list of top cell types for every cluster
    """
    # Annotate each cluster
    top = []
    for each_cluster in range(0, nclusters):
        try:
            gseapy.prerank("rnk/cluster_" + str(each_cluster) + ".rnk",
                           gmt_database + species + ".gmt",
                           'label_data/' + str(each_cluster) + '_folder')
        except:
            pass
        if os.path.isfile('label_data/' + str(each_cluster) + '_folder/gseapy.prerank.gene_sets.report.csv'):
            this_top = pd.read_csv(
                'label_data/' + str(each_cluster) + '_folder/' + 'gseapy.prerank.gene_sets.report.csv', header=0)
            label = this_top['Term'].iloc[0]
            top.append(str(label))
        else:
            top.append("unknown" + str(each_cluster))
    return top


if __name__ == "__main__":
    set_wd("/Users/sbunga/PycharmProjects/INDRA/scrna/"
           "test-project/annotation_script/test-outputs")
    df = to_df("/Users/sbunga/PycharmProjects/INDRA/scrna/"
           "test-project/annotation_script/seurat_clusters.txt", delim="\t")
    make_rank(df)
    top =  annotate_clusters(20, "mouse")
