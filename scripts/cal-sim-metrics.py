from sklearn.metrics import adjusted_rand_score
from sklearn.metrics import jaccard_score
from sklearn.metrics import normalized_mutual_info_score
from sklearn.cluster import AgglomerativeClustering
from sklearn.metrics import silhouette_score
import pandas as pd
import numpy as np


methods = ['true', 'raw', 'gamma', 'scimpute', 'saver', 'magic', 'dca']

methods = {
    'Raw': '../experiment/simulation-data/sim-true-counts.csv',
    '88.45% Dropout': '../experiment/simulation-data/sim-counts.csv',
    'DCA': '../experiment/simulation-result/counts.dca.tsv',
    'SAVER': '../experiment/simulation-result/counts.saver.csv',
    'scImpute': '../experiment/simulation-result/counts.scimpute.csv',
    'scImpute+': '../experiment/simulation-result/counts.gamma.csv',
    'I-Impute': '../experiment/simulation-result/counts.gamma-saver.1.csv'
}

groups_str = pd.read_csv(
    '../experiment/simulation-data/sim-cell-info.csv', index_col=0).values[:, 2]
groups = np.asarray([int(g.split('Group')[1]) for g in groups_str])
groups = groups - 1
print(groups)
with open('../experiment/simulation-result/metrics.csv', 'w') as output:
    output.write(
        'methods,ARI,NMI,SW\n')
    for method, infile in methods.items():
        sep = ','
        if method == 'DCA':
            sep = '\t'

        print('reading {0} data...'.format(method))
        print(infile)
        data = pd.read_csv(infile, index_col=0, sep=sep)
        data = data.fillna(0)
        
        clustering = AgglomerativeClustering(n_clusters=3).fit(data.T)
        labels = clustering.labels_
        print(labels)
        output.write(','.join([method, str(adjusted_rand_score(groups, labels)),
                               str(normalized_mutual_info_score(groups, labels, average_method='min')), 
                               str(silhouette_score(data.T, groups))])+'\n')
