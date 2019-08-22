from sklearn.metrics import adjusted_rand_score
from sklearn.metrics import jaccard_score
from sklearn.metrics import normalized_mutual_info_score
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
import pandas as pd
import numpy as np
import sys


methods = ['raw', 'gamma', 'scimpute', 'saver', 'magic', 'dca']
dataset = sys.argv[1]


groups = [int(float(c)) for c in open('../experiment/dataset-1000/' +
                                      dataset+'.raw.csv').readline().strip().split(',')[1:]]
with open('../experiment/dataset-1000-result/'+dataset+'.metrics.csv', 'w') as output:
    output.write(
        'methods,adjusted_rand_score,jaccard_similarity_score,normalized_mutual_info_score,silhouette_score\n')
    for method in methods:
        sep = ','
        if method == 'raw':
            infile = '../experiment/dataset-1000/' + dataset + '.' + method + '.csv'
        elif method == 'dca':
            infile = '../experiment/dataset-1000-result/' + \
                method + '/' + dataset + '.' + method + '.tsv'
            sep = '\t'
        else:
            infile = '../experiment/dataset-1000-result/' + \
                method + '/' + dataset + '.' + method + '.csv'

        print('reading {0} data...'.format(method))
        print(infile)
        data = pd.read_csv(infile, index_col=0, sep=sep)
        data = data.fillna(0)

        clustering = KMeans(n_clusters=5).fit(data.T)
        labels = clustering.labels_
        print(labels.shape)
        output.write(','.join([method, str(adjusted_rand_score(groups, labels)), str(jaccard_score(groups, labels, average='macro')),
                               str(normalized_mutual_info_score(groups, labels, average_method='min')), str(silhouette_score(data.T, groups))])+'\n')
