from sklearn.metrics import adjusted_rand_score
from sklearn.metrics import jaccard_score
from sklearn.metrics import normalized_mutual_info_score
from sklearn.cluster import AgglomerativeClustering
from sklearn.metrics import silhouette_score
import pandas as pd
import numpy as np
import sys



dataset = sys.argv[1]

methods = {
    'Raw': '../experiment/'+dataset+'/counts.raw.csv',
    'DCA': '../experiment/'+dataset+'/counts.dca.tsv',
    'SAVER': '../experiment/'+dataset+'/counts.saver.csv',
    'scImpute': '../experiment/'+dataset+'/counts.scimpute.csv',
    'scImpute+': '../experiment/'+dataset+'/counts.gamma.csv',
    'I-Impute': '../experiment/'+dataset+'/counts.gamma-saver.1.csv'
}

data = pd.read_csv(methods['Raw'], index_col=0)

groups = [i.split('.')[0] for i in data.columns]
groups = np.unique(groups, return_inverse=True)[1]

print(groups)
with open('../experiment/'+dataset+'/metrics.csv', 'w') as output:
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
        # if method == 'scImpute+' or method == 'I-Impute':
        # data = data/data.sum().replace(0, 1)

        clustering = AgglomerativeClustering(n_clusters=int(sys.argv[2])).fit(data.T)
        labels = clustering.labels_
        print(labels)
        output.write(','.join([method, str(adjusted_rand_score(groups, labels)),
                               str(normalized_mutual_info_score(
                                   groups, labels, average_method='min')),
                               str(silhouette_score(data.T, groups))])+'\n')
