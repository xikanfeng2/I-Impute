from sklearn.manifold import TSNE
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os
import sys

plt.switch_backend('agg')

methods = ['raw', 'gamma', 'scimpute', 'saver', 'magic', 'dca']
# ifile = sys.argv[1]
# gfile = sys.argv[2]
# ofile = sys.argv[3]

dataset = sys.argv[1]

tsnes = []

groups = [int(float(c)) for c in open('../experiment/dataset-1000/' +
                                      dataset+'.raw.csv').readline().strip().split(',')[1:]]

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
    data = pd.read_csv(infile, index_col=0, sep=sep)
    data = data.fillna(0)

    print('run tsne...')
    X_tsne = TSNE(n_components=2).fit_transform(data.T)
    tsnes.append(X_tsne)

print('drawing graph...')
plt.figure(figsize=(60, 10))
plt.title(dataset+' Dataset')
for index, tsne in enumerate(tsnes):
    plt_number = '16' + str(index+1)
    plt.subplot(int(plt_number))
    plt.scatter(tsne[:, 0], tsne[:, 1], c=groups, label=methods[index])
    plt.legend()
plt.savefig('../experiment/tsne-result/compared_result/' +
            dataset + '.png', dpi=300)
