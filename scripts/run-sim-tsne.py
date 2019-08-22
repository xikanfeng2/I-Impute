from sklearn.manifold import TSNE
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os
import sys

plt.switch_backend('agg')

methods = {
        'Raw': '../experiment/simulation-data/sim-true-counts.csv',
        '88.45% Dropout': '../experiment/simulation-data/sim-counts.csv',
        'DCA': '../experiment/simulation-result/counts.dca.tsv',
        'SAVER': '../experiment/simulation-result/counts.saver.csv',
        'scImpute': '../experiment/simulation-result/counts.scimpute.csv',
        'C-Impute': '../experiment/simulation-result/counts.gamma.csv',
        'I-Impute': '../experiment/simulation-result/counts.gamma-saver.10.csv'
    }
# methods = ['true', 'raw', 'gamma']
# ifile = sys.argv[1]
# gfile = sys.argv[2]
# ofile = sys.argv[3]


tsnes = []

groups_str = pd.read_csv('../experiment/simulation-data/sim-cell-info.csv', index_col=0).values[:,2]
groups = [int(g.split('Group')[1]) for g in groups_str]

for method, infile in methods.items():
    sep = ','
    if method == 'DCA':
        sep = '\t'

    print('reading {0} data...'.format(method))
    print(infile)
    data = pd.read_csv(infile, index_col=0, sep=sep)
    data = data.fillna(0)

    print('run tsne...')
    X_tsne = TSNE(n_components=2).fit_transform(data.T)
    tsnes.append(X_tsne)

print('drawing graph...')
plt.figure(figsize=(40, 20))
plt.title('Simulation Dataset')
for index, tsne in enumerate(tsnes):
    # if index < 3:
    #     plt_number = '24' + str(index+1)
    # else:
    plt_number = '24' + str(index+1)
    plt.subplot(int(plt_number))
    plt.scatter(tsne[:, 0], tsne[:, 1], c=groups)
    plt.title(list(methods.keys())[index], fontsize=30)
    plt.tick_params(labelsize=20)
    # plt.legend()

metrics_data = pd.read_csv('../experiment/simulation-result/metrics.csv', index_col=0)
plt.subplot(248)
ind = np.arange(metrics_data.shape[1])
print(ind)
width = 0.12
for i in range(0, metrics_data.shape[0]):
    plt.bar(ind + width * i,
            metrics_data.values[i], width, label=metrics_data.index[i])
plt.xticks(ind + width * 3.5, ('ARI', 'NMI', 'SW'))
plt.tick_params(labelsize=20)
plt.legend(loc='best', prop={'size': 13})
plt.title('Evaluation Metrics', fontsize=30)
plt.savefig('../experiment/tsne-result/simulation-result/sim.png', dpi=300)
