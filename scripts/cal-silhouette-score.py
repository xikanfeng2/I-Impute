


import sys
import os
import pandas as pd
import numpy as np
from sklearn.metrics import silhouette_score

methods = ['raw', 'gamma', 'scimpute', 'saver', 'magic' ,'dca']
# ifile = sys.argv[1]
# gfile = sys.argv[2]
# ofile = sys.argv[3]

dataset = sys.argv[1]

ss = []

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
    
    print('run silhouette_score...')
    
    ss.append(silhouette_score(data.T, groups))

print('writing result...')
with open('../experiment/compare-evaluation-metrics/silhouette_score_'+dataset+'.csv', 'w') as output:
    output.write('Method,Silhouette Score\n')
    for index, value in enumerate(ss):
        output.write(methods[index]+','+str(value)+'\n')