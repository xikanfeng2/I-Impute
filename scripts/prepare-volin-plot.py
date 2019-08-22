import pandas as pd
import numpy as np
import seaborn as sns

methods = {
    'Raw': '../experiment/simulation-data/sim-true-counts.csv',
    '88.45% Dropout': '../experiment/simulation-data/sim-counts.csv',
    'DCA': '../experiment/simulation-result/counts.dca.tsv',
    'SAVER': '../experiment/simulation-result/counts.saver.csv',
    'scImpute': '../experiment/simulation-result/counts.scimpute.csv',
    'C-Impute': '../experiment/simulation-result/counts.gamma.csv',
    'I-Impute': '../experiment/simulation-result/counts.gamma-saver.10.csv'
}

with open('../experiment/simulation-result/volin.csv', 'w') as output:
    output.write('Method,value\n')
    for method, infile in methods.items():
        sep = ','
        if method == 'DCA':
            sep = '\t'

        print('reading {0} data...'.format(method))
        print(infile)
        data = pd.read_csv(infile, index_col=0, sep=sep)
        data = data.fillna(0)

        for row in data.values:
            for i in row:
                output.write(method+','+str(i)+'\n')
