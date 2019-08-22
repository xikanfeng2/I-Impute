import pandas as pd
import numpy as np
import sys

data = pd.read_csv(sys.argv[1], sep='\t', index_col=0)
clusters = pd.read_csv(sys.argv[2], header=None, index_col=0, sep=' ')

sample = data.loc[:,clusters.index]
sample.to_csv(sys.argv[3])