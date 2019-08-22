import pandas as pd
import numpy as np
import h5py
import sys
import os

input_file = sys.argv[1]
output_dir = sys.argv[2]
random_count = int(sys.argv[3])

data_mat = h5py.File(input_file)
x = np.array(data_mat['X'])
y = np.array(data_mat['Y'])

dataset = pd.DataFrame(x)
dataset.index = y
print(np.unique(dataset.index, return_counts=True))

sample_dataset = dataset.sample(random_count)
sample_dataset.T.to_csv(os.path.join(output_dir, input_file.split('.')[0] + '_1000.csv'))
print(np.unique(sample_dataset.index, return_counts=True))
