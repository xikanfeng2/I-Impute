'''The goal of this demo is to show how to identify cell subpopulations based on latent
representations of gene expression learned by scScope.'''
import scscope as DeepImpute
import pandas as pd
import phenograph
import numpy as np
import sys

data = pd.read_csv(sys.argv[1], index_col=0)
latent_dim = 50
DI_model = DeepImpute.train(np.asarray(data.T), latent_dim, T=2, batch_size=64, max_epoch=300, num_gpus=4)
latent_code, imputed_val, _ = DeepImpute.predict(np.asarray(data.T), DI_model)
pd.DataFrame(imputed_val).T.to_csv(sys.argv[1]+'.scscope.csv')
