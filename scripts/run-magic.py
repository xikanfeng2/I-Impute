import magic
import pandas as pd
import sys

data = pd.read_csv(sys.argv[1], index_col=0)
magic_op = magic.MAGIC()
emt_magic = magic_op.fit_transform(data.T, genes='all_genes')
emt_magic.T.to_csv(sys.argv[2])
