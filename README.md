# I-Impute: a self-consistent method to impute single cell RNA sequencing data

I-Impute is a “self-consistent” method method to impute scRNA-seq data. I-Impute leverages continuous similarities and dropout probabilities and refines the data iteratively to make the final output "self-consistent". I-Impute exhibits robust imputation ability and follows the “self-consistency” principle. It offers perspicacity to uncover the underlying cell subtypes in real scRNA-Seq data. 

## Pre-requirements
* python3
* numpy>=1.16.1
* pandas>=0.23.4,<0.24
* scipy>=1.3.0
* scikit-learn>=0.21.1
* tasklogger>=0.4.0

### install requirements
```Bash
pip install -r requirements.txt
```

## Installation

### Installation with pip
To install with pip, run the following from a terminal:
```Bash
pip install i-impute
```

### Installation from Github
To clone the repository and install manually, run the following from a terminal:
```Bash
git clone https://github.com/xikanfeng2/I-Impute.git
cd I-Impute
python setup.py install
```

## Usage

### Quick start
The following code runs MAGIC on simulation data located in the I-Impute repository.

```Python
import iimpute
import pandas as pd

# read your reads count or RPKM or TPM data
data = pd.read_csv('simluation-data/sim-counts.csv', index_col=0)

# create I-Impute object
iimpute_operator = iimpute.IImpute(normalize=False)

# impute
imputed_data = iimpute_operator.impute(data)

# store result to a file
imputed_data.to_csv('your file name')

# iterative mode
iimpute_operator = iimpute.IImpute(normalize=False, iteration=True)

# impute
imputed_data = iimpute_operator.impute(data)

# store result to a file
imputed_data.to_csv('your file name')
```

### Parameters
```Python
IImpute(n=20, c_drop=0.5, p_pca=0.4, alpha=0.01, normalize=True, iteration=False, verbose=1)
```
Parameters

* n : int, optional, default: 20

    The nth of nearest neighbors on which to build kernel when calculating affinity matrix.

* c_drop : float, optional, default: 0.5

    Dropout event cutoff. For entry whose dropout probability is less than c_drop, we consider it as a real observation, its original value will remain. Otherwise, we conduct the imputation with the aid of information from similar cells.

* p_pca : float, optional, default: 0.4

    Percentage of variance explained by the selected components of PCA. It determines the nmumber of PCs used to calculate the distance between cells.

* alpha : float, optional, default: 0.01

    L1 penalty for Lasso regression.
    
* normalize : boolean, optional, default: True

    By default, I-Impute takes in an unnormalized matrix and performs library size normalization during the denoising step. However, if your data is already normalized or normalization is not desired, you can set normalize=False.

* iteration : boolean, optional, default: False

    The imputation process only performs once when False (it is equivalent to C-Impute described in our paper). The imputation process will iterate n times to achieve self-constistent imputation matrix.

* verbose : `int` or `boolean`, optional, default: 1

    If `True` or `> 0`, print status messages

## Cite us
Feng, X., Chen, L., Wang, Z., & Li, S. C. (2019). I-Impute: a self-consistent method to impute single cell RNA sequencing data. bioRxiv, 772723. doi: https://doi.org/10.1101/772723.

## Help
If you have any questions or require assistance using I-Impute, please contact us with xikanfeng2@gmail.com.