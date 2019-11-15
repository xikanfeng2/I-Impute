from sklearn.decomposition import PCA
from sklearn.metrics.pairwise import euclidean_distances
from sklearn.linear_model import Lasso
from scipy.stats import gamma, norm
from scipy.special import digamma
from scipy.optimize import brentq
import tasklogger
import pandas as pd
import numpy as np
import warnings
warnings.filterwarnings("ignore")

from . import utils


class IImpute:

    def __init__(self, n=20, c_drop=0.5, p_pca=0.4, alpha=0.01, normalize=True, iteration=False, verbose=1):
        self.ZERO_VALUE = np.log10(1.01)
        self.n = n
        self.c_drop = c_drop
        self.p_pca = p_pca
        self.alpha = alpha
        self.normalize = normalize
        self.iteration = iteration
        self._check_params()
        self.verbose = verbose
        if self.normalize:
            self.ZERO_VALUE = 0.01
        tasklogger.set_level(verbose)

    def _check_params(self):
        """Check I-Impute parameters

        This allows us to fail early - otherwise certain unacceptable
        parameter choices, such as n='10.5', would only fail after
        minutes of runtime.

        Raises
        ------
        ValueError : unacceptable choice of parameters
        """
        utils.check_positive(n=self.n, c_drop=self.c_drop, p_pca=self.p_pca, alpha=self.p_pca)
        utils.check_int(n=self.n)
        utils.check_between(v_min=0, v_max=1, c_drop=self.c_drop, p_pca=self.p_pca)
        utils.check_bool(normalize=self.normalize)
        utils.check_bool(iteration=self.iteration)

    def set_params(self, **params):
        """Set the parameters of I-Impute.

        Any parameters not given as named arguments will be left at their
        current value.

        Parameters
        ----------

        n : int, optional, default: 20
            nth of nearest neighbors on which to build kernel when calculating affinity matrix.

        c_drop : float, optional, default: 0.5
            Dropout event cutoff. For entry whose dropout probability is less than c_drop, we consider 
            it as a real observation, its original value will remain. Otherwise, we conduct the imputation 
            with the aid of information from similar cells.

        p_pca : float, optional, default: 0.4
            Percentage of variance explained by the selected components of PCA. It determines the nmumber of PCs 
            used to calculate the distance between cells.

        alpha : float, optional, default: 0.01
            L1 penalty for Lasso regression.
        
        normalize : boolean, optional, default: True
            By default, I-Impute takes in an unnormalized matrix and performs library size normalization during 
            the denoising step. However, if your data is already normalized or normalization is not desired, you 
            can set normalize=False.

        iteration : boolean, optional, default: False
            The imputation process only performs once when False (it is equivalent to C-Impute described in our paper).
            The imputation process will iterate n times to achieve self-constistent imputation matrix.

        verbose : `int` or `boolean`, optional, default: 1
            If `True` or `> 0`, print status messages

        Returns
        -------
        self
        """

        # kernel parameters
        if 'n' in params and params['n'] != self.n:
            self.n = params['n']
            del params['n']
        if 'c_drop' in params and params['c_drop'] != self.c_drop:
            self.c_drop = params['c_drop']
            del params['c_drop']
        if 'p_pca' in params and params['p_pca'] != self.p_pca:
            self.p_pca = params['p_pca']
            del params['p_pca']
        if 'alpha' in params and params['alpha'] != self.alpha:
            self.alpha = params['alpha']
            del params['alpha']
        if 'normalize' in params and params['normalize'] != self.normalize:
            self.normalize = params['normalize']
            del params['normalize']
        if 'iteration' in params and params['iteration'] != self.iteration:
            self.iteration = params['iteration']
            del params['iteration']
        if 'verbose' in params:
            self.verbose = params['verbose']
            tasklogger.set_level(self.verbose)
            del params['verbose']

        self._check_params()
        return self

    def _cal_explained_component_number(self, explained_variance_ratio_):
        """Calculates index for selected PCs

        Parameters
        ----------
        explained_variance_ratio_ : array, shape (n_components,)
            The amount of variance explained by each of the selected components.

        Returns
        -------
        index: int
            The index of selected components. For example, index = 2 means the first 
            three components of PCA will be used to calculate distance between cells.
        """
        sum = 0.0
        for index, r in enumerate(explained_variance_ratio_):
            sum += r
            if sum >= self.p_pca:
                return index + 1
        return 0


    def _detect_outliers(self, min_dist_vector):
        """Detect outlier cells

        Parameters
        ----------
        min_dist_vector : array, shape (n cells)
            The set of distances to the nearest beighbor of each cells            

        Returns
        -------
        index array: array
            The index array for all outlier cells
        """
        upper_quartile = np.percentile(min_dist_vector, 75)
        lower_quartile = np.percentile(min_dist_vector, 25)
        cut_off = upper_quartile + 1.5 * (upper_quartile - lower_quartile)
        return np.where(min_dist_vector > cut_off)


    def _root_function(self, alpha, target):
        """Calculate alpha for gamma distribution

        Parameters
        ----------
        alpha : float
            The target alpha value for gamma distribution

        target : float
            The target value     

        Returns
        -------
        alpha: float
            The target alpha value for gamma distribution
        """
        return np.log(alpha) - digamma(alpha) - target


    def _update_gamma_param(self, row, mix_p_gamma):
        """Update alpha, beta params for gamma distribution

        Parameters
        ----------
        row : array, shape (n genes)
            The reads counts array for a gene

        mix_p_gamma : object
            The mix gamma distribution

        Returns
        -------
        alpha, beta: float
            The alpha, beta value for gamma distribution
        """
        tp_s = np.sum(mix_p_gamma)
        tp_t = np.sum(mix_p_gamma * row)
        tp_u = np.sum(mix_p_gamma * np.log(row))
        tp_v = -tp_u / tp_s - np.log(tp_s / tp_t)
        if tp_v <= 0:
            alpha = 20
        else:
            alpha0 = (3 - tp_v + np.sqrt(np.power(tp_v - 3, 2) + 24 * tp_v)) / 12 / tp_v
            if alpha0 >= 20:
                alpha = 20
            else:
                alpha = brentq(self._root_function, 0.9 * alpha0, 1.1 *
                            alpha0, disp=False, args=(tp_v))

        ## need to solve log(x) - digamma(x) = tp_v
        ## We use this approximation to compute the initial value
        beta = tp_s / tp_t * alpha
        return alpha, beta


    def _EM_algorithm(self, remained_data):
        """EM algorithm for calculating the droup out matrix

        Parameters
        ----------
        remained_data : matrix, shape (m x n)
            The reads count matrix which passed the outlier detection

        Returns
        -------
        D: matrix, shape (m x n)
            The dropout probability matrix D
        """
        # params = np.zeros((remained_data.shape[0], 5))
        D = np.zeros(remained_data.shape)

        # remained_data = M X N, row is genes and columns is remained cells
        for index, row in enumerate(remained_data):
            row[np.isnan(row)] = self.ZERO_VALUE
            # init params (pi, a, b, mean, std)
            pi = np.sum(row == self.ZERO_VALUE) / row.shape[0]  # pi
            if pi == 0:
                pi = 0.01
            shape = 0.5  # shape
            rate = 1  # rate
            row_non_zreo = row[row > self.ZERO_VALUE]
            mean = np.mean(row_non_zreo)  # mean
            std = np.std(row_non_zreo)  # std
            if std == 0 or np.isnan(std):
                std = 0.01
            eps = 10
            iteration = 0
            loglik_old = 0

            # EM alorithm
            while(eps > 0.5):
                # calculate weight
                p_gamma = pi * gamma.pdf(row, a=shape, scale=1/rate)
                p_norm = (1 - pi) * norm.pdf(row, loc=mean, scale=std)
                p_norm[np.isnan(p_norm)] = self.ZERO_VALUE
                mix_p_gamma = p_gamma/(p_gamma + p_norm)
                mix_p_norm = 1 - mix_p_gamma

                # update params
                pi = np.sum(mix_p_gamma) / mix_p_gamma.shape[0]
                mean = np.sum(mix_p_norm * row) / np.sum(mix_p_norm)
                std = np.sqrt(
                    np.sum(mix_p_norm * np.power(row - mean, 2)) / np.sum(mix_p_norm))
                if std == 0 or np.isnan(std):
                    std = 0.01
                if np.sum(mix_p_gamma) == 0:
                    break
                shape, rate = self._update_gamma_param(row, mix_p_gamma)
                dmix = pi * gamma.pdf(row, a=shape, scale=1/rate) + \
                    (1 - pi) * norm.pdf(row, loc=mean, scale=std)
                loglik = np.sum(np.log10(dmix))
                eps = np.power(loglik - loglik_old, 2)
                loglik_old = loglik
                iteration += 1
                if iteration > 100:
                    break

            # params[index] = [pi, shape, rate, mean, std]
            p_gamma = pi * gamma.pdf(row, a=shape, scale=1/rate)
            p_norm = (1 - pi) * norm.pdf(row, loc=mean, scale=std)
            mix_p_gamma = p_gamma/(p_gamma + p_norm)
            D[index] = mix_p_gamma
        return D

    def _cimpute(self, data):
        """Main function of C-Impute

        Parameters
        ----------
        data : matrix, shape (m x n)
            The raw reads count matrix

        Returns
        -------
        imputed_data: matrix, shape (m x n)
            The imputed matrix
        """
        # tasklogger.log_info('reading data...')
        # read data
        raw_data = data.fillna(0)

        if self.normalize:
            tasklogger.log_info('normalizing data by library size...')
            # normalize by library size
            norm_data = (raw_data*np.power(10, 6)/raw_data.sum().replace(0, 1)).values
        else:
            norm_data = raw_data.values

        tasklogger.log_info('preprocessing data...')
        # remove zero sum genes and cells
        filtered_rows_indexes = np.where(
            np.all(norm_data == 0, axis=1))
        filtered_rows_data = np.delete(
            norm_data, filtered_rows_indexes[0], axis=0)
        filtered_columns_indexes = np.where(
            np.all(filtered_rows_data == 0, axis=0))
        filtered_data = np.delete(
            filtered_rows_data, filtered_columns_indexes[0], axis=1)

        # log(x + 1.01)
        if self.normalize:
            log_data = np.log10(filtered_data + 1.01)
        else:
            log_data = filtered_data + self.ZERO_VALUE

        tasklogger.log_info('performing pca...')
        # pca
        pca = PCA()
        pca_data = pca.fit_transform(log_data.T)
        selected_pca_data = pca_data[:, :self._cal_explained_component_number(
            pca.explained_variance_ratio_)].T

        tasklogger.log_info('detecting outlier cells...')
        # remove outlier cells
        # 1. cal distance matrix
        dist_matrix = euclidean_distances(selected_pca_data.T, selected_pca_data.T)
        dist_matrix[dist_matrix == 0.0] = np.inf
        # 2. get min distance vector for each cell
        min_dist_vector = dist_matrix.min(axis=0)
        # 3. find outlier cells
        outlier_indexes = self._detect_outliers(min_dist_vector)
        tmp_remained_data = np.delete(log_data, outlier_indexes[0], axis=1)
        remained_pca_data = np.delete(
            selected_pca_data, outlier_indexes[0], axis=1)

        # remove rows in which all values is zero
        all_zeros_rows_indexes = np.where(
            np.all(tmp_remained_data == self.ZERO_VALUE, axis=1))
        remained_data = np.delete(
            tmp_remained_data, all_zeros_rows_indexes[0], axis=0)

        tasklogger.log_info('calculating the affinity matrix...')
        # cal affinity matrix
        remained_dist_matrix = euclidean_distances(
            remained_pca_data.T, remained_pca_data.T)
        remained_dist_matrix[remained_dist_matrix == 0.0] = np.inf
        if self.n >= int(remained_dist_matrix.shape[0]/2):
            self.n = int(remained_dist_matrix.shape[0]/2)
        nth_smallest_dist_index = self.n - 1
        nth_smallest_dist_vector = np.partition(
            remained_dist_matrix, nth_smallest_dist_index, axis=1)[:, nth_smallest_dist_index]
        
        exp_matrix = np.exp(-remained_dist_matrix /
                            (2 * np.power(nth_smallest_dist_vector, 2)))
        larged_indexes = np.column_stack(
            np.where(remained_dist_matrix <= nth_smallest_dist_vector))
        remained_dist_matrix[remained_dist_matrix >
                             nth_smallest_dist_vector] = 0            
        for index in larged_indexes:
            remained_dist_matrix[index[0]][index[1]
                                        ] = exp_matrix[index[0]][index[1]]
        tasklogger.log_info('calculating the droupout probability matrix using EM algorithm...')
        # EM algorithm
        D = self._EM_algorithm(remained_data)
        
        tasklogger.log_info('calculating the weight matrix using non-negative least squares lasso regression...')
        imputed_data = np.zeros(remained_data.shape)
        for index, row in enumerate(remained_dist_matrix):
            tasklogger.log_info('imputing gene expression for cell ' + str(index))
            D_j = D[:, index]
            X_j = remained_data[:, index]
            Y = (1 - D_j) * X_j
            non_self_dist = np.delete(row, index)  # N-1 X 1
            non_self_data = np.delete(remained_data, index, axis=1)  # M X N-1
            non_self_D = np.delete(D, index, axis=1)

            diag = np.diag(non_self_dist)  # N-1 X N-1
            X = np.matmul(diag, non_self_data.T).T * (1 - non_self_D)  # M X N-1

            lasso = Lasso(alpha=self.alpha, positive=True, max_iter=3000)
            lasso.fit(X, Y)
            imputed_data[:, index] = lasso.predict(X)

        tasklogger.log_info('processing imputed data...')
        imputed_indexes = np.column_stack(np.where(D >= self.c_drop))

        for index in imputed_indexes:
            remained_data[index[0]][index[1]] = imputed_data[index[0]][index[1]]

        tasklogger.log_info('recovering unimputed data...')
        # recover deleted rows and columns
        # 1. insert rows
        for index in all_zeros_rows_indexes[0]:
            remained_data = np.insert(
                remained_data, index, tmp_remained_data[index], axis=0)
        # 2. insert columns
        for index in outlier_indexes[0]:
            remained_data = np.insert(
                remained_data, index, log_data[:, index], axis=1)
        # 3. recover original values for each value
        remained_data[remained_data < self.ZERO_VALUE] = self.ZERO_VALUE
        if self.normalize:
            remained_data = np.power(10, remained_data) - 1.01
        else:
            remained_data = remained_data - self.ZERO_VALUE
        # 4. insert columns
        for index in filtered_columns_indexes[0]:
            remained_data = np.insert(
                remained_data, index, filtered_rows_data[:, index], axis=1)
        # 5. insert rows
        for index in filtered_rows_indexes[0]:
            remained_data = np.insert(
                remained_data, index, norm_data[index], axis=0)
        # 6. recover original values
        if self.normalize:
            remained_data = remained_data * raw_data.sum().replace(0, 1).values / \
                np.power(10, 6)

        tasklogger.log_info('generating the final imputed matrix...')
        # recover indexes name and columns name
        final_imputed_data = pd.DataFrame(remained_data)
        final_imputed_data.index = raw_data.index
        final_imputed_data.columns = raw_data.columns

        # set percision
        final_imputed_data = round(final_imputed_data, 3)

        return final_imputed_data

    def impute(self, data):
        """Main function of I-Impute

        Parameters
        ----------
        data : matrix, shape (m x n)
            The raw reads count matrix

        Returns
        -------
        imputed_data: matrix, shape (m x n)
            The imputed matrix, pandas Dataframe object
        """
        tasklogger.log_start('I-Impute')
        imputed_data = None
        if self.iteration:
            exp_mse = 0.1
            mse = 100
            previous_imputed_data = data
            iteration = 1
            while mse > exp_mse:
                tasklogger.log_info('iteratively impute for the {0}th time'.format(iteration))
                current_imputed_data = self._cimpute(previous_imputed_data)
                dist_matrix = (current_imputed_data - previous_imputed_data)**2
                n_values = data.shape[0] * data.shape[1]
                mse = np.sqrt(dist_matrix.values.sum()/n_values)
                previous_imputed_data = current_imputed_data
                iteration += 1
            
            imputed_data = previous_imputed_data
        else:
            imputed_data = self._cimpute(data)
        
        tasklogger.log_complete('I-Impute')
        return imputed_data
        

    
