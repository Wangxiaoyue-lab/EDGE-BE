# bayesian_analysis/analysis.py

import pandas as pd
import numpy as np
import pymc as pm
import pytensor.tensor as tt
from scipy.stats import norm

class BayesianAnalysis:
    def __init__(self, data, ratio_col, time_points, gene='gene', ef=0, fc='fc_column', count=0, ratio=0.6):
        self.data = data
        self.ratio_col = ratio_col.split(',')
        self.time_points = [float(item.strip()) for item in time_points.split(',')]
        self.gene = gene
        self.ef = ef
        self.fc = fc
        self.count = count
        self.ratio = ratio
        self.epsilon = 5
        
    def preprocess_data(self):
        init = 'count-0_ratio'
        self.data[self.ratio_col] = self.data[self.ratio_col].apply(pd.to_numeric, errors='coerce')
        self.data = self.data[self.data[self.ratio_col].apply(lambda x: (x <= self.count).sum() / len(x) < self.ratio, axis=1)]
        
        self.data['ef'] = pd.to_numeric(self.data['ef'], errors='coerce')
        self.data = self.data[(self.data[init] > 0) & (self.data['ef'].notnull()) & (self.data['ef'] >= self.ef)]
        
        self.sample_names = self.data['sgRNA']
        self.sgRNA_names = self.data[self.gene]
        self.N0 = self.data[init].values.astype(float)
        self.P = self.data[self.ratio_col].values.astype(float)/30000000  # scale
        self.p = self.data['ef'].values.astype(float)  # editing efficiency
        
        self.log2_ratio = np.log2((self.data[self.fc].values.astype(float)+self.epsilon) / (self.data[init].values.astype(float)+self.epsilon))
    
    def non_zero_gmean(self, arr, axis=None):
        mask = arr != 0
        sum_values = np.sum(np.where(mask, arr, 0), axis=axis)
        count_non_zero = np.sum(mask, axis=axis)
        valid_mask = count_non_zero != 0
        result = np.empty_like(sum_values, dtype=float)
        result[valid_mask] = sum_values[valid_mask] / np.sqrt(count_non_zero[valid_mask])
        result[~valid_mask] = 0
        return result
    
    def calculate_means(self):
        unique_sgRNAs = self.sgRNA_names.unique()
        num_sgRNAs = len(unique_sgRNAs)
        time_points = self.P.shape[1]
        
        self.P_mean = np.zeros((num_sgRNAs, time_points))
        self.N0_mean = np.zeros(num_sgRNAs)
        self.p_mean = np.zeros(num_sgRNAs)
        self.k_initial = np.zeros(num_sgRNAs)
        
        for i, sgRNA in enumerate(unique_sgRNAs):
            idx = self.sgRNA_names == sgRNA
            self.P_mean[i, :] = self.non_zero_gmean(self.P[idx, :], axis=0)
            self.N0_mean[i] = self.non_zero_gmean(self.N0[idx])
            self.p_mean[i] = self.non_zero_gmean(self.p[idx])
            self.k_initial[i] = self.non_zero_gmean(self.log2_ratio[idx])
            
        self.log_N0 = np.log(self.N0_mean + self.epsilon)
        
    def run_bayesian_model(self):
        t = np.array(self.time_points, dtype=np.float64)
        k = self.k_initial
        
        with pm.Model() as model:
            k_prior = pm.Normal('k', mu=k, sigma=0.1, shape=len(self.k_initial), initval=k)
            
            exp_t = tt.exp(t)
            exp_kt = tt.exp(k_prior[:, np.newaxis] * t[np.newaxis, :])
            
            N_i = (1 - self.p_mean[:, np.newaxis]) * self.N0_mean[:, np.newaxis] * exp_kt + self.p_mean[:, np.newaxis] * self.N0_mean[:, np.newaxis] * exp_kt
            N_total = tt.sum(N_i, axis=0)
            P_hat = N_i / N_total
            
            P_obs = pm.Normal('P_obs', mu=P_hat, sigma=0.01, observed=self.P_mean)
            
            advi_optimizer = pm.adam(learning_rate=0.01)
            advi_fit = pm.fit(n=20000, method='advi', obj_optimizer=advi_optimizer)
            
            self.trace = advi_fit.sample(draws=10000)
            
    def post_process_results(self):
        k_estimates = self.trace.posterior['k'].mean(dim=['chain', 'draw']).values
        k_samples = self.trace.posterior['k'].values
        
        k_mean = np.mean(k_estimates)
        k_std = np.std(k_estimates)
        
        k_zscores = (k_estimates - k_mean) / k_std
        p_values = 1 * (1 - norm.cdf(np.abs(k_zscores)))
        p_posterior = (k_samples > -0.05).mean(axis=(0, 1))
        
        output = pd.DataFrame({
            'sgRNA': self.sgRNA_names.unique(),
            'k_estimate': k_estimates,
            'k_zscore': k_zscores,
            'p_post': p_posterior,
            'p_value': p_values
        })
        output.fillna('.', inplace=True)
        return output
