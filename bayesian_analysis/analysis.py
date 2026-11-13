"""Class-based interface for the LFCBayes model.

This module keeps the public workflow of the original GitHub implementation::

    analysis.preprocess_data()
    analysis.calculate_means()
    analysis.run_bayesian_model()
    result = analysis.post_process_results()

The statistical implementation follows ``lfcbayes_note2_final.py``.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pymc as pm
import pytensor.tensor as tt


DEFAULT_CONTROL_PATTERN = r"control|AAVS1|Non-target"


class BayesianAnalysis:
    """Fit the updated LFCBayes model to a merged count-ratio table.

    Parameters through ``ratio`` retain the order and defaults of the original
    GitHub class.  The remaining keyword arguments expose options added by the
    updated command-line implementation.
    """

    def __init__(
        self,
        data,
        ratio_col,
        time_points,
        gene="gene",
        ef=0,
        fc="fc_column",
        count=0,
        ratio=0.6,
        *,
        nt_pattern="Non-target",
        control_pattern=DEFAULT_CONTROL_PATTERN,
        scale_factor=30_000_000.0,
        epsilon=5.0,
        prior_sigma=0.1,
        obs_sigma=None,
        min_obs_sigma=1e-8,
        advi_steps=20_000,
        draws=10_000,
        learning_rate=0.01,
    ):
        self.data = data.copy()
        self.ratio_col = self._parse_csv(ratio_col)
        self.time_points = np.asarray(
            [float(item) for item in self._parse_csv(time_points)], dtype=np.float64
        )
        self.gene = gene
        self.ef = float(ef)
        self.fc = fc
        self.count = float(count)
        self.ratio = float(ratio)
        self.nt_pattern = nt_pattern
        self.control_pattern = control_pattern
        self.scale_factor = float(scale_factor)
        self.epsilon = float(epsilon)
        self.prior_sigma = float(prior_sigma)
        self.obs_sigma = None if obs_sigma is None else float(obs_sigma)
        self.min_obs_sigma = float(min_obs_sigma)
        self.advi_steps = int(advi_steps)
        self.draws = int(draws)
        self.learning_rate = float(learning_rate)

        if len(self.ratio_col) != len(self.time_points):
            raise ValueError("The number of ratio columns must match the number of time points.")

    @staticmethod
    def _parse_csv(value):
        if isinstance(value, str):
            return [item.strip() for item in value.split(",") if item.strip()]
        return list(value)

    @staticmethod
    def _infer_init_column(dataframe):
        for candidate in ("count-0_ratio", "count-A0_ratio", "count-C0_ratio"):
            if candidate in dataframe.columns:
                return candidate
        raise ValueError("Could not infer the T0 ratio column from the input table.")

    @staticmethod
    def _safe_std(values, fallback):
        values = np.asarray(values, dtype=float)
        values = values[np.isfinite(values)]
        if values.size <= 1:
            return fallback
        std = np.std(values, ddof=1)
        return fallback if not np.isfinite(std) or std <= 0 else float(std)

    @staticmethod
    def _empirical_one_sided_pvalue(beta_value, null_values, direction):
        null_values = np.asarray(null_values, dtype=float)
        null_values = null_values[np.isfinite(null_values)]
        if null_values.size == 0:
            return np.nan
        if direction == "negative":
            extreme = np.sum(null_values <= beta_value)
        elif direction == "positive":
            extreme = np.sum(null_values >= beta_value)
        else:
            raise ValueError("direction must be 'negative' or 'positive'.")
        return (extreme + 1.0) / (null_values.size + 1.0)

    def preprocess_data(self):
        """Validate, coerce, and filter the input rows."""
        self.init_column = self._infer_init_column(self.data)

        required = set(self.ratio_col) | {self.init_column, "ef", self.gene, self.fc}
        missing = sorted(required.difference(self.data.columns))
        if missing:
            raise ValueError(f"Missing required columns: {', '.join(missing)}")
        if self.fc not in self.ratio_col:
            raise ValueError("fc (the prior column) must be one of ratio_col.")

        numeric_columns = list(dict.fromkeys(self.ratio_col + [self.init_column, "ef", self.fc]))
        self.data[numeric_columns] = self.data[numeric_columns].apply(
            pd.to_numeric, errors="coerce"
        )
        keep = self.data[self.ratio_col].apply(
            lambda row: (row <= self.count).sum() / len(row) < self.ratio, axis=1
        )
        self.data = self.data.loc[keep].copy()
        self.data = self.data.loc[
            (self.data[self.init_column] > 0)
            & self.data["ef"].notna()
            & (self.data["ef"] >= self.ef)
        ].copy()

        if self.data.empty:
            raise ValueError("No rows remain after filtering.")
        return self

    def calculate_means(self):
        """Aggregate barcodes by the grouping column and construct the prior.

        The method name is retained for API compatibility.  The updated model
        uses robust medians for abundance and editing-efficiency aggregation.
        """
        if not hasattr(self, "init_column"):
            raise RuntimeError("Call preprocess_data() before calculate_means().")

        prior_source = self.data[self.fc].to_numpy(dtype=float)
        init_source = self.data[self.init_column].to_numpy(dtype=float)
        prior_lfc = pd.Series(
            np.log2((prior_source + self.epsilon) / (init_source + self.epsilon)),
            index=self.data.index,
            dtype=float,
        )

        output_names = []
        p0_median = []
        p_edit_median = []
        p_median = []
        grouped_ratios = []
        prior_mu = []

        for group_name, group in self.data.groupby(self.gene, sort=False):
            ratios = group[self.ratio_col].to_numpy(dtype=float) / self.scale_factor
            init_prop = group[self.init_column].to_numpy(dtype=float) / self.scale_factor
            edit_eff = group["ef"].to_numpy(dtype=float)

            output_names.append(group_name)
            p0_median.append(np.nanmedian(init_prop))
            p_edit_median.append(np.nanmedian(edit_eff))
            p_median.append(np.nanmedian(ratios, axis=0))
            grouped_ratios.append(ratios)
            # Mean log2 fold change equals the log2 geometric mean fold change.
            values = prior_lfc.loc[group.index].to_numpy(dtype=float)
            values = values[np.isfinite(values)]
            prior_mu.append(float(values.mean()) if values.size else 0.0)

        self.output_names = np.asarray(output_names, dtype=object)
        names = pd.Series(self.output_names).astype(str)
        self.is_nt = names.str.contains(
            self.nt_pattern, case=False, regex=True, na=False
        ).to_numpy(dtype=bool)
        self.is_control = names.str.contains(
            self.control_pattern, case=False, regex=True, na=False
        ).to_numpy(dtype=bool)
        if not np.any(self.is_nt):
            raise ValueError(f"No groups matched nt_pattern={self.nt_pattern!r}.")
        if not np.any(self.is_control):
            raise ValueError(f"No groups matched control_pattern={self.control_pattern!r}.")

        self.p0_median = np.asarray(p0_median, dtype=float)
        self.p_edit_median = np.asarray(p_edit_median, dtype=float)
        self.p_median = np.asarray(p_median, dtype=float)
        self.grouped_ratios = grouped_ratios
        raw_prior_mu = np.asarray(prior_mu, dtype=float)
        self.nt_geomean_center = float(np.nanmedian(raw_prior_mu[self.is_nt]))
        self.prior_mu = raw_prior_mu - self.nt_geomean_center
        self.beta_init = self.prior_mu.copy()
        self.prior_time = float(self.time_points[self.ratio_col.index(self.fc)])

        if not (
            np.all(np.isfinite(self.p0_median))
            and np.all(np.isfinite(self.p_edit_median))
            and np.all(np.isfinite(self.p_median))
            and np.all(np.isfinite(self.prior_mu))
        ):
            raise ValueError("Non-finite group summaries remain after aggregation.")

        if self.obs_sigma is None:
            residuals = [
                np.abs(group_ratio - group_median[np.newaxis, :]).ravel()
                for group_ratio, group_median in zip(self.grouped_ratios, self.p_median)
            ]
            residuals = np.concatenate(residuals) if residuals else np.array([], dtype=float)
            residuals = residuals[np.isfinite(residuals)]
            sigma = np.median(residuals) if residuals.size else self.min_obs_sigma
            self.obs_sigma_ = float(
                self.min_obs_sigma if not np.isfinite(sigma) or sigma < self.min_obs_sigma else sigma
            )
        else:
            self.obs_sigma_ = max(self.obs_sigma, self.min_obs_sigma)
        return self

    def run_bayesian_model(self):
        """Fit the updated editing-efficiency-aware model with ADVI."""
        if not hasattr(self, "prior_mu"):
            raise RuntimeError("Call calculate_means() before run_bayesian_model().")

        with pm.Model() as self.model:
            beta = pm.Normal(
                "beta",
                mu=self.prior_mu,
                sigma=np.full(self.prior_mu.shape, self.prior_sigma, dtype=float),
                shape=len(self.output_names),
                initval=self.beta_init,
            )
            growth_term = (
                1.0 - self.p_edit_median[:, np.newaxis]
            ) + self.p_edit_median[:, np.newaxis] * tt.exp(
                beta[:, np.newaxis] * self.time_points[np.newaxis, :]
            )
            numerator = self.p0_median[:, np.newaxis] * growth_term
            p_hat = numerator / tt.sum(numerator, axis=0, keepdims=True)
            pm.Normal("P_obs", mu=p_hat, sigma=self.obs_sigma_, observed=self.p_median)

            optimizer = pm.adam(learning_rate=self.learning_rate)
            self.advi_fit = pm.fit(
                n=self.advi_steps, method="advi", obj_optimizer=optimizer
            )
            self.trace = self.advi_fit.sample(draws=self.draws)
        return self.trace

    def post_process_results(self):
        """Return posterior estimates calibrated to all groups and NT controls."""
        if not hasattr(self, "trace"):
            raise RuntimeError("Call run_bayesian_model() before post_process_results().")

        beta_estimates = self.trace.posterior["beta"].mean(dim=["chain", "draw"]).values
        beta_samples = self.trace.posterior["beta"].values
        nt_beta = beta_estimates[self.is_nt]

        all_mean = float(np.mean(beta_estimates))
        all_sd = self._safe_std(beta_estimates, self.prior_sigma)
        nt_mean = float(np.mean(nt_beta))
        nt_sd = self._safe_std(nt_beta, self.prior_sigma)
        z_all = (beta_estimates - all_mean) / all_sd
        z_nt = (beta_estimates - nt_mean) / nt_sd

        p_all_upper = np.asarray([
            self._empirical_one_sided_pvalue(value, beta_estimates, "positive")
            for value in beta_estimates
        ])
        p_all_lower = np.asarray([
            self._empirical_one_sided_pvalue(value, beta_estimates, "negative")
            for value in beta_estimates
        ])
        p_nt_upper = np.asarray([
            self._empirical_one_sided_pvalue(value, nt_beta, "positive")
            for value in beta_estimates
        ])
        p_nt_lower = np.asarray([
            self._empirical_one_sided_pvalue(value, nt_beta, "negative")
            for value in beta_estimates
        ])
        post_up = (beta_samples > 0.0).mean(axis=(0, 1))
        post_down = (beta_samples < 0.0).mean(axis=(0, 1))

        output = pd.DataFrame({
            "sgRNA": self.output_names,
            "beta": beta_estimates,
            "post_prob_up": post_up,
            "post_prob_down": post_down,
            "z_all": z_all,
            "p_all_upper": p_all_upper,
            "p_all_lower": p_all_lower,
            "z_nt": z_nt,
            "p_nt_upper": p_nt_upper,
            "p_nt_lower": p_nt_lower,
            "control": self.is_control,
            "nt": self.is_nt,
            "p0": self.p0_median,
            "edit_eff": self.p_edit_median,
            "obs_sigma": self.obs_sigma_,
            "prior_mu": self.prior_mu,
            "prior_sigma": self.prior_sigma,
            "prior_col": self.fc,
            "prior_mode": "geomean_nt_centered",
            "nt_center": self.nt_geomean_center,
            "prior_time": self.prior_time,
            "time_points": ",".join(
                str(int(x)) if float(x).is_integer() else str(x) for x in self.time_points
            ),
        })

        # Backward-compatible aliases for consumers of the original class.
        output["k_estimate"] = output["beta"]
        output["k_zscore"] = output["z_all"]
        output["p_post"] = output["post_prob_up"]
        output["p_value"] = output[["p_all_upper", "p_all_lower"]].min(axis=1)
        return output

    def fit(self):
        """Run the full workflow and return the result table."""
        self.preprocess_data()
        self.calculate_means()
        self.run_bayesian_model()
        return self.post_process_results()
