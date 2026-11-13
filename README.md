# EDGE-BE

EDGE-BE is an editing-activity-adjusted Bayesian dynamic model for
multi-time-point CRISPR base-editing screens.

## Updated model

The current implementation:

- separates unedited and edited cell fractions in the growth model;
- aggregates barcode measurements with robust medians;
- builds a fold-change prior and centers it on Non-target controls;
- estimates observation noise from within-group barcode residuals;
- reports posterior up/down probabilities and empirical upper/lower-tail
  probabilities relative to both all groups and Non-target controls.

## Installation

Python 3.9 or newer is recommended.

```bash
pip install git+https://github.com/Wangxiaoyue-lab/EDGE-BE.git
```

For a local clone:

```bash
git clone https://github.com/Wangxiaoyue-lab/EDGE-BE.git
cd EDGE-BE
pip install .
```

## Command-line usage

```bash
bayesian-analysis \
  -f input.tsv \
  -o result.tsv \
  -n "count-D6_ratio,count-D14_ratio,count-D18_ratio" \
  -t "6,14,18" \
  -g gene \
  -e 0 \
  --prior-col count-D14_ratio \
  -a 0 \
  -r 0.6
```

Required arguments:

- `-f, --file`: input merged count-ratio TSV.
- `-o, --ofile`: output TSV.
- `-n, --name`: modeled ratio columns, comma separated.
- `-t, --time`: matching time points, comma separated.
- `-i, --fc, --prior-col`: one of the modeled ratio columns, used to build
  the fold-change prior.

Important optional arguments:

- `-g, --gene`: barcode grouping column; default `gene`.
- `-e, --ef`: minimum editing efficiency; default `0`.
- `-a, --count`: low-count threshold; default `0`.
- `-r, --ratio`: drop a row when the fraction of modeled values at or below
  `--count` reaches this value; default `0.6`.
- `--nt-pattern`: regex identifying the pure Non-target reference set.
- `--control-pattern`: regex annotating control groups.
- `--prior-sigma`: shared Normal prior standard deviation; default `0.1`.
- `--obs-sigma`: fixed observation sigma. If omitted, it is estimated from
  barcode residuals.
- `--advi-steps`, `--draws`, `--learning-rate`: ADVI settings.

Run `bayesian-analysis --help` for the complete option list.

The legacy options `-m/--smpmtd` and `-c/--cutoff` remain accepted so old
scripts do not fail. Only ADVI is supported, and `--cutoff` is deprecated
because the updated output reports probabilities on both sides of zero.

## Python API

```python
import pandas as pd
from bayesian_analysis import BayesianAnalysis

data = pd.read_csv("input.tsv", sep="\t")
analysis = BayesianAnalysis(
    data,
    ratio_col="count-D6_ratio,count-D14_ratio,count-D18_ratio",
    time_points="6,14,18",
    gene="gene",
    fc="count-D14_ratio",
)
result = analysis.fit()
result.to_csv("result.tsv", sep="\t", index=False)
```

The original staged API remains available:

```python
analysis.preprocess_data()
analysis.calculate_means()
analysis.run_bayesian_model()
result = analysis.post_process_results()
```

## Main output columns

- `beta`: posterior mean of the modeled effect.
- `post_prob_up`, `post_prob_down`: posterior sign probabilities.
- `z_nt`: score standardized using the Non-target reference set.
- `p_nt_upper`, `p_nt_lower`: empirical one-sided tail probabilities using
  Non-target guides.
- `z_all`, `p_all_upper`, `p_all_lower`: corresponding all-group statistics.
- `prior_mu`, `prior_sigma`, `obs_sigma`: model values retained for auditing.
- `control`, `nt`: control-set annotations.

For compatibility with the first GitHub class release, the result also
contains `k_estimate`, `k_zscore`, `p_post`, and `p_value` aliases. New code
should use the `beta`, posterior-probability, and explicit tail columns.
