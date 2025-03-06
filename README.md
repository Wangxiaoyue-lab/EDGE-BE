# EDGE-BE
An editing-activity-adjusted dynamic modeling method for multi-time-point screening data based on cellular growth

## Table of Contents

- [Features](#features)
- [Installation](#installation)
- [Usage](#usage)

## Features

- Perform Bayesian inference on genetic datasets.
- Analyze sgRNA efficiency and its impact on growth rates.
- Visualize the results with ELBO plots and trace plots.
- Export analysis results in a tab-separated values format.

## Installation

Ensure you have Python 3.7+ installed. You can install the package directly from GitHub:
pip install git+https://github.com/Wangxiaoyue-lab/EDGE-BE.git

Or clone the repository and install it locally:

git clone https://github.com/Wangxiaoyue-lab/EDGE-BE.git

cd EDGE-BE

pip install .

## Usage

bayesian-analysis -f path/to/input/file.tsv -o path/to/output/file.tsv -m advi -c -0.05 -n "count-1_ratio,count-2_ratio" -t "1,2,3" -g gene_column -e 0 -fc fc_column -a 0 -r 0.6

Parameters

-f, --file: Path to the input file (TSV format).

-o, --ofile: Path to save the output file.

-m, --smpmtd: Sampling method (advi by default).

-c, --cutoff: Cutoff value for functional score (default is -0.05).

-n, --name: Comma-separated list of column names representing count ratios.

-t, --time: Comma-separated list of time points.

-g, --gene: Column name containing gene information (default is gene).

-e, --ef: Cutoff value for editing efficiency (default is 0).

-fc, --fc: Column name for initial fold change.

-a, --count: Cutoff value for count ratio (default is 0).

-r, --ratio: Initial fold change value (default is 0.6).
