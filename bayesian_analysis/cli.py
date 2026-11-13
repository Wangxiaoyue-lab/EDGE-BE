"""Command-line interface for EDGE-BE Bayesian analysis."""

from pathlib import Path

import click
import pandas as pd

from .analysis import BayesianAnalysis


@click.command(context_settings={"help_option_names": ["-h", "--help"]})
@click.option("--file", "-f", required=True, type=click.Path(exists=True, dir_okay=False), help="Input merged count-ratio TSV file.")
@click.option("--ofile", "-o", required=True, type=click.Path(dir_okay=False), help="Output TSV file.")
@click.option("--name", "-n", required=True, help="Comma-separated modeled ratio columns.")
@click.option("--time", "-t", required=True, help="Comma-separated time points matching --name.")
@click.option("--gene", "-g", default="gene", show_default=True, help="Column used to group barcodes.")
@click.option("--ef", "-e", type=float, default=0.0, show_default=True, help="Minimum editing efficiency.")
@click.option("--prior-col", "--fc", "-i", required=True, help="Modeled ratio column used to construct the prior.")
@click.option("--count", "-a", type=float, default=0.0, show_default=True, help="Low-count ratio threshold.")
@click.option("--ratio", "-r", type=float, default=0.6, show_default=True, help="Maximum allowed fraction of low-count time points.")
@click.option("--nt-pattern", default="Non-target", show_default=True, help="Regex identifying Non-target groups.")
@click.option("--control-pattern", default=r"control|AAVS1|Non-target", show_default=True, help="Regex identifying control groups.")
@click.option("--scale-factor", type=float, default=30_000_000.0, show_default=True, help="Scale converting ratio columns to proportions.")
@click.option("--epsilon", type=float, default=5.0, show_default=True, help="Pseudocount used for prior fold change.")
@click.option("--prior-sigma", type=float, default=0.1, show_default=True, help="Normal prior standard deviation.")
@click.option("--obs-sigma", type=float, default=None, help="Observation sigma; estimated from barcode residuals when omitted.")
@click.option("--min-obs-sigma", type=float, default=1e-8, show_default=True, help="Lower bound for estimated observation sigma.")
@click.option("--advi-steps", type=int, default=20_000, show_default=True, help="Number of ADVI optimization steps.")
@click.option("--draws", type=int, default=10_000, show_default=True, help="Posterior draws from the variational approximation.")
@click.option("--learning-rate", type=float, default=0.01, show_default=True, help="ADVI learning rate.")
# Retain old command lines while making the updated model's fixed behavior clear.
@click.option("--smpmtd", "-m", default="advi", hidden=True)
@click.option("--cutoff", "-c", type=float, default=None, hidden=True)
def main(
    file,
    ofile,
    name,
    time,
    gene,
    ef,
    prior_col,
    count,
    ratio,
    nt_pattern,
    control_pattern,
    scale_factor,
    epsilon,
    prior_sigma,
    obs_sigma,
    min_obs_sigma,
    advi_steps,
    draws,
    learning_rate,
    smpmtd,
    cutoff,
):
    """Fit EDGE-BE and write guide/group-level posterior statistics."""
    if smpmtd.lower() != "advi":
        raise click.UsageError("The updated model currently supports only --smpmtd advi.")
    if cutoff is not None:
        click.echo(
            "Warning: --cutoff is deprecated; the output now reports posterior "
            "probabilities above and below zero.",
            err=True,
        )

    data = pd.read_csv(file, sep="\t")
    analysis = BayesianAnalysis(
        data=data,
        ratio_col=name,
        time_points=time,
        gene=gene,
        ef=ef,
        fc=prior_col,
        count=count,
        ratio=ratio,
        nt_pattern=nt_pattern,
        control_pattern=control_pattern,
        scale_factor=scale_factor,
        epsilon=epsilon,
        prior_sigma=prior_sigma,
        obs_sigma=obs_sigma,
        min_obs_sigma=min_obs_sigma,
        advi_steps=advi_steps,
        draws=draws,
        learning_rate=learning_rate,
    )
    result = analysis.fit()

    output_path = Path(ofile)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(output_path, sep="\t", index=False)
    click.echo(f"Wrote {len(result)} rows to {output_path}")


if __name__ == "__main__":
    main()
