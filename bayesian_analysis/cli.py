# bayesian_analysis/cli.py

import click
from .analysis import BayesianAnalysis

@click.command()
@click.option("--file","-f", help="merge stat file")
@click.option("--ofile","-o", help="output file")
@click.option("--smpmtd","-m", help="sampling method", default="advi")
@click.option("--cutoff","-c", type=float, default=-0.05)
@click.option("--name","-n", help="column name of count ratio, comma split")
@click.option("--time","-t", help="time point, comma split")
@click.option("--gene","-g", default="gene")
@click.option("--ef","-e", type=float, default=0)
@click.option("--fc","-i", help="initial fc of column")
@click.option("--count","-a", type=float, default=0)
@click.option("--ratio","-r", type=float, default=0.6)
def main(file, ofile, smpmtd, cutoff, name, time, gene, ef, fc, count, ratio):
    # Load data
    data = pd.read_csv(file, sep='\t')
    
    # Initialize BayesianAnalysis object
    ba = BayesianAnalysis(data, name.split(','), [float(t) for t in time.split(',')], gene=gene, ef=ef, fc=fc, count=count, ratio=ratio)
    
    # Preprocess data
    ba.preprocess_data()
    ba.calculate_means()
    
    # Run Bayesian model
    ba.run_bayesian_model()
    
    # Post-process results and save output
    ba.post_process_results().to_csv(ofile, sep='\t', index=False)

if __name__ == "__main__":
    main()
