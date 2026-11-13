from pathlib import Path

from setuptools import find_packages, setup


ROOT = Path(__file__).parent


setup(
    name="bayesian-analysis",
    version="0.2.0",
    packages=find_packages(),
    python_requires=">=3.9",
    install_requires=[
        "numpy",
        "pandas",
        "pymc>=5",
        "pytensor",
        "click",
    ],
    entry_points={
        "console_scripts": [
            "bayesian-analysis=bayesian_analysis.cli:main",
        ],
    },
    author="WangLab",
    description="Editing-activity-adjusted Bayesian modeling for CRISPR-BE screens.",
    long_description=(ROOT / "README.md").read_text(encoding="utf-8"),
    long_description_content_type="text/markdown",
)
