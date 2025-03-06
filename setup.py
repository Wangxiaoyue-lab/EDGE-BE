# setup.py

from setuptools import setup, find_packages

setup(
    name='bayesian-analysis',
    version='0.1',
    packages=find_packages(),
    install_requires=[
        'pandas',
        'numpy',
        'click',
        'pymc',
        'pytensor',
        'matplotlib',
        'arviz',
        'scipy'
    ],
    entry_points={
        'console_scripts': [
            'bayesian-analysis=bayesian_analysis.cli:main',
        ],
    },
    author='WangLab',
    description='A Bayesian analysis tool for CRISPR-BE screening data.',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
)
