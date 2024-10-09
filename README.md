Data and Scripts used to generate all the figures and tables for the ICPE25 paper "Evaluating Performance and Energy Efficiency of Emerging HPC Processors for a Coupled Scientific Simulation"

# Content
This repository contains:
1. Benchmarking results for runtime, energy and hardware events. Located in ``results/``.
2. Scripts to generate Figures and Tables. Located in ``scripts/``
3. Snakemake build system to automaticaly generate the figures and tables

# Requirements
- Poetry [[https://python-poetry.org/]]

# How to use
1. Download repository
2. Run ``poetry install``
3. Run ``poetry run snakemake -c 1`` 

Tables and Figures are generated and placed in ``artifacts/``
