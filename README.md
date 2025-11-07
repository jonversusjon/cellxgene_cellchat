# CellChat Analysis Scripts

This repository tracks the analysis scripts used for CellChat runs on CELLxGENE data.

## Run History

### Run: cellchat_runs_all_5k_disease_normal_1256434
- **Git Commit**: `bda1d8e` (Initial commit)
- **Date**: November 7, 2025
- **Status**: Running
- **Description**: CellChat analysis on multiple adult tissues with disease/normal comparison at 5k cells

## Repository Structure

- `scripts/`: Main analysis scripts (tracked by Git)
  - `run_cellchat_batch.R`: Main CellChat batch processing script
  - `cellchat_functions.R`: Helper functions for CellChat analysis
  - `meta_cell_types.R`: Cell type metadata processing
  - `run_cellchat_all.sbatch`: SLURM batch job script
  - `cellchat_all_5k_mutli_analysis_with_disease.yaml`: Configuration file

- `cellchat_runs_*/`: Output directories (NOT tracked by Git)
- `data_helpers/`: Data processing utilities
- `notebooks/`: Jupyter notebooks for exploration

## Connecting Code to Outputs

To link a specific output run to the code version used:
1. Note the git commit hash when starting a run: `git rev-parse HEAD`
2. Record it in this README or in the output directory
3. Use `git show <commit-hash>` to view the exact code used for that run

## Setup

To push this repository to GitHub:
```bash
# Create a new repository on GitHub, then:
git remote add origin https://github.com/YOUR_USERNAME/YOUR_REPO_NAME.git
git branch -M main
git push -u origin main
```
