# Cosmic Shear Analysis Using a Model-Independent Expansion Function
This repo contains CosmoSIS' original framework and my modifications to include a model-independent expansion function.

## Structure: 
```
.
├── cosmosis-standard-library/     # CosmoSIS original library, with some modifications
├── KiDS-1000_cosebis_config_me/   # Original KiDS COSEBIs configuration
├── mod/                           # Modules for model-independent analysis
├── thesis_config/                 
│   ├── config/                    # Main pipeline .ini files for CosmoSIS
│   │   ├── best_fit_values        # Best-fit values from Nelder-Mead
│   │   ├── evaluators             # Jupyter notebooks for analysis            
│   ├── inputs/                    # Input for the pipeline
│   ├── kcap_modules/              # KiDS' kcap modules
│   └── outputs/                   # Output chains from CosmoSIS runs
└── WnLog/                         # Precomputed COSEBIs weight functions
```

The root directory contains:
- `cosmosis-standard-library` contains the original code CosmoSIS library, with some modifications from me.
- `KiDS-1000_cosebis_config_me` contains the original configuration from KiDS
- `mod` contains some modules for model-independent analysis
- `thesis_config` contains my configurationss
- `WnLog` contains the calculated COSEBIs weight functions for theta covering (0.5, 300) and (60, 300). Use tools in https://github.com/KiDS-WL/Cat_to_Obs_K1000_P1.git and  https://github.com/maricool/2pt_stats.git to get a calculate the range that you want.

Inside `thesis_config`:
- `best_fit_values` contains best-fitting parameters from Nelder-Mead
- `config` contains pipeline files
- `evaluators` contains jupyter notebooks for various purposes

## Installation
Refer to https://cosmosis.readthedocs.io/en/latest/intro/installation.html for installing cosmosis

## Running the pipeline
Before running the pipeline the cosmosis configuration must be done:
- conda activate ./env
- source cosmosis-configure

To run the pipeline: `cosmosis <pipeline_file>.ini`

The important pipeline files are:
- `KiDS-1000_cosebis_config_me/config/pipeline_multinest.ini`
This contain's KiDS' original pipeline. The cosmosis-standard-library directory should be the original one, which is not in this repo (you can find it in cosmosis' repo).
- `thesis_config/config/pipeline_multinest.ini`
This is the pipeline for testing the LCDM model.
- `thesis_config/config/wl_pipe.ini`
This is the pipeline for using only weak lensing data.
- `thesis_config/config/joint_pipe.ini`
This is the pipeline for a joint analysis with other datasets (SN, BAO, CMB, QSO).

Modifications or parameters that enter the pipeline can be made in `values.ini` and types of priors can be set up in `priors.ini` files.

## Evaluation
After running the pipeline, the output is saved in `thesis_config/outputs/`. Go to the directory of the run (e.g. `thesis_config/outputs/run1`), copy `output_multinest.txt` and paste it to `thesis_config/config/tester/outputs/`. Use [`marginal_mode.ipynb`](thesis_config/config/evaluators/marginal_mode.ipynb) to obtain best-fitting values and uncertainties from the chain, and to plot the 2D parameter contours.

