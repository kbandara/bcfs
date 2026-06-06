## EEG analysis for a breaking continuous flash suppression (bCFS) dataset**

> ***Kav Bandara** University of Melbourne, 2025*

This repo contains the complete analysis pipeline for an EEG study of conscious visual perception using the **breaking continuous flash suppression (bCFS)** paradigm [(McFadyen et al., 2022)](https://www.frontiersin.org/articles/10.3389/fnbeh.2022.797119/full).

## Pipeline

The analysis has four main steps:

### 1. EEG Preprocessing (`s01_prepro`)
Preprocessing occurs in two sub-stages: automated ICA artefact rejection in MNE-Python (`s01_ica.py`), followed by the remainder of preprocessing in MATLAB (`s02_preprocessing.m`).

### 2. Source Reconstruction and ROI Localisation (`s02_source_analysis`)

### 3. SVM Decoding (`s03_decoding`)

### 4. DCM and PEB (`s04_dcm`)

This stage also includes a comparison of the standard DCM input function (`spm_erp_u_standard.m`) against a custom input function (`spm_erp_u_ramping.m`). This custom script replaced the `spm_erp_u.m` function in a copy of SPM which was added to the path where relevant, e.g.: 
```matlab
    if strcmp(input_type, 'standard')
        spm_root = '.../spm12_standard/spm12';
    else
        spm_root = '.../spm12_ramping/spm12';
    end
```

## Dependencies

This analysis pipeline was run using MNE-Python, matlab, and SPM12. 

SPM12 is freely available from the [Wellcome Centre for Human Neuroimaging](https://www.fil.ion.ucl.ac.uk/spm/software/spm12/).

## HPC Usage & Parallel Processing

Some scripts (ICA/DCM scripts) are setup to accept integer arguments to index the participant list and run in parallel using SLURM job arrays, e.g.:

```bash
# one job per participant
sbatch --array=1-30 ica.sh
```

 The MATLAB preprocessing and decoding scripts support parallel processing via `parfor` within a single multi-core job, which can be specified at the top of the script. 
 
A testing block at the top of each function (commented out by default) allows running a single participant interactively:

```matlab
% input_type = 'standard'; p_name = 'S01';
```
