# Simulation Data Repository

This repository contains molecular dynamics simulation trajectories and analyzed data.

## Data
### inp

The input files for CafeMol

#### AdK_force

Input files for AdK with force simulation, containing 2E8 MD steps

#### MFPT_c2o\MFPT_o2c

Input files for MFPT calculations, containing fewer MD steps

### turnover rate

Turnover rate values obtained from simulations. For details, please refer to README.md in turnover_rate/

### MFPT
The Mean First Passage Time (MFPT) calculated from short trajectories.
'o2c' represents trajectories starting from the open state with no ligand binding
'c2o' represents trajectories starting from the closed state with 2 ADP molecules bound
'phy' represents analysis focused only on conformational changes

### average_time
The average time calculated from trajectories
'o2c' represents the time required to form the catalytically competent state
'c2o' represents the time required for product release
'all' represents the time interval between two adjacent cycles

## Trajectories

Trajectories obtained from simulations

### MFPT_AdK_***

Short trajectories for MFPT calculations

## cafemol-ake-7_rmsd.tar.gz

The modified version of CafeMol
