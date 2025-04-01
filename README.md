# Simulation Data Repository

This repository contains molecular dynamics simulation trajectories and analyzed data.

## Data

Each subdirectory contains two versions of data files and corresponding plotting scripts:

1. Raw data files with names starting with 'raw_' and ending with '.dat'
2. Readable data files with names starting with 'readable_' and ending with '.md'

please refer to the readable data.

The dataset comprises three types of data:

### turnover rate

Turnover rate values obtained from simulations 
The unit is $\rm{ms^{-1}}$

### MFPT
The Mean First Passage Time (MFPT) calculated from short trajectories.
'o2c' represents trajectories starting from the open state with no ligand binding
'c2o' represents trajectories starting from the closed state with 2 ADP molecules bound
'phy' represents analysis focused only on conformational changes
The unit is ms

### average_time
The average time calculated from trajectories
'o2c' represents the time required to form the catalytically competent state
'c2o' represents the time required for product release
'all' represents the time interval between two adjacent cycles
The unit is ms
 

### script
Analysis scripts used for data processing.

## inp

The input files for CafeMol

### AdK_force

Input files for AdK with force simulation, containing 2E8 MD steps

### MFPT_c2o\MFPT_o2c

Input files for MFPT calculations, containing fewer MD steps

## Trajectories

Trajectories obtained from simulations

### MFPT_AdK_Pc028_m0300_s148_f10

Short trajectories for MFPT calculations

### rakesi_Pc028_m0300_f0

Trajecories without force applied with $P_{close} \approx 0.28$ (wild-type)

### rakesi_Pc099_m0300_f0

Trajecories without force applied with $P_{close} \approx 0.99$ (extremely closed)

### rakesi_Pc028_m0300_f0

Trajecories without force applied with $P_{close} \approx 0.01$ (extremely open)

### akesi_Pc028_m0300_f10

Trajecories with 10pN force applied with $P_{close} \approx 0.28$ (wild-type)

### akesi_dna_Pc028_m0300_bp70

Trajecories with 70bp DNA applied with $P_{close} \approx 0.28$ (wild-type)