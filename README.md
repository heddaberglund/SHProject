# Modelling PET Scanning for Breast Cancer Screening (Senior Honours Project)

All data and code I have produced/modified for my Senior Honours project on modelling PET scanning for breast cancer screening.
The PET Scanner simulation was performed with Geant4 using code adapted from [bwynneHEP/SimplePetScanner](https://github.com/bwynneHEP/SimplePetScanner).

## Reproducing the Geant4 Simulation

To reproduce the Geant4 simulation, you must first have the full [bwynneHEP/SimplePetScanner](https://github.com/bwynneHEP/SimplePetScanner) project downloaded and running correctly. This repository only provides the modified components, not the full simulation framework.

Once the base project is set up, replace **DetectorConstruction.cpp** and **LinearSourceAction.cpp** in your Geant4 project with the versions provided in the **Geant4/** directory.

**Note:**  
The current implementation of **LinearSourceAction.cpp** distributes decays uniformly throughout the entire phantom based on its volume.  
To restrict decays to specific regions of the phantom, remove the relevant conditional statements and adjust the numerical values defining the box dimensions or hemisphere radius.

# Data
All raw CSV files are located in the **data/** folder. 

# Code
To reproduce all figures and numerical outputs:

```bash
python3 data_analysis.py
