# PyCNO

## Description

Welcome to PyCNO - Python Computational Nuclear Oncology! This repository is centered around computational methods for predicitive dosimetry. Currently, this repository contains one prostate specific membrane antigen (PSMA) physiologically based pharmacokinetic (PBPK) model. We are building this framework to support multiple models and futher resources for simulating and fitting radiopharmaceutical biodistribution.

This repository contains a **python** implementation of physiologically based pharmacokinetic (PBPK) modeling of radiopharmaceutical therapies as developed in the following paper:

["Physiologically based radiopharmacokinetic (PBRPK) modeling to simulate and analyze radiopharmaceutical therapies: studies of non-linearities, multi-bolus injections, and albumin binding"
](https://link.springer.com/article/10.1186/s41181-023-00236-w)

EJNMMI Radiopharmacy and Chemistry, vol. 9, pp. 6, 2024

![SimBiology Implementation](/PBPK_model.png)

The original model was implemented in **MATLAB SimBiology** and can be found in the "matlab" directory. It has also been exported in the Systems Biology Markup Language (SBML) using the SBMLConversion.m script. The functions.py script contains some functions for using the SBML file in Python along with example usage in example_script.py. The motivation for having the model in Python is twofold:

1. Open Source - Can be run by anyone without a MATLAB Lisence
2. Easily integrated in workflows with other Python packages

That being said, we acknowledge the advantages of MATLAB SimBiology (namely the model GUI), and ourselves do a lot of model modifications in SimBiology, and then export the final product for integration into other workflows in Python. As such, we have included a short script in the matlab folder for seamless conversion from SimBiology to an SBML version compatible with PyCNO.

Finally, this is a new repository so we are making improvemnts. Check back for updates!

## Installation

1. Clone the repository:
   ```sh
   git clone https://github.com/qurit/PyCNO
   ```
2. Install:
   
   First navigate to cloned repository.
   ```sh
   conda create -n pycno_env python=3.12
   conda activate pycno_env
   pip install -e .
   ```

## Usage

To get started in **python**, see the examples folder.

To get started in **MATLAB Simbiology**, open the .sbproj file in the "matlab" directory. This will open a MATLAB based GUI where the model can be explored, edited, and run. If you wish to export your altered MATLAB Simbiology model to be run in python, use the export_model.m script as found in the "matlab" directory.

## **Citation**

If you use this repository, please cite the original paper:

> **Fele-Paranj, A., Saboury, B., Uribe, C., & Rahmim, A.**  
> *Physiologically based radiopharmacokinetic modeling to simulate and analyze radiopharmaceutical therapies: studies of non-linearities, multi-bolus injections, and albumin binding*.  
> **EJNMMI Radiopharmacy and Chemistry, vol. 9, pp. 6, 2024**  
> DOI: [10.1186/s41181-023-00236-w](https://doi.org/10.1186/s41181-023-00236-w) 
