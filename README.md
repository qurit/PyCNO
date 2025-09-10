# PyPBPK

## Description

This repository contains a **python** implementation of physiologically based pharmacokinetic (PBPK) modeling of radiopharmaceutical therapies as developed in the following paper:

["Physiologically based radiopharmacokinetic (PBRPK) modeling to simulate and analyze radiopharmaceutical therapies: studies of non-linearities, multi-bolus injections, and albumin binding"
](https://link.springer.com/article/10.1186/s41181-023-00236-w)

EJNMMI Radiopharmacy and Chemistry, vol. 9, pp. 6, 2024

The original model was implemented in **MATLAB SimBiology** and can be found in the "matlab" directory. It has also been exported in the Systems Biology Markup Language (SBML) using the SBMLConversion.m script. The functions.py script contains some functions for using the SBML file in Python along with example usage in example_script.py. The motivation for having the model in Python is twofold:

1. Open Source - Can be run by Anyone wihtout a MATLAB Lisence
2. Easily integrated in workflows with other Python packages

That being said, we acknowledge the advantages of MATLAB SimBiology (namely the model GUI), and ourselves do a lot of model modifications in SimBiology, and then export the final product for integration into other workflows in Python.

## Installation

1. Clone the repository:
   ```sh
   git clone https://github.com/your-username/PyPBPK.git
   ```
2. Install dependencies:
   ```sh
   pip install -r requirements.txt
   ```

## Usage

To get started in **python**, run the `runPBPK.py` script. This is the most basic implementation that will generate time activity curves (TACs) for a single virtual patient. The default parameters approximate a population average. These parameters can be altered at the bottom of the script. To generate multiple patients by varying these parameters, use the `PBPKSweep.py` script and enter a range of values in the "args" list.

To get started in **MATLAB Simbiology**, open the .sbproj file in the "matlab" directory. This will open a MATLAB based GUI where the model can be explored, edited, and run. If you wish to export your altered MATLAB Simbiology model to be run in python, use the export_model.m script as found in the "matlab" directory.

## **Citation**

If you use this repository, please cite the original paper:

> **Fele-Paranj, A., Saboury, B., Uribe, C., & Rahmim, A.**  
> *Physiologically based radiopharmacokinetic modeling to simulate and analyze radiopharmaceutical therapies: studies of non-linearities, multi-bolus injections, and albumin binding*.  
> **EJNMMI Radiopharmacy and Chemistry, vol. 9, pp. 6, 2024**  
> DOI: [10.1186/s41181-023-00236-w](https://doi.org/10.1186/s41181-023-00236-w) 
