#%%
import os
import numpy as np
import matplotlib.pyplot as plt
from pyCNO import runPBPK
swept_parameters = ['Tumor1Volume', 'Rden_Tumor1']
swept_values = [(Tumor1Volume,Rden_Tumor1)
                for Tumor1Volume in np.linspace(0.1, 0.5, 6)
                for Rden_Tumor1 in np.arange(0.1, 0.5, 0.1)
                ]

model_path = os.path.join(os.path.dirname(__file__), "PBPK_177Lu_PSMA.sbml")

time, TACs = runPBPK(
    model_path=model_path,
    stop=10000,
    swept_parameters=swept_parameters,
    swept_values=swept_values,
)
plt.plot(time, TACs[0])
plt.show()
# %%
import libsbml

num_curves = 1
reader = libsbml.SBMLReader()
document = reader.readSBML("PBPK_177Lu_PSMA.sbml")
sbml_model = document.getModel()

for parameter in sbml_model.getListOfParameters():
    for parameter_name in list(swept_parameters.keys()):
        if parameter.getName() == parameter_name:
            parameter.setValue(float(swept_parameters.pop(parameter_name)))
if swept_parameters:
    [print(f"Parameter {parameter_name} not found in the model.") for parameter_name in swept_parameters.keys()]