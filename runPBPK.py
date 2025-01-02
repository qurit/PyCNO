#%%
import libsbml
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import roadrunner
# %%
def set_parameter_values(sbml_model, parameter_list):
    for parameter in sbml_model.getListOfParameters():
        for paramter_name in list(parameter_list.keys()):
            if parameter.getName() == paramter_name:
                parameter.setValue(parameter_list.pop(paramter_name))
    if parameter_list:
        [print(f"Parameter {paramter_name} not found in the model.") for paramter_name in parameter_list.keys()]

def set_species_values(sbml_model, species_list):
    for species_name in list(species_list.keys()):
        for compartment in sbml_model.getListOfCompartments():
            if compartment.getName() == species_name.split('.')[0]:
                break
        for species in sbml_model.getListOfSpecies():
            if species.getCompartment() == compartment.getId() and species.getName() == species_name.split('.')[1]:
                species.setInitialAmount(species_list.pop(species_name))
                break
    if species_list:
        [print(f"Species {species_name} not found in the model.") for species_name in species_list.keys()]


def get_species(species_name, result, rr, sbml_model):
    for compartment in sbml_model.getListOfCompartments():
            if compartment.getName() == species_name.split('.')[0]:
                compartment_size = getattr(rr, compartment.getId())
                break
    for species in sbml_model.getListOfSpecies():
        if species.getName() == species_name.split('.')[1]:
            return (result[f"[{species.getId()}]"] * compartment_size * NMOL2MBQ)

# %%
if __name__ == '__main__':
    model_path = 'model/PSMAModel.sbml'
    reader = libsbml.SBMLReader()
    document = reader.readSBML(model_path)
    sbml_model = document.getModel()
    DECAY_CONSTANT = 0.00632
    NMOL2MBQ = DECAY_CONSTANT / 60 * 6.022e23 / 10**9 / 10**6


    parameter_list = {'lambdaPhys': DECAY_CONSTANT,
                     'Tumor1Volume': 0.01,
                     'lambdaRel_Tumor1': 1.5e-4,
                     'Tumor1VolumeCoeff': 1,
                     'Rden_Tumor2': 0,
                     'lambdaRel_Tumor2': 1.5e-4,
                     'lambdaRel_SG': 4.2e-4,
                     'lambdaRel_Kidney': 2.9e-4,
                     'TER_Kidney': 0.2,
                     'bodySurfaceArea': 1.9,
                     'bodyWeight': 100,
                     'bodyHeight': 160,
                     'f_SG': 0.074,
                     'R0_TumorRest': 13,
                     'kPSAlb_Tumor1': 0,
                     'kPSAlb_Tumor2': 0,
                     'kPSAlb_TumorRest': 0,
                     'k_on_toAlb': 0,
                     'k_off_toAlb': 0,
                     'AlbuminDen': 0,
                     'K_D_Alb': 1,
                     'Rden_Tumor1': 100,
                     'Rden_SG': 38,
                     'Rden_Kidney': 14,
                     }
    set_parameter_values(sbml_model, parameter_list)

    species_list = {'Vein.Hot': 0.005,
                    'Vein.Cold': 3.995}
    set_species_values(sbml_model, species_list)

    sbml_string = libsbml.writeSBMLToString(document)
    rr = roadrunner.RoadRunner(sbml_string)
    result = rr.simulate(0, 500, 1000)

# %%

def sum_region(region, species_name, result, rr, sbml_model):
    for compartment in sbml_model.getListOfCompartments():
            if region in compartment.getName():
                compartment_id = compartment.getId()
                compartment_size = getattr(rr, compartment_id)
                for species in sbml_model.getListOfSpecies():
                    if species.getCompartment() == compartment_id and species_name in species.getName():
                        total = total + (result[f"[{species.getId()}]"] * compartment_size * NMOL2MBQ)

time = result['time']
tumor1 = sum_region('Tumor1', 'Hot', result, rr, sbml_model)
# %%