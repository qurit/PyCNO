import libsbml
import numpy as np
import roadrunner
from pathlib import Path
from tqdm.contrib.concurrent import process_map

def set_parameter_values(sbml_model, parameter_list):
    for parameter in sbml_model.getListOfParameters():
        for parameter_name in list(parameter_list.keys()):
            if parameter.getName() == parameter_name:
                parameter.setValue(float(parameter_list.pop(parameter_name)))
    if parameter_list:
        [print(f"Parameter {parameter_name} not found in the model.") for parameter_name in parameter_list.keys()]

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

def get_parameter(sbml_model, parameter_name):
    for parameter in sbml_model.getListOfParameters():
        if parameter.getName() == parameter_name:
                return parameter.getValue()
    print(f"Parameter {parameter_name} not found in the model.")

def get_parameter_id(sbml_model, parameter_name):
    for parameter in sbml_model.getListOfParameters():
        if parameter.getName() == parameter_name:
                return parameter.getId()
    print(f"Parameter {parameter_name} not found in the model.")

def get_species(species_name, result, rr, sbml_model):
    NMOL2MBQ = get_parameter(sbml_model, 'lambdaPhys') / 60 * 6.022e23 / 10**9 / 10**6
    for compartment in sbml_model.getListOfCompartments():
            if compartment.getName() == species_name.split('.')[0]:
                compartment_size = compartment.getSize()
                break
    for species in sbml_model.getListOfSpecies():
        if species.getName() == species_name.split('.')[1]:
            return (result[f"[{species.getId()}]"] * compartment_size * NMOL2MBQ)
        
def sum_region(region, species_name, result, sbml_model):
    NMOL2MBQ = get_parameter(sbml_model, 'lambdaPhys') / 60 * 6.022e23 / 10**9 / 10**6
    total = np.zeros(result.shape[0])
    for compartment in sbml_model.getListOfCompartments():
            if region in compartment.getName():
                compartment_id = compartment.getId()
                compartment_size = compartment.getSize()
                for species in sbml_model.getListOfSpecies():
                    if species.getCompartment() == compartment_id and species_name in species.getName():
                        total += (result[f"[{species.getId()}]"] * compartment_size * NMOL2MBQ)
                        print(compartment_size)

    return total

def parameter_sweep(sbml_string, start, stop, steps, parameter_ids, swept_values):
    document = libsbml.readSBMLFromString(sbml_string)
    sbml_model = document.getModel()
    for i, id in enumerate(parameter_ids):
        sbml_model.getParameter(id).setValue(swept_values[i])

    sbml_string = libsbml.writeSBMLToString(document)
    rr = roadrunner.RoadRunner(sbml_string)
    result = rr.simulate(start, stop, steps)
    return result

def multicore_parameter_sweep(args):
        return parameter_sweep(*args)

def run_model(model_name: str = None,
            start: int = 0, 
            stop: int = 60, 
            steps: int = 100, 
            hotamount: float = 10, 
            coldamount: float = 100, 
            parameters: dict = None,
            observables: list = None,
            swept_parameters: list = None,
            swept_values: list = None
            ):
    
    num_curves = 1
    module_path = Path(__file__).parent
    model_path = module_path / "models" / f"{model_name}.sbml"
    reader = libsbml.SBMLReader()
    document = reader.readSBML(model_path)
    sbml_model = document.getModel()

    # Set parameters based on matlab variants (temporary fix)
    varinat_values = {'Rden_Tumor1': 57,
                      'Rden_Tumor2': 19,
                      'Rden_SG': 38,
                      'Rden_Kidney': 14,
                      'lambdaRel_Tumor1': 1.5e-04,
                      'lambdaRel_Tumor2': 1.5e-04,
                      'lambdaRel_SG': 4.2e-04,
                      'lambdaRel_Kidney': 2.9e-04,
                      'TER_Kidney': 0.2,
                      'bodySurfaceArea': 1.9,
                      'bodyWeight': 100,
                      'bodyHeight': 160,
                      'f_SG': 0.0740,
                      'R0_TumorRest': 13,
                      
                      'kPSAlb_Tumor1': 0,
                      'kPSAlb_Tumor2': 0,
                      'kPSAlb_TumorRest': 0,
                      'k_on_toAlb': 0,
                      'k_off_toAlb': 0,
                      'AlbuminDen': 0,
                      'K_D_Alb': 1}
    set_parameter_values(sbml_model, varinat_values)

    if parameters:
        set_parameter_values(sbml_model, parameters)

    if hotamount:
        if coldamount:
            species_list = {'Vein.Hot': hotamount,
                            'Vein.Cold': coldamount}
        else:
            species_list = {'Vein.Hot': hotamount}
        set_species_values(sbml_model, species_list)

    if swept_parameters:
        parameter_ids = []
        for parameter in swept_parameters:
            parameter_ids.append(get_parameter_id(sbml_model, parameter))
        sbml_string = libsbml.writeSBMLToString(document)

        args = [(sbml_string, start, stop, steps, parameter_ids, values) for values in swept_values]
        num_curves = len(args)
        result = process_map(multicore_parameter_sweep, args)

    else:
        sbml_string = libsbml.writeSBMLToString(document)
        rr = roadrunner.RoadRunner(sbml_string)
        result = []
        result.append(rr.simulate(start, stop, steps))

    if not observables:
        observables = ['Tumor1', 'Tumor2', 'Kidney', 'Heart', 'SG', 'Bone', 'TumorRest', 'Spleen', 'Liver', 'Prostate', 'GI', 'Rest', 'Skin', 'Muscle', 'Brain', 'RedMarrow', 'Lungs', 'Adipose']
    
    time = np.linspace(start, stop, steps)
    TACs = np.zeros((num_curves, len(time), len(observables)))
    for curve in range(num_curves):
        for i, observable in enumerate(observables):
            TACs[curve,:,i] = (sum_region(observable, 'Hot', result[curve], sbml_model))

    return time, TACs