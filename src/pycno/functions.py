import libsbml
import numpy as np
import roadrunner
from pathlib import Path
import os
import multiprocessing as mp
mp.set_start_method("spawn", force=True)
from tqdm.contrib.concurrent import process_map



# NEW: module-level worker to avoid pickling Model (and SWIG objects)
def parameter_sweep_worker(args):
    sbml_string, stop, steps, parameter_ids, observables, time, swept_values = args
    document = libsbml.readSBMLFromString(sbml_string)
    sbml_model = document.getModel()
    for i, pid in enumerate(parameter_ids):
        sbml_model.getParameter(pid).setValue(float(swept_values[i]))

    sbml_string = libsbml.writeSBMLToString(document)
    rr = roadrunner.RoadRunner(sbml_string)
    result = rr.simulate(0, stop, steps)

    TAC = np.zeros((len(time), len(observables)))
    for i, observable in enumerate(observables):
        TAC[:, i] = (sum_region(observable, 'Hot', result, sbml_model, rr))
    return TAC

class ModelError(Exception):
    pass

class Model():
    """
    Initializes SBML model.

    Args:
        model_name (str): Name of model in PyCNO or path to SBML file
        parameters (dict): Parameter input values
        compartment_volumes (dict): Compartment volumes in L
        hotamount (float): Hot ligand amount in nmol
        coldamount (float): Cold ligand amount in nmol
    """
    _watched_attrs = {"model_name", "hotamount", "coldamount", "parameters", "compartment_volumes"}

    def __init__(self, model_name, hotamount=None, coldamount=None, parameters=None, compartment_volumes=None):
        self._initializing = True  # prevent recursion during init

        self.model_name = model_name
        self.hotamount = hotamount
        self.coldamount = coldamount
        self.parameters = parameters
        self.compartment_volumes = compartment_volumes

        self.sbml_model = None
        self.document = None

        self._initializing = False
        self.initialize_sbml_model()

    def __setattr__(self, name, value):
        super().__setattr__(name, value)
        if getattr(self, "_initializing", False):
            return
        if name in getattr(self, "_watched_attrs", set()):
            self.initialize_sbml_model()

    def initialize_sbml_model(self):
        if os.path.exists(self.model_name):
            model_path = self.model_name
        else:
            module_path = Path(__file__).parent
            model_path = module_path / "models" / f"{self.model_name}.sbml"
            if os.path.exists(model_path) is False:
                raise FileNotFoundError(f"Model {self.model_name} not found.")

        reader = libsbml.SBMLReader()
        document = reader.readSBML(model_path)
        sbml_model = document.getModel()

        if self.parameters:
            set_parameter_values(sbml_model, self.parameters)

        if self.compartment_volumes:
            set_compartment_sizes(sbml_model, self.compartment_volumes)

        if self.hotamount is not None:
            species_list = None
            for comp in sbml_model.getListOfCompartments():
                if 'Vein' == comp.getName():
                    species_list = {'Vein.Hot': self.hotamount,
                                    'Vein.Cold': self.coldamount}
                    break
                elif 'Blood' == comp.getName():
                    species_list = {'Blood.Hot': self.hotamount,
                                    'Blood.Cold': self.coldamount}
                    break
            if species_list is None:
                raise ModelError(f"No blood or vein compartment found for injection.")
                    
            set_species_values(sbml_model, species_list)

        self.sbml_model = sbml_model
        self.document = document
    
    def save_sbml(self, path):
        """
        Saves current SBML model to file.

        Args:
            path (str): Path to save SBML file
        """
        libsbml.writeSBMLToFile(self.document, path)

    def simulate(self,
                stop: int = 60, 
                steps: int = 100,
                observables: list = None,
                swept_parameters: list = None,
                swept_values: list = None
                ):
        """
        Simulates SBML model.

        Args:
            model_name (str): Name of model in PyCNO or path to SBML file
            observables (list): Regions to output
            stop (int): Simulation end time in minutes
            steps (int): Number of simulation steps
            hotamount (float): Hot ligand amount in nmol
            coldamount (float): Cold ligand amount in nmol
            parameters (dict): Parameter input values
            compartment_volumes (dict): Compartment volumes in L
            swept_parameters (list): Parameters to sweep over
            swept_values (list): Values to sweep over
            
        Returns:
            TACs[n_curves, n_steps, n_observables]: Time activity curves in units of MBq.
        
        """

        num_curves = 1
        time = np.linspace(0, stop, steps)

        if swept_parameters:
            parameter_ids = []
            for parameter in swept_parameters:
                parameter_ids.append(get_parameter_id(self.sbml_model, parameter))
            sbml_string = libsbml.writeSBMLToString(self.document)

            args = [(sbml_string, stop, steps, parameter_ids, observables, time, values) for values in swept_values]
            # Use threads to avoid os.fork with JAX's multithreading
            results = process_map(parameter_sweep_worker, args, max_workers=os.cpu_count(), chunksize=1)
            TACs = np.stack(results, axis=0)
        else:
            sbml_string = libsbml.writeSBMLToString(self.document)
            rr = roadrunner.RoadRunner(sbml_string)
            result = []
            result.append(rr.simulate(0, stop, steps))
            TACs = np.zeros((num_curves, len(time), len(observables)))
            for i, observable in enumerate(observables):
                TACs[0,:,i] = (sum_region(observable, 'Hot', result[0], self.sbml_model, rr))

        return time, TACs


def set_parameter_values(sbml_model, parameter_dict_in):
    parameter_dict = parameter_dict_in.copy() 
    for parameter in sbml_model.getListOfParameters():
        for parameter_name in list(parameter_dict.keys()):
            if parameter.getName() == parameter_name:
                parameter.setValue(float(parameter_dict.pop(parameter_name)))
    if parameter_dict:
        [print(f"Parameter {parameter_name} not found in the model.") for parameter_name in parameter_dict.keys()]

def set_species_values(sbml_model, species_dict_in):
    species_dict = species_dict_in.copy()
    for species_name in list(species_dict.keys()):
        for compartment in sbml_model.getListOfCompartments():
            if compartment.getName() == species_name.split('.')[0]:
                break
        for species in sbml_model.getListOfSpecies():
            if species.getCompartment() == compartment.getId() and species.getName() == species_name.split('.')[1]:
                species.setInitialAmount(species_dict.pop(species_name))
                break
    if species_dict:
        [print(f"Species {species_name} not found in the model.") for species_name in species_dict.keys()]

def set_compartment_sizes(sbml_model, compartment_dict_in):
    compartment_dict = compartment_dict_in.copy()
    for comp in sbml_model.getListOfCompartments():
        if comp.getName() in compartment_dict:
            comp.setSize(float(compartment_dict[comp.getName()]))
            compartment_dict.pop(comp.getName())
    if compartment_dict:
        for missing in compartment_dict.keys():
            print(f"Compartment {missing} not found in the model.")

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
        
def sum_region(region, species_name, result, sbml_model, rr):
    NMOL2MBQ = get_parameter(sbml_model, 'lambdaPhys') / 60 * 6.022e23 / 10**9 / 10**6
    total = np.zeros(result.shape[0])
    for compartment in sbml_model.getListOfCompartments():
            if compartment.getName().startswith(region):
                compartment_id = compartment.getId()
                compartment_size = rr[f"{compartment_id}"]
                for species in sbml_model.getListOfSpecies():
                    if species.getCompartment() == compartment_id and species_name in species.getName():
                        total += (result[f"[{species.getId()}]"] * compartment_size * NMOL2MBQ)
    return total