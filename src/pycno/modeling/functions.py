from tqdm.contrib.concurrent import process_map
import libsbml
import numpy as np
import roadrunner
from pathlib import Path
import os
import pickle
from scipy import sparse
import multiprocessing as mp
mp.set_start_method("spawn", force=True)
import diffrax
import equinox as eqx
import inspect
import jax
import jax.numpy as jnp
import pandas as pd
from pycno.utils.jax_conversion import convert_model_to_jax
from dataclasses import dataclass

# TODO: Compartment units
# TODO: Warning for parameters altered at runtime
# TODO: Error for non found observables
# TODO: Add docstrings
# TODO: Add cold species return
# TODO: Add all possible inputs for sweep
# TODO: Update watched attrs
# TODO: Dose error exceptions
# TODO: Fix FIM
# TODO: Add multidosing to jax

class ModelError(Exception):
    pass
class SimulationError(Exception):
    pass

@dataclass
class Dose():
    times: list
    targets: dict
    def set_ids(self, sbml_model):
        self.ids = []
        for species_name in list(self.targets.keys()):
            for compartment in sbml_model.getListOfCompartments():
                if compartment.getName() == species_name.split('.')[0]:
                    break
            for species in sbml_model.getListOfSpecies():
                if species.getCompartment() == compartment.getId() and species.getName() == species_name.split('.')[1]:
                    self.ids.append(species.getId())
                    break

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
    _watched_attrs = {"model_name", "hotamount",
                      "coldamount", "parameters", "compartment_volumes"}

    def __init__(self,
                 model_name,
                 parameters=None,
                 compartment_volumes=None):

        self._initializing = True

        self.model_name = model_name
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

    def __getstate__(self):
        state = self.__dict__.copy()
        for key, value in list(state.items()):
            try:
                pickle.dumps(value)
            except (TypeError, pickle.PicklingError):
                del state[key]
        return state

    def __setstate__(self, state):
        self.__dict__.update(state)

    def initialize_sbml_model(self):
        if os.path.exists(self.model_name):
            model_path = self.model_name
        else:
            module_path = Path(__file__).parent.parent
            model_path = module_path / "models" / f"{self.model_name}.sbml"
            if os.path.exists(model_path) is False:
                model_path = module_path / "models" / f"{self.model_name}"
                if os.path.exists(model_path) is False:
                    raise FileNotFoundError(f"Model {self.model_name} not found.")

        reader = libsbml.SBMLReader()
        document = reader.readSBML(model_path)
        sbml_model = document.getModel()

        if self.parameters:
            set_parameter_values(sbml_model, self.parameters)

        if self.compartment_volumes:
            set_compartment_sizes(sbml_model, self.compartment_volumes)

        self.sbml_model = sbml_model
        self.document = document

        self.NMOL2MBQ = get_parameter(self.sbml_model, 'lambdaPhys') / \
        60 * 6.022e23 / 10**9 / 10**6

    def simulate(self,
                 dose: Dose,
                 stop: int = 60,
                 steps: int = 100,
                 output_compartments: list = None,
                 return_all_subcompartments: bool = False,
                 output_parameters: list = None,
                 return_all_parameters: bool = False,
                 swept_parameters: list = None,
                 swept_values: list = None,
                 disable_progress_bar: bool = False,
                 maximum_integrator_steps: int = 20000
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
        self.dose = dose
        self.dose.set_ids(self.sbml_model)

        self.stop = stop
        self.steps = steps
        self.output_compartments = output_compartments
        self.output_parameters = output_parameters
        self.maximum_integrator_steps = maximum_integrator_steps


        if self.output_compartments and return_all_subcompartments:
            raise ValueError("Cannot specify output_compartments when return_all_subcompartments is True.")
        if self.output_parameters and return_all_parameters:
            raise ValueError("Cannot specify output_paramters when return_all_parameters is True")

        if return_all_subcompartments:
            self.output_compartments = self.get_compartments() + self.get_subcompartments()
        elif not self.output_compartments:
            self.output_compartments = self.get_compartments()

        if return_all_parameters:
            self.output_parameters = self.get_parameters(return_values=False)

        self.ids_to_return = self.get_return_ids()

        sbml_string = libsbml.writeSBMLToString(self.document)

        self.TACs_masks = sparse.csr_array(self.get_masks())

        self.start_times = np.array(self.dose.times)
        self.end_times = np.append(np.array(self.dose.times)[1:], self.stop)
        self.cycle_steps = np.rint((self.end_times - self.start_times) * self.steps / self.stop + 1).astype(int)

        all_times = []
        for start, end, n_steps in zip(self.start_times, self.end_times, self.cycle_steps):
            n_steps = int(n_steps)
            t_cycle = np.linspace(start, end, n_steps)
            all_times.append(t_cycle[:-1])

        self.time = np.concatenate(all_times)

        if swept_parameters:
            parameter_ids = []
            for parameter in swept_parameters:
                parameter_ids.append(get_parameter_id(
                    self.sbml_model, parameter))

            args = [(sbml_string, parameter_ids, values) for values in swept_values]
        else:
            args = [(sbml_string, None, None)]
            disable_progress_bar = True
        results = process_map(self.run_single_simulation, args, max_workers=os.cpu_count(
        ), chunksize=1, disable=disable_progress_bar)
        if self.output_parameters is not None:
            TACs, PARAMS = [np.stack(res, axis=0) for res in zip(*results)]
        else:
            TACs = np.stack(results, axis=0)

        if self.output_parameters is not None:
            return self.time, TACs, PARAMS
        else:
            return self.time, TACs

    def run_single_simulation(self, args):
        sbml_string, parameter_ids, swept_values = args
        document = libsbml.readSBMLFromString(sbml_string)
        sbml_model = document.getModel()
        if parameter_ids is not None:
            for i, pid in enumerate(parameter_ids):
                sbml_model.getParameter(pid).setValue(float(swept_values[i]))
            sbml_string = libsbml.writeSBMLToString(document)

        rr = roadrunner.RoadRunner(sbml_string)
        rr.integrator.maximum_num_steps = self.maximum_integrator_steps
        rr.timeCourseSelections = self.ids_to_return

        num_cycles = len(self.dose.times)
        all_time_segments = []
        all_results_segments = []
        for cycle in range(num_cycles):
            for index, id in enumerate(self.dose.ids):
                rr[f'{id}'] += list(self.dose.targets.values())[index][cycle]

            result_segment = rr.simulate(float(self.start_times[cycle]), float(self.end_times[cycle]), int(self.cycle_steps[cycle]))
            time_segment = np.linspace(self.start_times[cycle], self.end_times[cycle], self.cycle_steps[cycle])
            all_time_segments.append(time_segment[:-1])
            all_results_segments.append(result_segment[:-1,:])

        result = np.concatenate(all_results_segments)

        TAC = np.zeros((len(self.time), len(self.output_compartments)))
        TAC = result @ self.TACs_masks.T * self.NMOL2MBQ

        if self.output_parameters is not None:
            PARAMS = np.zeros((len(self.time), len(self.output_parameters)))
            PARAMS = np.stack(result[:,-len(self.output_parameters):])
            return TAC, PARAMS
        else:
            return TAC

    def get_return_ids(self):
        ids = []
        for region in self.output_compartments:
            if '+' in region:
                terms = [term.strip() for term in region.split("+")]
            else:
                terms = []
                terms.append(region)
            for term in terms:
                for compartment in self.sbml_model.getListOfCompartments():
                    if compartment.getName() == term:
                        compartment_id = compartment.getId()
                        for other_comp in self.sbml_model.getListOfCompartments():
                            if other_comp.getOutside() == compartment_id:
                                other_comp_id = other_comp.getId()
                                for species in self.sbml_model.getListOfSpecies():
                                    if species.getCompartment() == other_comp_id:
                                        ids.append(species.getId())
                        for species in self.sbml_model.getListOfSpecies():
                            if species.getCompartment() == compartment_id:
                                ids.append(species.getId())
                        break
        if self.output_parameters is not None:
            [ids.append(get_parameter_id(self.sbml_model, p)) for p in self.output_parameters]

        ids = list(dict.fromkeys(ids))
        return ids

    def save_sbml(self, path):
        """
        Saves current SBML model to file.

        Args:
            path (str): Path to save SBML file
        """
        libsbml.writeSBMLToFile(self.document, path)

    def get_compartments(self):
        compartment_list = []
        for compartment in self.sbml_model.getListOfCompartments():
            if not compartment.getOutside():
                compartment_list.append(compartment.getName())
        return compartment_list

    def get_subcompartments(self):
        subcompartment_list = []
        for compartment in self.sbml_model.getListOfCompartments():
            if compartment.getOutside():
                subcompartment_list.append(compartment.getName())
        return subcompartment_list

    def get_parameters(self, return_values=True):
        params = self.sbml_model.getListOfParameters()
        if return_values:
            return [(p.getName(), p.getValue()) for p in params]
        else:
            return [p.getName() for p in params]

    def get_masks(self):
        masks = np.zeros((len(self.output_compartments), len(self.ids_to_return)), dtype=bool)
        all_tags = self.get_tags('Hot')
        for idx, tags in enumerate(all_tags):
            masks[idx, :] = np.isin(self.ids_to_return, tags)
        return masks

    def get_tags(self, species_name):
        all_tags = []
        for region in self.output_compartments:
            tags = []
            if '+' in region:
                terms = [term.strip() for term in region.split("+")]
            else:
                terms = []
                terms.append(region)
            for term in terms:
                for compartment in self.sbml_model.getListOfCompartments():
                    if compartment.getName() == term:
                        compartment_id = compartment.getId()
                        for other_comp in self.sbml_model.getListOfCompartments():
                            if other_comp.getOutside() == compartment_id:
                                other_comp_id = other_comp.getId()
                                for species in self.sbml_model.getListOfSpecies():
                                    if species.getCompartment() == other_comp_id and species_name in species.getName():
                                        tags.append(species.getId())
                        for species in self.sbml_model.getListOfSpecies():
                            if species.getCompartment() == compartment_id and species_name in species.getName():
                                    tags.append(species.getId())
                        break
            all_tags.append(tags)
        return all_tags

    def create_jax_model(self, dose):
        for idx, tag in enumerate(dose.ids):
            species = self.sbml_model.getSpecies(tag)
            species.setInitialAmount(list(dose.targets.values())[idx][0])

        sbml_string = libsbml.writeSBMLToString(self.document)
        return convert_model_to_jax(sbml_string)

    # @eqx.filter_jit
    def compute_sensitivities(self, dose: Dose, t, output_compartments):
            dose.set_ids(self.sbml_model)
            rollout, name_list_y, _, name_list_c, y0, c = self.create_jax_model(dose)
            sig = inspect.signature(rollout)
            y0 = sig.parameters["y0"].default
            c = sig.parameters["c0"].default

            stepsize_controller = diffrax.PIDController(atol=1e-10, rtol=1e-3)
            def func(c, region_indices):
                ys, _, _, c_updated = rollout(
                    t1=0, ts=jnp.array([t-1, t]), deltaT=0.1, y0=y0, c0=c,
                    stepsize_controller=stepsize_controller,
                    max_steps=1_000_000,  # increased budget
                )
                tacs = jnp.zeros((len(region_indices), ys.shape[1]))
                for i, idxs in enumerate(region_indices):
                    tacs = tacs.at[i].set((ys.take(idxs, axis=0).sum(axis=0) * self.NMOL2MBQ))
                return tacs, (tacs, c_updated) #for returning gradient and values

            region_indices = []
            for region in output_compartments:
                if '+' in region:
                    terms = [term.strip() for term in region.split("+")]
                else:
                    terms = []
                    terms.append(region)
                temp_indices = jnp.array([])
                for term in terms:
                    temp_indices = jnp.concatenate([temp_indices, get_indices(term, name_list_y)], dtype=jnp.int32)
                region_indices.append(temp_indices)

            grads, (tacs, c_updated) = eqx.filter_jacrev(func, has_aux=True)(c, region_indices)
            grads = grads[:,-1,:]
            sens =  (grads*c_updated).T/tacs[:,-1]
            sens = sens / jnp.abs(jnp.max(sens))
            return pd.DataFrame(sens, index=name_list_c, columns=output_compartments)


def get_indices(region, compartment_list):
    matches = [i for i, c in enumerate(compartment_list) if c.startswith(f'Hot{region}')]
    if not matches:
        print(f"No compartments found for '{region}'")
        return jnp.zeros(1)
    return jnp.array(matches)

def set_parameter_values(sbml_model, parameter_dict_in):
    parameter_dict = parameter_dict_in.copy()
    for parameter in sbml_model.getListOfParameters():
        for parameter_name in list(parameter_dict.keys()):
            if parameter.getName() == parameter_name:
                parameter.setValue(float(parameter_dict.pop(parameter_name)))
    if parameter_dict:
        [print(f"Parameter {parameter_name} not found in the model.")
         for parameter_name in parameter_dict.keys()]

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

def get_species(species_name, result, sbml_model):
    NMOL2MBQ = get_parameter(sbml_model, 'lambdaPhys') / \
        60 * 6.022e23 / 10**9 / 10**6
    for species in sbml_model.getListOfSpecies():
        if species.getName() == species_name.split('.')[1]:
            return (result[f"{species.getId()}"] * NMOL2MBQ)
