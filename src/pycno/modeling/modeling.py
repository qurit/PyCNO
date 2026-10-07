import inspect
import multiprocessing as mp
import os
import pickle
from dataclasses import dataclass
from pathlib import Path
from typing import ClassVar

import diffrax
import equinox as eqx
import jax.numpy as jnp
import libsbml
import numpy as np
import pandas as pd
import roadrunner
from tqdm.auto import tqdm
from tqdm.contrib.concurrent import process_map

from pycno import Dose
from pycno.utils.jax_conversion import convert_model_to_jax

from ..exceptions import ModelError, SimulationError

mp.set_start_method("spawn", force=True)

# TODO: Warning for parameters altered at runtime
# TODO: Add cold species return
# TODO: Add all possible inputs for sweep
# TODO: Update watched attrs
# TODO: Fix FIM
# TODO: Add multidosing to jax
# TODO: Height/Weight Scaling


@dataclass
class SimulationResult:
    time: np.ndarray
    tacs: np.ndarray
    output_compartments: list
    parameters: np.ndarray | None = None


class Model:
    """
    Initializes SBML model.

    Args:
        model_name (str): Name of model in PyCNO or path to SBML file
        parameters (dict): Parameter input values
        compartment_volumes (dict): Compartment volumes in L
        hotamount (float): Hot ligand amount in nmol
        coldamount (float): Cold ligand amount in nmol
    """

    _watched_attrs: ClassVar[set[str]] = {
        "model_name",
        "parameters",
        "compartment_volumes",
    }

    def __init__(self, model_name, parameters=None, compartment_volumes=None):

        self._initializing = True

        self.model_name = model_name
        self.parameters = parameters
        self.compartment_volumes = compartment_volumes

        self.sbml_model: libsbml.Model
        self.document: libsbml.SBMLDocument
        self.output_parameters = None
        self.output_compartments = None

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
            except (TypeError, pickle.PicklingError):  # noqa: PERF203
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
            self.set_parameter_values(sbml_model, self.parameters)

        if self.compartment_volumes:
            self.set_compartment_sizes(sbml_model, self.compartment_volumes)

        self.sbml_model = sbml_model
        self.document = document

        self.NMOL2MBQ = (
            self.get_parameter(self.sbml_model, "lambdaPhys")
            / 60
            * 6.022e23
            / 10**9
            / 10**6
        )

        self.sbml_string = libsbml.writeSBMLToString(self.document)

        if "InitialAssignment" in self.sbml_string:
            self.sim_func = self.simulate_single_initial_assignments
            self.simulate_initial_assignments = True
            print(
                "Warning: InitialAssignments are used in model. Consider switching to AssignmentRules for faster simulation over paramter sweeps."
            )
        else:
            self.sim_func = self.simulate_single
            self.simulate_initial_assignments = False
            self.rr = roadrunner.RoadRunner(self.sbml_string)

    def simulate(
        self,
        dose: Dose,
        stop: int = 60,
        steps: int = 100,
        breakpoints: list | None = None,
        return_dose_times: bool = True,
        output_compartments: list | None = None,
        return_all_subcompartments: bool = False,
        output_parameters: list | None = None,
        return_all_parameters: bool = False,
        swept_parameters: list | None = None,
        swept_values: list | None = None,
        disable_progress_bar: bool = False,
        maximum_integrator_steps: int = 20000,
        parallel: bool = False,
        chunksize: int | None = None,
    ):
        """
        Simulates SBML model.

        Args:
            model_name (str): Name of model in PyCNO or path to SBML file
            observables (list): Regions to output
            stop (int): Simulation end time in minutes
            steps (int): Number of simulation steps
            breakpoints (list, optional): Extra segment-boundary times
                (minutes), in addition to any dose times, at which the local
                resolution of ``time`` changes (e.g. the fine/coarse
                boundary from :func:`build_hybrid_time`). Ignored if
                ``time`` is not given.
            parameters (dict): Parameter input values
            compartment_volumes (dict): Compartment volumes in L
            swept_parameters (list): Parameters to sweep over
            swept_values (list): Values to sweep over

        Returns:
            TACs[n_curves, n_steps, n_observables]: Time activity curves in units of MBq.

        """
        time = None  # TODO add explicit time grid
        if time is not None:
            self.time = np.asarray(time, dtype=np.float64)
            if self.time.ndim != 1 or self.time.size < 2:
                raise SimulationError("`time` must be a 1-D array with >= 2 points.")
            if self.time[0] != 0:
                raise SimulationError("`time` must start at 0.")
            if np.any(np.diff(self.time) <= 0):
                raise SimulationError("`time` must be strictly increasing.")
            self.stop = float(self.time[-1])
            self.steps = int(self.time.size)
        else:
            self.stop = float(stop)
            self.steps = int(steps)
            self.time = np.linspace(0, self.stop, self.steps)

        self._resolution_breakpoints = (
            np.asarray(breakpoints, dtype=np.float64)
            if (time is not None and breakpoints is not None and len(breakpoints) > 0)
            else np.array([], dtype=np.float64)
        )

        self.dose = dose
        self.dose.initialize_dose(self.stop, self.steps)
        self.dose.set_ids(self.sbml_model)

        self.num_cycles = len(self.dose.times)
        self.return_dose_times = return_dose_times

        if self.stop <= 0:
            raise SimulationError("Simulation stop time must be greater than 0.")

        (
            self.starts,
            self.stops,
            self.indv_steps,
            self.dose_sim_mask,
            self.time_indices,
        ) = self.build_simulation_segments()

        self.output_compartments = output_compartments
        self.output_parameters = output_parameters
        self.maximum_integrator_steps = maximum_integrator_steps

        if self.output_compartments and return_all_subcompartments:
            raise ValueError(
                "Cannot specify output_compartments when return_all_subcompartments is True."
            )
        if self.output_parameters and return_all_parameters:
            raise ValueError(
                "Cannot specify output_paramters when return_all_parameters is True"
            )

        if return_all_subcompartments:
            self.output_compartments = (
                self.get_compartments() + self.get_subcompartments()
            )
        elif not self.output_compartments:
            self.output_compartments = self.get_compartments()

        if return_all_parameters:
            self.output_parameters = self.get_parameters(return_values=False)

        self.ids_to_return = self.get_return_ids()
        self.TACs_masks = self.get_masks()

        if swept_parameters:
            if swept_values is None:
                raise ValueError(
                    "swept_values must be provided when swept_parameters is specified."
                )
            parameter_ids = [
                self.get_parameter_id(self.sbml_model, parameter)
                for parameter in swept_parameters
            ]

        if self.simulate_initial_assignments:
            if chunksize is None:
                chunksize = 1
            if swept_parameters:
                assert swept_values is not None
                args = [
                    (
                        self.sbml_string,
                        self.maximum_integrator_steps,
                        self.ids_to_return,
                        list(self.dose.targets.values()),
                        self.dose.ids,
                        self.starts,
                        self.stops,
                        self.indv_steps,
                        self.dose_sim_mask,
                        parameter_ids,
                        swept_vals,
                    )
                    for swept_vals in swept_values
                ]
            else:
                disable_progress_bar = True
                args = [
                    (
                        self.sbml_string,
                        self.maximum_integrator_steps,
                        self.ids_to_return,
                        list(self.dose.targets.values()),
                        self.dose.ids,
                        self.starts,
                        self.stops,
                        self.indv_steps,
                        self.dose_sim_mask,
                        None,
                        None,
                    )
                ]
        else:
            if chunksize is None:
                chunksize = 100
            self.rr.integrator.maximum_num_steps = self.maximum_integrator_steps
            self.rr.timeCourseSelections = self.ids_to_return
            if swept_parameters:
                assert swept_values is not None
                args = [
                    (
                        self.rr,
                        list(self.dose.targets.values()),
                        self.dose.ids,
                        self.starts,
                        self.stops,
                        self.indv_steps,
                        self.dose_sim_mask,
                        parameter_ids,
                        swept_vals,
                    )
                    for swept_vals in swept_values
                ]
            else:
                disable_progress_bar = True
                args = [
                    (
                        self.rr,
                        list(self.dose.targets.values()),
                        self.dose.ids,
                        self.starts,
                        self.stops,
                        self.indv_steps,
                        self.dose_sim_mask,
                        None,
                        None,
                    )
                ]

        if parallel:
            results = process_map(
                self.sim_func,
                args,
                max_workers=os.cpu_count(),
                chunksize=chunksize,
                disable=disable_progress_bar,
            )

        else:
            results = []
            if disable_progress_bar:
                for arg in args:
                    results.append(self.sim_func(arg))
            else:
                for arg in tqdm(args):
                    results.append(self.sim_func(arg))

        TAC = np.zeros((len(self.time), len(self.output_compartments)))
        TAC = np.stack(results) @ self.TACs_masks.T * self.NMOL2MBQ
        TAC = TAC[:, self.time_indices, :]

        if self.output_parameters is not None:
            PARAMS = np.zeros((len(self.time), len(self.output_parameters)))
            PARAMS = np.stack(results)[:, :, -len(self.output_parameters) :]
            PARAMS = PARAMS[:, self.time_indices, :]

        return SimulationResult(
            time=self.time,
            tacs=TAC,
            output_compartments=self.output_compartments,
            parameters=PARAMS if self.output_parameters is not None else None,
        )

    def get_return_ids(self):
        assert self.output_compartments is not None, (
            "output_compartments must be set before calling get_return_ids()"
        )
        ids = []
        for region in self.output_compartments:
            if "+" in region:
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
                                ids.extend(
                                    species.getId()
                                    for species in self.sbml_model.getListOfSpecies()
                                    if species.getCompartment() == other_comp_id
                                )
                        ids.extend(
                            species.getId()
                            for species in self.sbml_model.getListOfSpecies()
                            if species.getCompartment() == compartment_id
                        )
                        break
        if self.output_parameters is not None:
            [
                ids.append(self.get_parameter_id(self.sbml_model, p))
                for p in self.output_parameters
            ]

        ids = list(dict.fromkeys(ids))
        return ids

    def save_sbml(self, path):
        """
        Saves current SBML model to file.

        Args:
            path (str): Path to save SBML file
        """
        libsbml.writeSBMLToFile(self.document, str(path))

    def get_compartments(self):
        return [
            compartment.getName()
            for compartment in self.sbml_model.getListOfCompartments()
            if not compartment.getOutside()
        ]

    def get_subcompartments(self):
        return [
            compartment.getName()
            for compartment in self.sbml_model.getListOfCompartments()
            if compartment.getOutside()
        ]

    def get_parameters(self, return_values=True):
        params = self.sbml_model.getListOfParameters()
        if return_values:
            return [(p.getName(), p.getValue()) for p in params]
        else:
            return [p.getName() for p in params]

    def get_masks(self):
        assert self.output_compartments is not None, (
            "output_compartments must be set before calling get_masks()"
        )
        masks = np.zeros(
            (len(self.output_compartments), len(self.ids_to_return)), dtype=bool
        )
        all_tags = self.get_tags("Hot")
        for idx, tags in enumerate(all_tags):
            masks[idx, :] = np.isin(self.ids_to_return, tags)
        return masks

    def get_tags(self, species_name):
        assert self.output_compartments is not None, (
            "output_compartments must be set before calling get_tags()"
        )
        all_tags = []
        for region in self.output_compartments:
            tags = []
            if "+" in region:
                terms = [term.strip() for term in region.split("+")]
            else:
                terms = [region]
            for term in terms:
                for compartment in self.sbml_model.getListOfCompartments():
                    if compartment.getName() == term:
                        compartment_id = compartment.getId()
                        for other_comp in self.sbml_model.getListOfCompartments():
                            if other_comp.getOutside() == compartment_id:
                                other_comp_id = other_comp.getId()
                                tags.extend(
                                    species.getId()
                                    for species in self.sbml_model.getListOfSpecies()
                                    if species.getCompartment() == other_comp_id
                                    and species_name in species.getName()
                                )
                        tags.extend(
                            species.getId()
                            for species in self.sbml_model.getListOfSpecies()
                            if species.getCompartment() == compartment_id
                            and species_name in species.getName()
                        )
                        break
                else:
                    raise ModelError(f"Compartment {term} not found in model.")
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
                t1=0,
                ts=jnp.array([t - 1, t]),
                deltaT=0.1,
                y0=y0,
                c0=c,
                stepsize_controller=stepsize_controller,
                max_steps=1_000_000,  # increased budget
            )
            tacs = jnp.zeros((len(region_indices), ys.shape[1]))
            for i, idxs in enumerate(region_indices):
                tacs = tacs.at[i].set(ys.take(idxs, axis=0).sum(axis=0) * self.NMOL2MBQ)
            return tacs, (tacs, c_updated)  # for returning gradient and values

        region_indices = []
        for region in output_compartments:
            if "+" in region:
                terms = [term.strip() for term in region.split("+")]
            else:
                terms = []
                terms.append(region)
            temp_indices = jnp.array([])
            for term in terms:
                temp_indices = jnp.concatenate(
                    [temp_indices, self.get_indices(term, name_list_y)], dtype=jnp.int32
                )
            region_indices.append(temp_indices)

        grads, (tacs, c_updated) = eqx.filter_jacrev(func, has_aux=True)(
            c, region_indices
        )
        grads = grads[:, -1, :]
        sens = (grads * c_updated).T / tacs[:, -1]
        sens = sens / jnp.abs(jnp.max(sens))
        return pd.DataFrame(sens, index=name_list_c, columns=output_compartments)

    def build_simulation_segments(self):
        """
        Build contiguous integration segments that respect dose times, the
        requested stop time, and t=0, so doses always fall on a segment
        boundary and never get skipped over by the integrator.

        When ``self._resolution_breakpoints`` is non-empty (set by
        ``simulate(..., time=..., breakpoints=...)``), delegates to
        :meth:`_build_segments_with_breakpoints` instead, which additionally
        splits segments at those breakpoints and sizes each segment from how
        many points of the (possibly non-uniform) ``self.time`` actually fall
        inside it — this is what makes piecewise resolution (e.g. 1-second
        early, 1-minute late) actually take effect, since a single RoadRunner
        segment is always internally uniform.

        Returns:
            starts, ends: arrays of segment boundaries (len = n_segments)
            seg_steps: int array, number of integrator points per segment
            dose_at_start: bool array, True where a dose should be applied
                        at the *start* of that segment (post-dose state
                        is what gets returned for any output time that
                        coincides with this start)
        """
        resolution_breaks = getattr(
            self, "_resolution_breakpoints", np.array([], dtype=np.float64)
        )
        if resolution_breaks.size > 0:
            return self._build_segments_with_breakpoints(resolution_breaks)

        all_time = np.unique(np.round(np.concatenate((self.time, self.dose.times)), 4))

        if self.num_cycles > 1 or self.dose.times[0] != 0:
            starts = np.unique(
                np.concatenate(
                    (
                        [0],
                        self.dose.times,
                        all_time[np.where(np.isin(all_time, self.dose.times))[0] + 1],
                    )
                )
            ).astype(np.float64)

            stops = np.unique(
                np.concatenate(
                    (
                        [self.time[-1]],
                        self.dose.times[1:],
                        all_time[np.where(np.isin(all_time, self.dose.times))[0] + 1],
                    )
                )
            )

            if self.dose.times[0] == 0:
                starts = np.unique(
                    np.concatenate(
                        (
                            all_time[
                                np.where(np.isin(all_time, self.dose.times))[0] - 1
                            ][1:],
                            starts,
                        )
                    )
                ).astype(np.float64)

                stops = np.unique(
                    np.concatenate(
                        (
                            self.dose.times[1:],
                            all_time[
                                np.where(np.isin(all_time, self.dose.times))[0] - 1
                            ][1:],
                            stops,
                        )
                    )
                ).astype(np.float64)
            else:
                starts = np.unique(
                    np.concatenate(
                        (
                            all_time[
                                np.where(np.isin(all_time, self.dose.times))[0] - 1
                            ],
                            starts,
                        )
                    )
                ).astype(np.float64)

                stops = np.unique(
                    np.concatenate(
                        (
                            self.dose.times,
                            all_time[
                                np.where(np.isin(all_time, self.dose.times))[0] - 1
                            ],
                            stops,
                        )
                    )
                ).astype(np.float64)

            indv_steps = np.rint((stops - starts) / (self.stop / (self.steps - 1))) + 1
            indv_steps[np.isin(starts, self.dose.times)] = 2
            indv_steps[np.isin(stops, self.dose.times)] = 2
            indv_steps = indv_steps.astype(np.int64)

            dose_sim_mask = np.isin(starts, self.dose.times)

            time_indices = np.concatenate(
                [
                    np.r_[np.ones(n - 1), 0] if i < len(indv_steps) - 1 else np.ones(n)
                    for i, n in enumerate(indv_steps)
                ]
            ).astype(np.bool)

            if self.return_dose_times:
                self.time = all_time
            else:
                time_indices[time_indices] = np.isin(all_time, self.time)

        else:
            starts = np.array([0], dtype=np.float64)
            stops = np.array([self.time[-1]])
            indv_steps = np.array([self.steps], dtype=np.int64)
            dose_sim_mask = np.array([True])
            time_indices = np.ones_like(self.time).astype(np.bool)

        return starts, stops, indv_steps, dose_sim_mask, time_indices

    def _build_segments_with_breakpoints(self, resolution_breaks: np.ndarray):
        """
        Generalized segment builder used when an explicit (possibly
        non-uniform) ``self.time`` was supplied along with resolution
        ``breakpoints`` (see :meth:`simulate`).

        Segment boundaries are the union of ``{0, stop}``, any dose times,
        and ``resolution_breaks``. Each segment's integrator point count is
        read directly from how many points of ``self.time`` fall inside it
        (inclusive of both ends), so the caller's intended local resolution
        (e.g. 1 s early / 1 min late) is reproduced exactly, rather than
        being derived from one global average rate as in the dose-only path.

        Consecutive segments share their boundary sample (RoadRunner's
        ``simulate(start, stop, n)`` always includes both endpoints); the
        duplicate is dropped via ``time_indices``, same convention as
        :meth:`build_simulation_segments`.
        """
        t_end = float(self.time[-1])
        dose_times = np.asarray(self.dose.times, dtype=np.float64)

        boundaries = np.unique(
            np.round(np.concatenate(([0.0, t_end], dose_times, resolution_breaks)), 4)
        )
        boundaries = boundaries[(boundaries >= 0) & (boundaries <= t_end)]
        if boundaries.size < 2:
            boundaries = np.array([0.0, t_end])

        starts = boundaries[:-1]
        stops = boundaries[1:]

        t_round = np.round(self.time, 4)
        indv_steps = np.empty(len(starts), dtype=np.int64)
        for i, (s, e) in enumerate(zip(starts, stops)):
            n = int(np.sum((t_round >= s) & (t_round <= e)))
            indv_steps[i] = max(n, 2)

        dose_sim_mask = np.isin(starts, np.round(dose_times, 4))

        time_indices = np.concatenate(
            [
                np.r_[np.ones(n - 1), 0] if i < len(indv_steps) - 1 else np.ones(n)
                for i, n in enumerate(indv_steps)
            ]
        ).astype(np.bool_)

        return starts, stops, indv_steps, dose_sim_mask, time_indices

    @staticmethod
    def get_indices(region, compartment_list):
        matches = [
            i
            for i, c in enumerate(compartment_list)
            if c.startswith(f"{region}") and c.endswith(".Hot")
        ]
        if not matches:
            print(f"No compartments found for '{region}'")
            return jnp.zeros(1)
        return jnp.array(matches)

    @staticmethod
    def set_parameter_values(sbml_model, parameter_dict_in):
        parameter_dict = parameter_dict_in.copy()
        for parameter in sbml_model.getListOfParameters():
            for parameter_name in list(parameter_dict.keys()):
                if parameter.getName() == parameter_name:
                    parameter.setValue(float(parameter_dict.pop(parameter_name)))
        if parameter_dict:
            [
                print(f"Parameter {parameter_name} not found in the model.")
                for parameter_name in parameter_dict
            ]
            raise ModelError(
                "Some parameters were not found in the model. Check parameter names."
            )

    @staticmethod
    def set_compartment_sizes(sbml_model, compartment_dict_in):
        compartment_dict = compartment_dict_in.copy()
        for comp in sbml_model.getListOfCompartments():
            if comp.getName() in compartment_dict:
                comp.setSize(float(compartment_dict[comp.getName()]))
                compartment_dict.pop(comp.getName())
        if compartment_dict:
            for missing in compartment_dict:
                print(f"Compartment {missing} not found in the model.")
                raise ModelError(
                    "Some compartments were not found in the model. Check compartment names."
                )

    @staticmethod
    def get_parameter(sbml_model, parameter_name):
        for parameter in sbml_model.getListOfParameters():
            if parameter.getName() == parameter_name:
                return parameter.getValue()
        raise ModelError(f"Parameter {parameter_name} not found in the model.")

    @staticmethod
    def get_parameter_id(sbml_model, parameter_name):
        for parameter in sbml_model.getListOfParameters():
            if parameter.getName() == parameter_name:
                return parameter.getId()
        raise ModelError(f"Parameter {parameter_name} not found in the model.")

    @staticmethod
    def get_species(species_name, result, sbml_model):
        NMOL2MBQ = (
            Model.get_parameter(sbml_model, "lambdaPhys")
            / 60
            * 6.022e23
            / 10**9
            / 10**6
        )
        for species in sbml_model.getListOfSpecies():
            if species.getName() == species_name.split(".")[1]:
                return result[f"{species.getId()}"] * NMOL2MBQ

    @staticmethod
    def simulate_single(args):
        (
            rr,
            dose_values,
            dose_ids,
            starts,
            stops,
            indv_steps,
            dose_sim_mask,
            parameter_ids,
            swept_values,
        ) = args
        rr.reset()
        result_segments = []
        cycle = 0
        for sub_sim, dose_true in enumerate(dose_sim_mask):
            if dose_true:
                for index, id in enumerate(dose_ids):
                    rr[f"{id}"] += dose_values[index][cycle]
                cycle += 1
            if parameter_ids is not None:
                for i, pid in enumerate(parameter_ids):
                    rr.setValue(pid, float(swept_values[i]))
            result_segments.append(
                rr.simulate(starts[sub_sim], stops[sub_sim], int(indv_steps[sub_sim]))
            )
        return np.concatenate(result_segments)

    @staticmethod
    def simulate_single_initial_assignments(args):
        (
            sbml_string,
            maximum_integrator_steps,
            ids_to_return,
            dose_values,
            dose_ids,
            starts,
            stops,
            indv_steps,
            dose_sim_mask,
            parameter_ids,
            swept_values,
        ) = args
        rr = roadrunner.RoadRunner(sbml_string)
        rr.integrator.maximum_num_steps = maximum_integrator_steps
        rr.timeCourseSelections = ids_to_return
        result_segments = []
        cycle = 0
        for sub_sim, dose_true in enumerate(dose_sim_mask):
            if dose_true:
                for index, id in enumerate(dose_ids):
                    rr[f"{id}"] += dose_values[index][cycle]
                cycle += 1
            if parameter_ids is not None:
                for i, pid in enumerate(parameter_ids):
                    rr.setValue(pid, float(swept_values[i]))
            result_segments.append(
                rr.simulate(starts[sub_sim], stops[sub_sim], int(indv_steps[sub_sim]))
            )
        return np.concatenate(result_segments)
