from dataclasses import dataclass

import numpy as np

from ..exceptions import DoseError, SimulationError


@dataclass
class Dose:
    times: list | np.ndarray
    targets: dict

    def initialize_dose(self, stop, steps):
        self.times = np.round(np.array(self.times, dtype=float), 4)
        if np.any(self.times < 0):
            raise DoseError("Dose times must be non-negative.")
        if np.any(self.times > stop):
            raise DoseError("Dose time exceeds simulation stop time.")
        if not np.all(np.diff(self.times) >= 0):
            raise DoseError("Dose times must be sorted (non-decreasing).")
        if np.any(np.diff(self.times) < np.round(3.1 * stop / steps, 2)):
            raise SimulationError(
                "Not enough simulation steps for dose times. Increase steps."
            )
        if self.times[0] < 2.1 * stop / steps and self.times[0] != 0:
            raise SimulationError(
                "First dose time is too close to t=0. Increase steps."
            )
        if self.times[-1] > stop - 2.1 * stop / steps:
            raise SimulationError(
                "Last dose time is too close to simulation stop time. Increase steps or stop."
            )

        if len(self.times) > 1:
            for target in self.targets.items():
                if len(target[1]) == 1:
                    self.targets[target[0]] = [target[1][0]] * len(self.times)
                elif len(target[1]) != len(self.times):
                    raise DoseError(
                        "Dose target values must either match number of dose times or be a single value to be used for all times."
                    )

    def set_ids(self, sbml_model):
        self.ids = []
        for species_name in list(self.targets.keys()):
            for compartment in sbml_model.getListOfCompartments():
                if compartment.getName() == species_name.split(".")[0]:
                    break
                else:
                    compartment = None
            if compartment is None:
                comp = species_name.split(".")[0]
                raise DoseError(
                    f"{comp} compartment not found in model. Change dose target."
                )
            for species in sbml_model.getListOfSpecies():
                if (
                    species.getCompartment() == compartment.getId()
                    and species.getName() == species_name.split(".")[1]
                ):
                    self.ids.append(species.getId())
                    break
