import libsbml
import numpy as np
import pytest

from pycno import Dose

from ..exceptions import DoseError, SimulationError


def test_initialize_dose_valid():
    dose = Dose(
        times=[0, 10, 20],
        targets={
            "Blood.Hot": [10, 20, 30],
            "Blood.Cold": [100, 200, 300],
        },
    )

    dose.initialize_dose(stop=60, steps=100)

    np.testing.assert_array_equal(dose.times, [0, 10, 20])
    assert dose.targets["Blood.Hot"] == [10, 20, 30]
    assert dose.targets["Blood.Cold"] == [100, 200, 300]


def test_initialize_dose_converts_and_rounds_times():
    dose = Dose(
        times=[2.123456, 10.987654],
        targets={"Blood.Hot": [10, 20]},
    )

    dose.initialize_dose(stop=60, steps=100)

    np.testing.assert_array_equal(dose.times, [2.1235, 10.9877])


def test_initialize_dose_rejects_negative_time():
    dose = Dose(
        times=[-1],
        targets={"Blood.Hot": [10]},
    )

    with pytest.raises(DoseError, match="non-negative"):
        dose.initialize_dose(stop=60, steps=100)


def test_initialize_dose_rejects_time_after_stop():
    dose = Dose(
        times=[61],
        targets={"Blood.Hot": [10]},
    )

    with pytest.raises(DoseError, match="exceeds simulation stop"):
        dose.initialize_dose(stop=60, steps=100)


def test_initialize_dose_rejects_unsorted_times():
    dose = Dose(
        times=[20, 10],
        targets={"Blood.Hot": [10, 20]},
    )

    with pytest.raises(DoseError, match="sorted"):
        dose.initialize_dose(stop=60, steps=100)


def test_initialize_dose_rejects_doses_too_close_together():
    # 3.1 * 60 / 100 = 1.86 min
    dose = Dose(
        times=[10, 11],
        targets={"Blood.Hot": [10, 20]},
    )

    with pytest.raises(SimulationError, match="Not enough simulation steps"):
        dose.initialize_dose(stop=60, steps=100)


def test_initialize_dose_rejects_first_dose_too_close_to_zero():
    # 2.1 * 60 / 100 = 1.26 min
    dose = Dose(
        times=[1],
        targets={"Blood.Hot": [10]},
    )

    with pytest.raises(SimulationError, match="too close to t=0"):
        dose.initialize_dose(stop=60, steps=100)


def test_initialize_dose_allows_dose_at_zero():
    dose = Dose(
        times=[0],
        targets={"Blood.Hot": [10]},
    )

    dose.initialize_dose(stop=60, steps=100)

    np.testing.assert_array_equal(dose.times, [0])


def test_initialize_dose_rejects_last_dose_too_close_to_stop():
    # 60 - 2.1 * 60 / 100 = 58.74 min
    dose = Dose(
        times=[59],
        targets={"Blood.Hot": [10]},
    )

    with pytest.raises(
        SimulationError,
        match="too close to simulation stop time",
    ):
        dose.initialize_dose(stop=60, steps=100)


def test_initialize_dose_broadcasts_single_target_value():
    dose = Dose(
        times=[0, 10, 20],
        targets={"Blood.Hot": [10]},
    )

    dose.initialize_dose(stop=60, steps=100)

    assert dose.targets["Blood.Hot"] == [10, 10, 10]


def test_initialize_dose_rejects_wrong_target_length():
    dose = Dose(
        times=[0, 10, 20],
        targets={"Blood.Hot": [10, 20]},
    )

    with pytest.raises(
        DoseError,
        match="must either match number of dose times",
    ):
        dose.initialize_dose(stop=60, steps=100)


def test_single_dose_does_not_require_target_broadcasting():
    dose = Dose(
        times=[0],
        targets={"Blood.Hot": [10]},
    )

    dose.initialize_dose(stop=60, steps=100)

    assert dose.targets["Blood.Hot"] == [10]


def test_set_ids_maps_target_names_to_sbml_ids():
    doc = libsbml.readSBML("src/pycno/models/PSMA.sbml")
    sbml_model = doc.getModel()

    dose = Dose(
        times=[0],
        targets={"Art.Hot": [10]},
    )

    dose.set_ids(sbml_model)

    species = sbml_model.getElementBySId("mwef7a6362_1ff0_46f5_a04d_3f766fe34a1f")

    assert dose.ids == [species.getId()]


def test_set_ids_maps_multiple_targets_in_order():
    doc = libsbml.readSBML("src/pycno/models/PSMA.sbml")
    sbml_model = doc.getModel()

    dose = Dose(
        times=[0],
        targets={
            "Art.Hot": [10],
            "Art.Cold": [20],
            "Tumor1.R": [30],
        },
    )

    dose.set_ids(sbml_model)

    assert dose.ids == [
        "mwef7a6362_1ff0_46f5_a04d_3f766fe34a1f",
        "mw4e9d0b6d_b2b7_46ca_9e70_41a875cdf5dd",
        "mwf964af22_c053_4c6a_8f95_368a588ea3ec",
    ]


def test_set_ids_rejects_unknown_compartment():
    doc = libsbml.readSBML("src/pycno/models/PSMA.sbml")
    sbml_model = doc.getModel()

    dose = Dose(
        times=[0],
        targets={"NotACompartment.Hot": [10]},
    )

    with pytest.raises(DoseError, match="compartment not found"):
        dose.set_ids(sbml_model)


def test_set_ids_does_not_add_id_for_unknown_species():
    doc = libsbml.readSBML("src/pycno/models/PSMA.sbml")
    sbml_model = doc.getModel()

    dose = Dose(
        times=[0],
        targets={"Art.NotASpecies": [10]},
    )

    dose.set_ids(sbml_model)

    assert dose.ids == []
