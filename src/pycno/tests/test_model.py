import numpy as np
import pytest

from pycno.modeling.modeling import Dose, Model, ModelError, SimulationResult


def test_model_loads_by_name():
    model = Model("PSMA")

    assert model.sbml_model is not None
    assert model.document is not None
    assert model.sbml_string is not None
    assert model.NMOL2MBQ > 0


def test_model_invalid_name():
    with pytest.raises(FileNotFoundError, match=r"Model .* not found"):
        Model("NotARealModel")


def test_get_compartments():
    model = Model("PSMA")
    assert model.sbml_model is not None

    compartments = model.get_compartments()

    assert isinstance(compartments, list)
    assert "Art" in compartments
    assert "Tumor1" in compartments


def test_get_subcompartments():
    model = Model("PSMA")
    assert model.sbml_model is not None

    compartments = model.get_subcompartments()

    assert isinstance(compartments, list)
    assert len(compartments) > 0


def test_get_parameters():
    model = Model("PSMA")
    assert model.sbml_model is not None

    parameters = model.get_parameters()

    assert isinstance(parameters, list)
    assert len(parameters) > 0
    assert all(len(parameter) == 2 for parameter in parameters)


def test_get_parameter():
    model = Model("PSMA")
    assert model.sbml_model is not None

    value = model.get_parameter(model.sbml_model, "lambdaPhys")

    assert isinstance(value, float)
    assert value > 0


def test_get_parameter_invalid_name():
    model = Model("PSMA")
    assert model.sbml_model is not None

    with pytest.raises(ModelError, match="Parameter NotAParameter not found"):
        model.get_parameter(model.sbml_model, "NotAParameter")


def test_get_parameter_id():
    model = Model("PSMA")
    assert model.sbml_model is not None

    parameter_id = model.get_parameter_id(model.sbml_model, "lambdaPhys")

    assert isinstance(parameter_id, str)
    assert parameter_id != ""


def test_get_parameter_id_invalid_name():
    model = Model("PSMA")
    assert model.sbml_model is not None

    with pytest.raises(ModelError, match="Parameter NotAParameter not found"):
        model.get_parameter_id(model.sbml_model, "NotAParameter")


def test_set_parameter_values():
    model = Model("PSMA")
    assert model.sbml_model is not None

    original = model.get_parameter(model.sbml_model, "lambdaPhys")

    model.set_parameter_values(
        model.sbml_model,
        {"lambdaPhys": original * 2},
    )

    updated = model.get_parameter(model.sbml_model, "lambdaPhys")

    assert updated == pytest.approx(original * 2)


def test_get_parameters_without_values():
    model = Model("PSMA")
    assert model.sbml_model is not None

    parameters = model.get_parameters(return_values=False)

    assert isinstance(parameters, list)
    assert len(parameters) > 0
    assert "lambdaPhys" in parameters


def test_set_compartment_sizes():
    model = Model("PSMA")
    assert model.sbml_model is not None

    compartment = next(
        comp
        for comp in model.sbml_model.getListOfCompartments()
        if comp.getName() == "Art"
    )
    original = compartment.getSize()

    model.set_compartment_sizes(
        model.sbml_model,
        {"Art": original * 2},
    )

    assert compartment.getSize() == pytest.approx(original * 2)


def test_set_compartment_sizes_invalid_name():
    model = Model("PSMA")
    assert model.sbml_model is not None

    with pytest.raises(ModelError, match="Some compartments were not found"):
        model.set_compartment_sizes(
            model.sbml_model,
            {"NotACompartment": 1.0},
        )


def test_model_parameter_override():
    model = Model(
        "PSMA",
        parameters={"lambdaPhys": 2.0},
    )
    assert model.sbml_model is not None

    assert model.get_parameter(model.sbml_model, "lambdaPhys") == pytest.approx(2.0)


def test_model_compartment_volume_override():
    model = Model(
        "PSMA",
        compartment_volumes={"Art": 2.0},
    )
    assert model.sbml_model is not None

    compartment = next(
        comp
        for comp in model.sbml_model.getListOfCompartments()
        if comp.getName() == "Art"
    )

    assert compartment.getSize() == pytest.approx(2.0)


def test_get_indices():
    model = Model("PSMA")
    assert model.sbml_model is not None

    compartment_list = [
        "Art.Hot",
        "Art.Cold",
        "Tumor1.Hot",
        "Tumor1.Cold",
    ]

    indices = model.get_indices("Art", compartment_list)

    assert np.array_equal(indices, np.array([0]))


def test_get_indices_no_match(capsys):
    model = Model("PSMA")
    assert model.sbml_model is not None

    compartment_list = [
        "Art.Hot",
        "Art.Cold",
        "Tumor1.Hot",
        "Tumor1.Cold",
    ]

    indices = model.get_indices("NotACompartment", compartment_list)

    assert np.array_equal(indices, np.zeros(1))
    assert "No compartments found" in capsys.readouterr().out


def test_get_return_ids():
    model = Model("PSMA")
    assert model.sbml_model is not None

    model.output_compartments = ["Art"]

    ids = model.get_return_ids()

    assert isinstance(ids, list)
    assert len(ids) > 0

    species_ids = {species.getId() for species in model.sbml_model.getListOfSpecies()}

    assert set(ids).issubset(species_ids)


def test_get_return_ids_combined_compartments():
    model = Model("PSMA")
    assert model.sbml_model is not None

    model.output_compartments = ["Art + Tumor1"]

    ids = model.get_return_ids()

    assert isinstance(ids, list)
    assert len(ids) > 0

    names = {
        species.getName()
        for species in model.sbml_model.getListOfSpecies()
        if species.getId() in ids
    }

    assert "Hot" in names


def test_get_tags():
    model = Model("PSMA")
    assert model.sbml_model is not None

    model.output_compartments = ["Art"]

    tags = model.get_tags("Hot")

    assert isinstance(tags, list)
    assert len(tags) == 1
    assert len(tags[0]) > 0

    species_ids = {species.getId() for species in model.sbml_model.getListOfSpecies()}

    assert set(tags[0]).issubset(species_ids)


def test_get_masks():
    model = Model("PSMA")
    assert model.sbml_model is not None

    model.output_compartments = ["Art", "Tumor1"]
    model.ids_to_return = model.get_return_ids()

    masks = model.get_masks()

    assert isinstance(masks, np.ndarray)
    assert masks.dtype == bool
    assert masks.shape == (
        len(model.output_compartments),
        len(model.ids_to_return),
    )


def test_get_species():
    model = Model("PSMA")
    assert model.sbml_model is not None

    species = next(
        s for s in model.sbml_model.getListOfSpecies() if s.getName() == "Hot"
    )

    result = {species.getId(): 10.0}

    value = model.get_species(
        "Art.Hot",
        result,
        model.sbml_model,
    )

    expected = 10.0 * model.NMOL2MBQ

    assert value == pytest.approx(expected)


def test_get_species_invalid_name():
    model = Model("PSMA")
    assert model.sbml_model is not None

    with pytest.raises((IndexError, KeyError, ModelError)):
        model.get_species(
            "NotARealSpecies",
            {},
            model.sbml_model,
        )


def test_model_getstate():
    model = Model("PSMA")
    assert model.sbml_model is not None

    state = model.__getstate__()

    assert isinstance(state, dict)
    assert "model_name" in state
    assert "parameters" in state
    assert "compartment_volumes" in state


def test_model_setstate():
    model = Model("PSMA")
    assert model.sbml_model is not None
    state = model.__getstate__()

    restored = Model.__new__(Model)
    restored.__setstate__(state)

    assert restored.model_name == model.model_name
    assert restored.parameters == model.parameters
    assert restored.compartment_volumes == model.compartment_volumes


def test_model_pickle():
    import pickle

    model = Model("PSMA")
    assert model.sbml_model is not None

    serialized = pickle.dumps(model)
    restored = pickle.loads(serialized)

    assert restored.model_name == model.model_name
    assert restored.parameters == model.parameters
    assert restored.compartment_volumes == model.compartment_volumes


def test_save_sbml(tmp_path):
    model = Model("PSMA")
    assert model.sbml_model is not None

    output_path = tmp_path / "test_model.sbml"
    model.save_sbml(output_path)

    assert output_path.exists()
    assert output_path.stat().st_size > 0


def test_save_sbml_is_valid(tmp_path):
    import libsbml

    model = Model("PSMA")
    assert model.sbml_model is not None

    output_path = tmp_path / "test_model.sbml"
    model.save_sbml(output_path)

    document = libsbml.readSBML(str(output_path))

    assert document.getModel() is not None
    assert document.getModel().getNumCompartments() > 0
    assert document.getModel().getNumSpecies() > 0


def test_create_jax_model():
    model = Model("PSMA")
    assert model.sbml_model is not None

    dose = Dose(
        times=[1.0],
        targets={
            "Art.Hot": [10.0],
        },
    )

    dose.initialize_dose(stop=10.0, steps=100)
    dose.set_ids(model.sbml_model)

    jax_model = model.create_jax_model(dose)

    assert jax_model is not None

    species = model.sbml_model.getSpecies(dose.ids[0])
    assert species.getInitialAmount() == pytest.approx(10.0)


def test_build_simulation_segments_no_doses():
    model = Model("PSMA")
    assert model.sbml_model is not None

    model.time = np.linspace(0, 10, 11)
    model.stop = 10.0
    model.steps = 11
    model.dose = Dose(times=[0], targets={"Art.Hot": [1.0]})
    model.num_cycles = 1
    model.return_dose_times = True
    model._resolution_breakpoints = np.array([])

    starts, stops, steps, dose_mask, time_indices = model.build_simulation_segments()

    assert np.array_equal(starts, [0.0])
    assert np.array_equal(stops, [10.0])
    assert np.array_equal(steps, [11])
    assert np.array_equal(dose_mask, [True])
    assert np.all(time_indices)


def test_build_simulation_segments_with_dose():
    model = Model("PSMA")
    assert model.sbml_model is not None

    model.time = np.linspace(0, 10, 11)
    model.stop = 10.0
    model.steps = 11
    model.dose = Dose(times=[2.0], targets={"Art.Hot": [1.0]})
    model.num_cycles = 1
    model.return_dose_times = True
    model._resolution_breakpoints = np.array([])

    starts, stops, steps, dose_mask, time_indices = model.build_simulation_segments()

    assert np.array_equal(starts, [0.0, 1.0, 2.0, 3.0])
    assert np.array_equal(stops, [1.0, 2.0, 3.0, 10.0])
    assert np.array_equal(steps, [2, 2, 2, 8])
    assert np.array_equal(dose_mask, [False, False, True, False])
    assert isinstance(time_indices, np.ndarray)


def test_build_simulation_segments_multiple_doses():
    model = Model("PSMA")
    assert model.sbml_model is not None

    model.time = np.linspace(0, 10, 11)
    model.stop = 10.0
    model.steps = 11
    model.dose = Dose(times=[2.0, 6.0], targets={"Art.Hot": [1.0, 1.0]})
    model.num_cycles = 2
    model.return_dose_times = True
    model._resolution_breakpoints = np.array([])

    starts, stops, steps, dose_mask, time_indices = model.build_simulation_segments()

    assert np.array_equal(starts, [0.0, 1.0, 2.0, 3.0, 5.0, 6.0, 7.0])
    assert np.array_equal(stops, [1.0, 2.0, 3.0, 5.0, 6.0, 7.0, 10.0])
    assert np.array_equal(steps, [2, 2, 2, 3, 2, 2, 4])
    assert np.array_equal(dose_mask, [False, False, True, False, False, True, False])
    assert isinstance(time_indices, np.ndarray)


def test_build_simulation_segments_with_breakpoints():
    model = Model("PSMA")
    assert model.sbml_model is not None

    model.time = np.array([0, 1, 2, 3, 4, 5, 10], dtype=float)
    model.stop = 10.0
    model.steps = len(model.time)
    model.dose = Dose(times=[5.0], targets={"Art.Hot": [1.0]})
    model.num_cycles = 1
    model.return_dose_times = True
    model._resolution_breakpoints = np.array([3.0])

    starts, stops, steps, dose_mask, time_indices = model.build_simulation_segments()

    assert np.array_equal(starts, [0.0, 3.0, 5.0])
    assert np.array_equal(stops, [3.0, 5.0, 10.0])
    assert np.array_equal(steps, [4, 3, 2])
    assert np.array_equal(dose_mask, [False, False, True])
    assert isinstance(time_indices, np.ndarray)


def test_simulate_basic():
    model = Model("PSMA")
    assert model.sbml_model is not None

    dose = Dose(
        times=[2.0],
        targets={"Art.Hot": [10.0]},
    )

    result = model.simulate(
        dose=dose,
        stop=10,
        steps=11,
        disable_progress_bar=True,
    )

    assert isinstance(result, SimulationResult)
    assert result.time.shape == (11,)
    assert result.tacs.shape[1] == 11
    assert result.tacs.shape[2] == len(result.output_compartments)


def test_simulate_dose_affects_activity():
    model = Model("PSMA")
    assert model.sbml_model is not None

    dose = Dose(
        times=[2.0],
        targets={"Art.Hot": [10.0]},
    )

    result = model.simulate(
        dose=dose,
        stop=10,
        steps=11,
        output_compartments=["Art"],
        disable_progress_bar=True,
    )

    assert np.all(result.tacs >= 0)
    assert np.any(result.tacs[0, 2:] > 0)


def test_simulate_output_parameters():
    model = Model("PSMA")
    assert model.sbml_model is not None

    dose = Dose(
        times=[2.0],
        targets={"Art.Hot": [10.0]},
    )

    result = model.simulate(
        dose=dose,
        stop=10,
        steps=11,
        output_compartments=["Art"],
        output_parameters=["lambdaPhys"],
        disable_progress_bar=True,
    )

    assert result.parameters is not None
    assert result.parameters.shape == (1, 11, 1)
    assert np.all(result.parameters > 0)
