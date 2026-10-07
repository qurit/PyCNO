import importlib
import inspect
import sys
import tempfile
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Any

from sbmltoodejax.modulegeneration import GenerateModel
from sbmltoodejax.parse import ParseSBMLFile


def convert_model_to_jax(model_string):
    model_name = "model_temp"
    model = ParseSBMLFile(model_string)

    with tempfile.TemporaryDirectory() as tmpdir:
        module_file = Path(tmpdir) / f"{model_name}.py"
        GenerateModel(model, module_file)
        sys.path.insert(0, tmpdir)
        try:
            model_file = importlib.import_module(model_name)
            rollout = model_file.ModelRollout
        finally:
            sys.path.pop(0)
    name_list_y, name_list_w, name_list_c = get_rollout_names(model_string, rollout())

    sig = inspect.signature(rollout.__call__)
    y0 = sig.parameters["y0"].default
    c = sig.parameters["c0"].default

    return rollout(), name_list_y, name_list_w, name_list_c, y0, c


def get_rollout_names(
    model_string: str, rollout: Any
) -> tuple[list[str], list[str], list[str]]:
    """Parse an SBML file and map rollout IDs to human-readable names.

    Returns:
        (name_list_y, name_list_w, name_list_c)

    """
    root = ET.fromstring(model_string)
    ns = {"sbml": root.tag.split("}")[0].strip("{")}

    # Separate dictionaries for compartments, parameters, species
    id_to_name_compartments: dict[str, str] = {}
    for c in root.findall(".//sbml:compartment", ns):
        cid = c.get("id")
        if cid is None:
            raise ValueError("SBML compartment is missing an ID")
        id_to_name_compartments[cid] = c.get("name") or cid

    id_to_name_parameters: dict[str, str] = {}
    for p in root.findall(".//sbml:parameter", ns):
        pid = p.get("id")
        if pid is None:
            raise ValueError("SBML parameter is missing an ID")
        id_to_name_parameters[pid] = p.get("name") or pid

    id_to_name_species = {}
    for sp in root.findall(".//sbml:species", ns):
        sid = sp.get("id")
        if sid is None:
            raise ValueError("SBML species is missing an ID")
        sname = sp.get("name") or sid
        comp_id = sp.get("compartment")
        if comp_id in id_to_name_compartments:
            sname = f"{id_to_name_compartments[comp_id]}.{sname}"
        id_to_name_species[sid] = sname

    # Merge dictionaries in the same order as second snippet
    id_to_name = {}
    id_to_name.update(id_to_name_compartments)
    id_to_name.update(id_to_name_parameters)
    id_to_name.update(id_to_name_species)

    # Preserve rollout key order
    name_list_y = [id_to_name[k] for k in rollout.y_indexes if k in id_to_name]
    name_list_w = [id_to_name[k] for k in rollout.w_indexes if k in id_to_name]
    name_list_c = [id_to_name[k] for k in rollout.c_indexes if k in id_to_name]

    return name_list_y, name_list_w, name_list_c
