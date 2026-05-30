import sbmltoodejax
import inspect
import tempfile
import importlib
import sys
import os
import xml.etree.ElementTree as ET
from typing import Any, Tuple, List

def convert_model_to_jax(model_string):
    model_name = 'model_temp'
    model = sbmltoodejax.parse.ParseSBMLFile(model_string)

    with tempfile.TemporaryDirectory() as tmpdir:
        module_file = os.path.join(tmpdir, f'{model_name}.py')
        sbmltoodejax.modulegeneration.GenerateModel(model, module_file)
        sys.path.insert(0, tmpdir)
        try:
            model_file = importlib.import_module(model_name)
            rollout = model_file.ModelRollout
        finally:
            sys.path.pop(0)
    name_list_y, name_list_w, name_list_c = get_rollout_names(model_string, rollout())

    sig = inspect.signature(rollout.__call__)
    y0 = sig.parameters["y0"].default
    c  = sig.parameters["c0"].default

    return rollout(), name_list_y, name_list_w, name_list_c, y0, c

def get_rollout_names(model_string: str, rollout: Any) -> Tuple[List[str], List[str], List[str]]:
    """
    Parse an SBML file and map rollout IDs to human-readable names.

    Returns:
        (name_list_y, name_list_w, name_list_c)
    """
    root = ET.fromstring(model_string)
    ns = {"sbml": root.tag.split('}')[0].strip('{')}

    # Separate dictionaries for compartments, parameters, species
    id_to_name_compartments = {
        c.get("id"): c.get("name") or c.get("id")
        for c in root.findall(".//sbml:compartment", ns)
    }
    id_to_name_parameters = {
        p.get("id"): p.get("name") or p.get("id")
        for p in root.findall(".//sbml:parameter", ns)
    }
    id_to_name_species = {}
    for sp in root.findall(".//sbml:species", ns):
        sid = sp.get("id")
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
    name_list_y = [id_to_name[k] for k in rollout.y_indexes.keys() if k in id_to_name]
    name_list_w = [id_to_name[k] for k in rollout.w_indexes.keys() if k in id_to_name]
    name_list_c = [id_to_name[k] for k in rollout.c_indexes.keys() if k in id_to_name]

    return name_list_y, name_list_w, name_list_c
