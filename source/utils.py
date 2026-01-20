def replace_line(filepath: str, search_string: str, new_line: str):
    """Replace a line in a file that contains a specific search string."""
    with open(filepath) as f:
        lines = f.readlines()

    with open(filepath, "w") as f:
        for i, line in enumerate(lines):
            if search_string in line:
                lines[i] = new_line + "\n"
        f.writelines(lines)


def from_comp_dict_to_str(comp_dict: dict) -> str:
    """Convert a composition dictionary to a formatted string."""
    return " ".join([f"{k} {v}" for k, v in comp_dict.items()])


def convert_species_name(species_name: str, model) -> str:
    """Handles the conversion of species like MOIST to the appropriate name in the older model."""
    if species_name == "MOIST" and (model == "Polimi_1805" or model == "Polimi_1402"):
        return "ACQUA"

    if species_name in ["XYHW", "GMSW", "XYLW"] and model == "Polimi_1402":
        return "HCE"

    # If no conversion is needed, return the original name
    return species_name


def convert_comp_dict(comp_dict: dict, model: str) -> dict:
    """Convert the keys of a composition dictionary based on the model."""
    # if TANN or TGL are present, they need to be removed and composition re-normalized
    renormalize = False
    if ("TANN" in comp_dict or "TGL" in comp_dict) and model == "Polimi_1402":
        renormalize = True
        if "TANN" in comp_dict:
            comp_dict.pop("TANN")
        if "TGL" in comp_dict:
            comp_dict.pop("TGL")

    # Special handling for Polimi_dummy model
    if (model == "Polimi_dummy"):
        biomass_fraction = 0.0
        for species in list(comp_dict.keys()):
            # exclude CHAR, ASH, MOIST from biomass fraction
            if species not in ["CHAR", "ASH", "MOIST"]:
                biomass_fraction += comp_dict.pop(species)
        return {"BIOMASS": biomass_fraction, **comp_dict}

    converted_dict = {}
    for species, fraction in comp_dict.items():
        converted_name = convert_species_name(species, model)
        converted_dict[converted_name] = fraction

    if renormalize:
        total_fraction = sum(converted_dict.values())
        for species in converted_dict:
            converted_dict[species] /= total_fraction

    return converted_dict
