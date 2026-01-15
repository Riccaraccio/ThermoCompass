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
