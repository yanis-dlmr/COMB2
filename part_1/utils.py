def species_name_to_label(name: str) -> str:
    """
    Convert a species name to a LaTeX label.
    For example, 'CH4' becomes 'CH$_4$'.
    """
    label = ''
    for c in name:
        if c.isdigit():
            label += f'$_{c}$'
        else:
            label += c
    return label

def generate_X (species: list) -> str:
    """
    Generate a string of the form "CH4:1, N2:7.52, O2:2" from a list of dictionnaries
    """
    X: str = ""
    for i in range(len(species)):
        X += species[i]['name'] + ":" + str(species[i]['X'])
        if i != len(species)-1:
            X += ", "
    return X