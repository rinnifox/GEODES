import pandas as pd

from geodes import utils


def _calc_charge_clamp_dist(chain, charge_clamps):
    """Calculation of distance between charge clamp residues."""

    # extracting vectors of coordinates for every residue in charge clamp
    clamp_vectors = {}

    for res in chain:
        if res.id[1] == charge_clamps[0] or res.id[1] == charge_clamps[1] or res.id[1] == charge_clamps[2]:
            for atom in res:
                if atom.get_name() == 'CA':
                    clamp_vectors[res.id[1]] = atom.get_vector()
    table = {}
    for elem in range(len(clamp_vectors)):
        table[f'{charge_clamps[elem]}-{charge_clamps[elem-1]}'] = (
            clamp_vectors[charge_clamps[elem]]
            - clamp_vectors[charge_clamps[elem-1]]
        )

    # calculation of distance between charge clamp residues
    dist = {}

    for line in table:
        dist[line] = (table[line][0]**2 + table[line][1]**2 + table[line][2]**2) ** 0.5

    return dist


def calc_charge_clamp_dist(pdb_file, charge_clamps):
    """
    Calculation of distance between charge clamp residues.

    Parameters
    ----------
    pdb_file: str
        Filename of .pdb file used for calculation.
    charge_clamps: list of ints
        Charge clamp residues list.

    Returns
    -------
    dict of pairwise distances between charge clamp residues.

    """

    _, _, _, chain, _ = utils.get_model_and_structure(pdb_file)
    if not isinstance(charge_clamps, list):
        if charge_clamps is None:
            raise ValueError("Charge clamp residues list is None!")
        else:
            raise ValueError(f"Unexpected type for Charge clamp: {type(charge_clamps)}")
    return _calc_charge_clamp_dist(chain, charge_clamps)


def charge_clamp_dist_to_pandas(pdb_file, clamp_resid, protein_name=None, **kwargs):
    """
    Putting distance between charge clamp residues in pandas dataframe.

    Parameters
    ----------
    pdb_file: str
        Filename of .pdb file used for calculation.
    clamp_resid: list of ints
        Charge clamp residues list.
    protein_name: str, default=None
        Protein name to be added to the resulting dataframe.

    Returns
    -------
    pandas.DataFrame with calculated descriptor.

    """
    cols_cl_dist = ['prot_name'] + [f'Dist clamp{elem}' for elem in range(1, 4)]
    #df_cl_dist = pd.DataFrame(columns=cols_cl_dist)
    clamp_dist = None
    try:
        clamp_dist = calc_charge_clamp_dist(pdb_file, clamp_resid)
    except KeyError:
        if protein_name:
            print(f'{protein_name}: KeyError while calculating clamp dist')
        else:
            print('KeyError while calculating clamp dist')

    except ValueError as e:
        if protein_name:
            print(f'{protein_name}: {e}')
        else:
            print(e)

    cl_dist = [protein_name]
    if clamp_dist is not None:
        for elem in clamp_dist:
            cl_dist.append(clamp_dist[elem])
    df_cl_dist = pd.DataFrame([cl_dist], columns=cols_cl_dist)
    #df_cl_dist = pd.concat([df_cl_dist, pd.Series(cl_dist, index=cols_cl_dist[0:len(cl_dist)])], ignore_index=True)
    return df_cl_dist
