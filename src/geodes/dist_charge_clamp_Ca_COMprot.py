import numpy as np
import pandas as pd

from geodes import utils
from geodes.COM_protein import _calc_COM_protein


def _calc_COM_clamp(chain, atom_struct, charge_clamps):
    """Calculate distances between protein's center of mass and every charge clamps."""

    atom_coord = []
    com = np.array(_calc_COM_protein(atom_struct))  # We need here to calculate protein's center of mass

    # Find charge clamps from user input in PDB structure and get it coordinates
    for res in chain:
        if res.id[1] == charge_clamps[0] or res.id[1] == charge_clamps[1] or res.id[1] == charge_clamps[2]:
            for atom in res:
                if atom.get_name() == 'CA':
                    atom_coord.append(atom.get_coord())

    # Calculate distance between center of mass and every charge

    ch_clamp_dist = []
    for elem in atom_coord:
        vect = np.array(elem) - np.array(com)
        ch_clamp_dist.append(np.sqrt(np.sum(vect ** 2)))

    return ch_clamp_dist


def calc_COM_clamp(pdb_file, charge_clamps):
    """
    Calculate distances between protein's center of mass and every charge clamps.

    Parameters
    ----------
    pdb_file: str
        Filename of .pdb file used for calculation.
    charge_clamps: list of ints
        Charge clamp residues list.

    Returns
    -------
    list of distances between protein's center of mass and every charge clamps.

    """
    _, _, _, chain, atom_struct = utils.get_model_and_structure(pdb_file)

    if not isinstance(charge_clamps, list):
        if charge_clamps is None:
            raise ValueError("Charge clamp residues list is None!")
        else:
            raise ValueError(f"Unexpected type for Charge clamp: {type(charge_clamps)}")

    return _calc_COM_clamp(chain, atom_struct, charge_clamps)


def COM_clamp_to_pandas(pdb_file, clamp_resid, protein_name=None, **kwargs):
    """
    Putting distances between protein's center of mass and every charge clamps in pandas dataframe.

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
    cols_comclampdist = ['prot_name'] + [f'Dist COM-clamp{i}' for i in range(1, 4)]
    #df_clamps = pd.DataFrame(columns=cols_comclampdist)
    clamps = None

    try:
        clamps = calc_COM_clamp(pdb_file, clamp_resid)
    except KeyError:
        if protein_name:
            print(f'{protein_name}: KeyError while calculating COM-clamp')
        else:
            print('KeyError while calculating COM-clamp')

    except ValueError as e:
        if protein_name:
            print(f'{protein_name}: {e}')
        else:
            print(e)

    data_clamps = [protein_name]
    if clamps is not None:
        for dist in clamps:
            data_clamps.append(dist)
    df_clamps = pd.DataFrame([data_clamps], columns=cols_comclampdist)
    #df_clamps = pd.concat([df_clamps, pd.Series(data_clamps, index=cols_comclampdist[0:len(data_clamps)])], ignore_index=True)
    return df_clamps
