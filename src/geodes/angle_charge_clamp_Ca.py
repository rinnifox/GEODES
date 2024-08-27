from Bio import PDB
import numpy as np
import pandas as pd

import utils


def _calc_charge_clamp_angles(chain, charge_clamps):
    """Calculation of angles between charge clamp residues."""

    # extracting vectors with coordinates of every charge clamp residue
    clamp_vectors = {}

    for res in chain:
        if res.id[1] == charge_clamps[0] or res.id[1] == charge_clamps[1] or res.id[1] == charge_clamps[2]:
            for atom in res:
                if atom.get_name() == 'CA':
                    clamp_vectors[res.id[1]] = atom.get_vector()

    # calculation angles in triangle formed by charge clamp residues
    angles = {}
    for elem in range(len(clamp_vectors)):
        angles[f'{charge_clamps[elem]}-{charge_clamps[elem - 1]}-{charge_clamps[elem - 2]}'] = (
            np.degrees(
                PDB.vectors.calc_angle(
                    clamp_vectors[charge_clamps[elem]],
                    clamp_vectors[charge_clamps[elem - 1]],
                    clamp_vectors[charge_clamps[elem - 2]]
                )
            )
        )

    return angles


def calc_charge_clamp_angles(pdb_file, charge_clamps):
    """
    Calculation of angles between charge clamp residues.

    Parameters
    ----------
    pdb_file: str
        Filename of .pdb file used for calculation.
    charge_clamps: list of ints
        Charge clamp residues list.

    Returns
    -------
    dict of angles.

    """

    _, _, _, chain, _ = utils.get_model_and_structure(pdb_file)
    if not isinstance(charge_clamps, list):
        if charge_clamps is None:
            raise ValueError("Charge clamp residues list is None!")
        else:
            raise ValueError(f"Unexpected type for Charge clamp: {type(charge_clamps)}")
    return _calc_charge_clamp_angles(chain, charge_clamps)


def charge_clamp_angles_to_pandas(pdb_file, clamp_resid, protein_name=None, **kwargs):
    """
    Putting angles between charge clamp residues in pandas dataframe.

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
    cols_cl_angle = ['prot_name'] + [f'Angle clamp{elem}-{el}' for elem in range(1, 3) for el in range(elem+1, 4)]
    df_cl_angles = pd.DataFrame(columns=cols_cl_angle)
    clamp_angle = None
    try:
        clamp_angle = calc_charge_clamp_angles(pdb_file, clamp_resid)
    except KeyError:
        if protein_name:
            print(f'{protein_name}: KeyError while calculating clamp angle')
        else:
            print('KeyError while calculating clamp angle')

    except ValueError as e:
        if protein_name:
            print(f'{protein_name}: {e}')
        else:
            print(e)
    cl_angle = [protein_name]
    if clamp_angle is not None:
        for elem in clamp_angle:
            cl_angle.append(clamp_angle[elem])
    df_cl_angles = df_cl_angles.append(pd.Series(cl_angle, index=cols_cl_angle[0:len(cl_angle)]), ignore_index=True)
    return df_cl_angles
