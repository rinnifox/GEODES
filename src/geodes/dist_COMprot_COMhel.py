import numpy as np
import pandas as pd

from COM_helix import _calc_COM_helix
from COM_protein import _calc_COM_protein
import utils


def _calc_prot_hel_dist(chain, atom_struct, ref):
    """Calculate distance between protein's center of mass and between every helix's center of mass."""

    # Calculate protein's center of mass and centers of helices masses
    hel_COMs = _calc_COM_helix(chain, ref)
    prot_COM = _calc_COM_protein(atom_struct)
    # Calculate distance between protein's COM and helices COMs as vector distance
    prot_hel_dists = []
    for i in range(0, len(hel_COMs)):
        vect = np.array(hel_COMs[i]) - np.array(prot_COM)
        prot_hel_dists.append(np.sqrt(np.sum(vect ** 2)))
    return prot_hel_dists


def calc_prot_hel_dist(pdb_file, ref):
    """
    Calculate distance between protein's center of mass and between every helix's center of mass.

    Parameters
    ----------
    pdb_file: str
        Filename of .pdb file used for calculation.
    ref: list of ints
        List of amino acid numbers pairs (start, end) for each helix.

    Returns
    -------
    list of  distances lists between helices and protein center of mass.

    """
    _, _, _, chain, atom_struct = utils.get_model_and_structure(pdb_file)

    if not isinstance(ref, list):
        if ref is None:
            raise ValueError("Ref list is None!")
        else:
            raise ValueError(f"Unexpected type for ref: {type(ref)}")

    return _calc_prot_hel_dist(chain, atom_struct, ref)


def prot_hel_dist_to_pandas(pdb_file, ref, protein_name=None, **kwargs):
    """
    Putting distance between protein's center of mass and between every helix's center of mass in pandas dataframe.

    Parameters
    ----------
    pdb_file: str
        Filename of .pdb file used for calculation.
    ref: list of ints
        List of amino acid numbers pairs (start, end) for each helix.
    protein_name: str, default=None
        Protein name to be added to the resulting dataframe.

    Returns
    -------
    pandas.DataFrame with calculated descriptor.

    """
    cols_protheldist = ['prot_name'] + ['Dist prot-H' + str(elem) for elem in range(1, 14)]
    df_prothel = pd.DataFrame(columns=cols_protheldist)
    prothel = None

    try:
        prothel = calc_prot_hel_dist(pdb_file, ref)
    except KeyError:
        if protein_name:
            print(f'{protein_name}: Error while calculating prot-helix distance')
        else:
            print('Error while calculating prot-helix distance')

    except ValueError as e:
        if protein_name:
            print(f'{protein_name}: {e}')
        else:
            print(e)

    data_prothel = [protein_name]
    if prothel is not None:
        for elem in prothel:
            data_prothel.append(elem)
    df_prothel = df_prothel.append(pd.Series(data_prothel, index=cols_protheldist[0:len(data_prothel)]),
                                   ignore_index=True)
    return df_prothel
