import numpy as np
import pandas as pd

import utils


def _calc_angles_between_hel(chain, ref):
    """Calculation of angles between all helices in structure."""

    helix_borders = ref
    vect_helices = []
    for i in range(0, len(helix_borders)):
        vect_helices.append((chain[helix_borders[i][1]]['CA'].get_coord()-chain[helix_borders[i][0]]['CA'].get_coord()))

    cos_between_hel = []
    for j in range(0, len(vect_helices) - 1):
        inter_result = []
        for i in range(j + 1, len(vect_helices)):
            OHel1 = np.linalg.norm(vect_helices[j])
            OHel2 = np.linalg.norm(vect_helices[i])
            OHels = np.dot(vect_helices[i], vect_helices[j])
            inter_result.append(OHels/(OHel1*OHel2))

        cos_between_hel.append((np.degrees(np.arccos(inter_result))).tolist())

    return cos_between_hel


def calc_angles_between_hel(pdb_file, ref):
    """Calculation of angles between all helices in structure.

    Parameters
    ----------
    pdb_file: str
        Filename of .pdb file used for calculation.
    ref: list of ints
        List of amino acid numbers pairs (start, end) for each helix.

    Returns
    -------
    list of pairwise angles lists between helices.

    """
    _, _, _, chain, _ = utils.get_model_and_structure(pdb_file)
    if not isinstance(ref, list):
        if ref is None:
            raise ValueError("Ref list is None!")
        else:
            raise ValueError(f"Unexpected type for ref: {type(ref)}")
    return _calc_angles_between_hel(chain, ref)


def angles_between_hel_to_pandas(pdb_file, ref, protein_name=None, **kwargs):
    """Putting angles between all helices in pandas dataframe.

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
    cols_cos = ['prot_name'] + [f'Angle H{i}-H{j}' for i in range(1, 13) for j in range(i+1, 14)]
    df_cos = pd.DataFrame(columns=cols_cos)
    cos = None
    try:
        cos = calc_angles_between_hel(pdb_file, ref)

    except KeyError:
        if protein_name:
            print(f'{protein_name}: KeyError while calculating cos')
        else:
            print('KeyError while calculating cos')

    except ValueError as e:
        if protein_name:
            print(f'{protein_name}: {e}')
        else:
            print(e)

    data_cos = [protein_name]
    if cos is not None:
        for elem in cos:
            for angle in elem:
                data_cos.append(angle)
    df_cos = df_cos.append(pd.Series(data_cos, index=cols_cos[0:len(data_cos)]), ignore_index=True)
    return df_cos
