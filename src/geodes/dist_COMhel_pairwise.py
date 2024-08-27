import itertools
import numpy as np
import pandas as pd

from COM_helix import _calc_COM_helix
import utils


def _calc_pairwise_sep_dist(chain, ref):
    """Calculate separation distance between every helix."""

    # Calculate centers of mass of every helix
    hel_COMs = list(itertools.chain(*_calc_COM_helix(chain, ref)))
    # Calculate pairwise separations between every COM's of every helix as vector distance
    pairwise_seps = []
    for j in range(0, len(hel_COMs) - 1):
        result = []
        for i in range(j + 1, len(hel_COMs)):
            vect = np.array(hel_COMs[i]) - np.array(hel_COMs[j])
            result.append(np.sqrt(np.sum(vect ** 2)))
        pairwise_seps.append(result)
    return pairwise_seps


def calc_pairwise_sep_dist(pdb_file, ref):
    """
    Calculate separation distance between every helix.

    Parameters
    ----------
    pdb_file: str
        Filename of .pdb file used for calculation.
    ref: list of ints
        List of amino acid numbers pairs (start, end) for each helix.

    Returns
    -------
    list of pairwise separation distances lists between helices.

    """
    _, _, _, chain, _ = utils.get_model_and_structure(pdb_file)

    if not isinstance(ref, list):
        if ref is None:
            raise ValueError("Ref list is None!")
        else:
            raise ValueError(f"Unexpected type for ref: {type(ref)}")

    return _calc_pairwise_sep_dist(chain, ref)


def pairwise_sep_dist_to_pandas(pdb_file, ref, protein_name=None, **kwargs):
    """
    Putting separation distance between every helix in pandas dataframe.

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
    cols_pairwise = (
        ['prot_name']
        + ['PairwiseSep H' + str(i) + '-H' + str(j)
            for i in range(1, 14)
            for j in range(i+1, 14)]
    )
    df_pairseps = pd.DataFrame(columns=cols_pairwise)
    pairseps = None

    try:
        pairseps = calc_pairwise_sep_dist(pdb_file, ref)
    except KeyError:
        if protein_name:
            print(f'{protein_name}: KeyError while calculating pairwise sep dist')
        else:
            print('KeyError while calculating pairwise sep dist')

    except ValueError as e:
        if protein_name:
            print(f'{protein_name}: {e}')
        else:
            print(e)

    data_pairseps = [protein_name]
    if pairseps is not None:
        for elem in pairseps:
            for dist in elem:
                data_pairseps.append(dist)
    df_pairseps = df_pairseps.append(pd.Series(data_pairseps, index=cols_pairwise[0:len(data_pairseps)]),
                                     ignore_index=True)
    return df_pairseps
