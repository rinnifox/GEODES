from Bio import PDB
import pandas as pd

import utils


def _calc_sse_content(dssp):
    """Calculation of secondary structure content."""

    # preparing for extruction sse from dssp structure
    resamount = len(dssp.keys()) + 1
    dssp_structures = ['H', 'B', 'E', 'G', 'I', 'T', 'S']
    sses = list()

    # extracting sse from dssp
    for i in range(0, len(dssp.keys())):
        if dssp[list(dssp.keys())[i]][2] in dssp_structures:
            sses.append(dssp[list(dssp.keys())[i]][2])

    # making dict of all possible secondary structures and counting their percentage
    sse = {'Helix': sses.count('H') / resamount * 100,
           'Beta bridge': sses.count('B') / resamount * 100,
           'Strand': sses.count('E') / resamount * 100,
           'Helix-3': sses.count('G') / resamount * 100,
           'Helix-5': sses.count('I') / resamount * 100,
           'Turn': sses.count('T') / resamount * 100,
           'Bend': sses.count('S') / resamount * 100,
           'Other': (resamount - len(sses)) / resamount * 100}

    return sse


def calc_sse_content(pdb_file):
    """
    Calculation of secondary structure content.

    Parameters
    ----------
    pdb_file: str
        Filename of .pdb file used for calculation.

    Returns
    -------
    dict of all possible secondary structures and counting their percentage.

    """
    _, _, model, _, _ = utils.get_model_and_structure(pdb_file)
    dssp = PDB.DSSP(model, pdb_file)
    return _calc_sse_content(dssp)


def sse_content_to_pandas(pdb_file, protein_name=None, **kwargs):
    """
    Putting secondary structure content in pandas dataframe.

    Parameters
    ----------
    pdb_file: str
        Filename of .pdb file used for calculation.
    protein_name: str, default=None
        Protein name to be added to the resulting dataframe.

    Returns
    -------
    pandas.DataFrame with calculated descriptor.

    """
    cols_sse = ['prot_name', 'SSE Helix', 'SSE Beta bridge',
                'SSE Strand', 'SSE Helix-3', 'SSE Helix-5',
                'SSE Turn', 'SSE Bend', 'SSE Other']
    df_sse = pd.DataFrame(columns=cols_sse)
    sse = None
    try:
        sse = calc_sse_content(pdb_file)
    except KeyError:
        if protein_name:
            print(f'{protein_name}: KeyError while calculating sse')
        else:
            print('KeyError while calculating sse')

    except ValueError as e:
        if protein_name:
            print(f'{protein_name}: {e}')
        else:
            print(e)

    data_sse = [protein_name]
    if sse is not None:
        for struct in sse:
            data_sse.append(sse[struct])
    df_sse = df_sse.append(pd.Series(data_sse, index=cols_sse[0:len(data_sse)]), ignore_index=True)
    return df_sse
