from Bio import PDB
import pandas as pd

from geodes import constraints, utils


def _calc_acc_per_hel(dssp, ref):
    """Calculating of solvent-accessibility area per helix, Å2."""

    helix_borders = []

    res_num = utils.getResidues(dssp)
    for i in ref:

        start = utils.getNum(i[0], res_num)
        end = utils.getNum(i[1], res_num)

        helix_borders.append([start, end])

    helices = dict()
    acc = dict()

    # for i in range(len(dssp.keys())):
    #    if dssp[list(dssp.keys())[i]][1] == 'X':
    #        print(i)
    #        print(dssp.keys()[i])
    #        print(dssp[list(dssp.keys())[i]][1])

    for i in range(len(helix_borders)):
        helices[i] = {
            el: dssp[list(dssp.keys())[el]][3] * constraints.RES_MAX_ACC[dssp[list(dssp.keys())[el]][1]]
            for el in range(helix_borders[i][0], helix_borders[i][1]+1)
        }

    # calculation of solvent accessibility per helix according to data from dssp-file
    for key in helices.keys():
        for res in helices[key].keys():
            if key in acc:
                acc[key] += helices[key][res]
            else:
                acc[key] = helices[key][res]

    for key, acc_sum in acc.items():
        acc[key] = acc_sum/(len(helices[key]))

    # return dict with accessibility for every alpha-helix in structure
    return acc


def calc_acc_per_hel(pdb_file, ref):
    """
    Calculating of solvent-accessibility area per helix, Å2. Requires DSSP module.

    Parameters
    ----------
    pdb_file: str
        Filename of .pdb file used for calculation.
    ref: list of ints
        List of amino acid numbers pairs (start, end) for each helix.
    Returns
    -------
    dict with accessibility for every alpha-helix in structure.
"""

    _, _, model, _, _ = utils.get_model_and_structure(pdb_file)
    dssp = PDB.DSSP(model, pdb_file)
    if not isinstance(ref, list):
        if ref is None:
            raise ValueError("Ref list is None!")
        else:
            raise ValueError(f"Unexpected type for ref: {type(ref)}")
    return _calc_acc_per_hel(dssp, ref)


def acc_per_hel_to_pandas(pdb_file, ref, protein_name=None, **kwargs):
    """Putting solvent-accessibility area per helix in pandas dataframe.

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
    cols_acc = ['prot_name'] + ['ACC H' + str(elem) for elem in range(1, len(ref)+1)]
    #df_acc = pd.DataFrame(columns=cols_acc)

    acc_hels = None
    try:
        acc_hels = calc_acc_per_hel(pdb_file, ref)
    except KeyError:
        if protein_name:
            print(f'{protein_name}: KeyError while calculating acc')
        else:
            print('KeyError while calculating acc')

    except ValueError as e:
        if protein_name:
            print(f'{protein_name}: {e}')
        else:
            print(e)

    data_acc = [protein_name]
    if acc_hels is not None:
        for acc in acc_hels:
            data_acc.append(acc_hels[acc])

    df_acc = pd.DataFrame([data_acc], columns=cols_acc)
    #df_acc = pd.concat([df_acc, pd.Series(data_acc, index=cols_acc[0:len(data_acc)])], ignore_index=True)
    return df_acc
