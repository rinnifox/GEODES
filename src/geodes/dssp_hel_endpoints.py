from Bio import PDB
import numpy as np
import pandas as pd

from geodes import utils


def _calc_dssp_hel(dssp, ref):
    """Calculate differences with dssp module."""
#     TODO: Split function into smaller functions
    chainA = [key for key in dssp.keys() if key[0] == 'A']
    helix_map = np.zeros([1, len(chainA)])
    res_num = utils.getResidues(dssp)
    dssp_start = 0
    dssp_end = 0
    result = []

    # print(res_num)

    for i in range(len(ref)):

        # print(ref[i][0])

        start = utils.getNum(ref[i][0], res_num)
        end = utils.getNum(ref[i][1], res_num)

        # finding starting point
        start_longer_counter = 0
        start_shorter_counter = 0

# TODO: wrap in single func

        if dssp[list(dssp.keys())[start]][2] in ['H', 'G']:
            # check the first iteration
            while dssp[list(dssp.keys())[start-1]][2] in ['H', 'G'] and utils.getRes(start-1, res_num) != dssp_end:
                start_longer_counter += 1
                start -= 1
            missing = False
        else:
            missing_counter = 0
            missing = True
            while missing_counter < (end-start):
                start += 1
                start_shorter_counter += 1
                if dssp[list(dssp.keys())[start]][2] in ['H', 'G']:
                    missing = False
                    break
                else:
                    missing_counter += 1

        # finding endpoint
        if missing is False:
            end_longer_counter = 0
            end_shorter_counter = 0
            if dssp[list(dssp.keys())[end]][2] in ['H', 'G']:
                if i != (len(ref)-1):

                    while ((dssp[list(dssp.keys())[end+1]][2] in ['H', 'G'])
                            and (end+1 != utils.getNum(ref[i+1][0], res_num))):
                        end_longer_counter += 1
                        end += 1
                else:
                    while dssp[list(dssp.keys())[end+1]][2] in ['H', 'G']:
                        end_longer_counter += 1
                        end += 1

                        try:
                            dssp[list(dssp.keys())[end+1]][2] in ['H', 'G']
                        except IndexError:
                            break

            else:

                while dssp[list(dssp.keys())[end]][2] not in ['H', 'G']:
                    end -= 1
                    end_shorter_counter += 1

            if start_shorter_counter > 0:
                dssp_start = ref[i][0] + start_shorter_counter
            else:
                dssp_start = ref[i][0] - start_longer_counter

            if end_shorter_counter > 0:
                dssp_end = ref[i][1] - end_shorter_counter
            else:
                dssp_end = ref[i][1] + end_longer_counter

            result.append([dssp_start, dssp_end])

            for i in range(start, end+1):
                helix_map[0][i] = 1

        else:
            result.append([0, 0])

    extras = []
    map_elem = 0

# TODO: wrap
    while map_elem < helix_map.shape[1]:
        if helix_map[0][map_elem] == 0:
            if dssp[list(dssp.keys())[map_elem]][2] in ['H', 'G']:
                extra_counter = map_elem

                while dssp[list(dssp.keys())[extra_counter+1]][2] in ['H', 'G']:
                    extra_counter += 1

                extras.append([utils.getRes(map_elem, res_num), utils.getRes(extra_counter, res_num)])

                if map_elem == extra_counter:
                    map_elem += 1
                else:
                    map_elem = extra_counter + 1
            else:
                map_elem += 1
        else:
            map_elem += 1

        n_res = 0
        for e in extras:
            n_res += e[1]-e[0]+1

    return result, n_res


def calc_dssp_hel(pdb_file, ref):
    """
    Calculates differences with DSSP output.

    Parameters
    ----------
    pdb_file: str
        Filename of .pdb file used for calculation.
    ref: list of lists (int, int)
        List of amino acid numbers pairs (start, end) for each helix.

    Returns
    -------
    tuple(list, int): first item (list) contains boundaries of helices
    second item (int) contains number of aminoacids that are from helices based on dssp,
    but are not considered as a part of helices in our reference.

    """
    _, _, model, _, _ = utils.get_model_and_structure(pdb_file)
    dssp = PDB.DSSP(model, pdb_file)
    if not isinstance(ref, list):
        if ref is None:
            raise ValueError("Ref list is None!")
        else:
            raise ValueError(f"Unexpected type for ref: {type(ref)}")
    return _calc_dssp_hel(dssp, ref)


def dssp_hel_to_pandas(pdb_file, ref, protein_name=None, **kwargs):
    """
    Putting differences in structure with dssp module in pandas dataframe.

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
    cols_dssp = (['prot_name']
                 + ['DSSP start_H' + str(elem) for elem in range(1, len(ref)+1)]
                 + ['DSSP end_H' + str(elem) for elem in range(1, len(ref)+1)])
    #df_dssp = pd.DataFrame(columns=cols_dssp)
    dssp_hels = None

    try:
        dssp_hels = calc_dssp_hel(pdb_file, ref)
    except KeyError:
        if protein_name:
            print(f'{protein_name}: KeyError while calculating dssp')
        else:
            print('KeyError while calculating dssp')

    except ValueError as e:
        if protein_name:
            print(f'{protein_name}: {e}')
        else:
            print(e)

    data_dssp_hels = [protein_name]
    if dssp_hels is not None:
        for hel in dssp_hels[0]:
            data_dssp_hels.append(hel[0])
            data_dssp_hels.append(hel[1])

    df_dssp = pd.DataFrame([data_dssp_hels], columns=cols_dssp)
    #df_dssp = pd.concat([df_dssp, pd.Series(data_dssp_hels, index=cols_dssp[0:len(data_dssp_hels)])], ignore_index=True)
    return df_dssp


def dssp_extra_to_pandas(pdb_file, ref, protein_name=None, **kwargs):
    """
    Putting differences with DSSP in pandas dataframe (extra).

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
    cols_extra_res = ['prot_name', 'N_res extra helical']
    #df_extra = pd.DataFrame(columns=cols_extra_res)
    dssp_hels = None

    try:
        dssp_hels = calc_dssp_hel(pdb_file, ref)
    except KeyError:
        if protein_name:
            print(f'{protein_name}: KeyError while calculating dssp')
        else:
            print('KeyError while calculating dssp')

    except ValueError as e:
        if protein_name:
            print(f'{protein_name}: {e}')
        else:
            print(e)

    data_extra_hels = [protein_name]
    if dssp_hels is not None:
        data_extra_hels.append(dssp_hels[1])
    df_extra = pd.DataFrame([data_extra_hels], columns=cols_extra_res)
    #df_extra = pd.concat([df_extra, pd.Series(data_extra_hels, index=cols_extra_res[0:len(data_extra_hels)])],
    #                           ignore_index=True)
    return df_extra
