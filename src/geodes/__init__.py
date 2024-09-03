'''
This repository contains GEODES python package.
This tool provides the set of functions for calculating various geometrical descriptors for proteins
given their structure as a PDB-file and outputting them in a form of a dataframe or CSV-table.

To use this package, clone this repository and enter the following command from the location of this project

`pip install . -r requirements.txt`

Several descriptors require DSSP module to be installed, this can be done using:

`sudo apt-get install dssp`

Additionally one may need KPAX software for providing user settings for protein helices borders: http://kpax.loria.fr/

To run the calculations one needs to have the data in PDB format and its secondary structure annotation in the
list format.
Charge clamp residues are optional, one can set up a personal configuration file to turn off the calculation of these
descriptors.

Detailed examples of usage and the following analysis are located at `usage_example.ipynb` and `analysis_example.ipynb`.

'''

from geodes.dssp_acc_hel import calc_acc_per_hel, acc_per_hel_to_pandas
from geodes.angle_hel_pairwise import calc_angles_between_hel, angles_between_hel_to_pandas
from geodes.angle_charge_clamp_Ca import calc_charge_clamp_angles, charge_clamp_angles_to_pandas
from geodes.dist_charge_clamp_Ca import calc_charge_clamp_dist, charge_clamp_dist_to_pandas
from geodes.angle_COMprot_hel_Ca_endpoints import calc_COM_Calpha_angles, COM_Calpha_angles_to_pandas
from geodes.dist_charge_clamp_Ca_COMprot import calc_COM_clamp, COM_clamp_to_pandas
from geodes.COM_helix import calc_COM_helix, COM_helix_to_pandas
from geodes.COM_protein import calc_COM_protein, COM_protein_to_pandas
from geodes.dssp_hel_endpoints import calc_dssp_hel, dssp_hel_to_pandas, dssp_extra_to_pandas
from geodes.dist_hel_Ca_endpoints import calc_len_of_hel, len_of_hel_to_pandas
from geodes.dist_COMhel_pairwise import calc_pairwise_sep_dist, pairwise_sep_dist_to_pandas
from geodes.dist_COMprot_COMhel import calc_prot_hel_dist, prot_hel_dist_to_pandas
from geodes.dssp_sse_content import calc_sse_content, sse_content_to_pandas
from geodes import constraints, utils
from geodes.main import DescCalculator


__all__ = [
    "calc_acc_per_hel", 'acc_per_hel_to_pandas',
    'calc_angles_between_hel', 'angles_between_hel_to_pandas',
    'calc_charge_clamp_angles', 'charge_clamp_angles_to_pandas',
    'calc_charge_clamp_dist', 'charge_clamp_dist_to_pandas',
    'calc_COM_Calpha_angles', 'COM_Calpha_angles_to_pandas',
    'calc_COM_clamp', 'COM_clamp_to_pandas',
    'calc_COM_helix', 'COM_helix_to_pandas',
    'calc_COM_protein', 'COM_protein_to_pandas',
    'calc_dssp_hel', 'dssp_hel_to_pandas', 'dssp_extra_to_pandas',
    'calc_len_of_hel', 'len_of_hel_to_pandas',
    'calc_pairwise_sep_dist', 'pairwise_sep_dist_to_pandas',
    'calc_prot_hel_dist', 'prot_hel_dist_to_pandas',
    'calc_sse_content', 'sse_content_to_pandas',
    'constraints', 'utils', 'DescCalculator',
]
