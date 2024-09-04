# GEODES 
This repository contains **GEODES** Python package. This tool is based upon the collection of Python modules for calculating various geometrical descriptors for proteins given their structure as a PDB-file and outputting them in a form of a dataframe or CSV-table.


## Installation

Several descriptors require DSSP module to be installed, this can be done using:

`conda install salilab::dssp`

The alternative installation using:

`sudo apt-get install dssp`

should be chosen with caution as currently the clean installation of DSSP 4.0.4 (intended for Ubuntu 22.04) does not work properly. It might be working on other versions though.

Additionally one may need KPAX software for providing user settings for protein helices borders: http://kpax.loria.fr/

## Usage

To see an example of usage of this package, see [`usage_example.ipynb`](/examples/usage_example.ipynb)


## Authors
- Karina Pats, Visiting Research Assistant, Nazarbayev University; PhD candidate, ITMO University;
- Igor Glukhov, MSc student (2022-2023), ITMO University
- Elizaveta Vinogradova, MSc student (2021-2022), ITMO University


## Acknowledgements
- Dr. Ferdinand Molnár, Associate Professor, Nazarbayev University
- Dr. Marie-Dominique Devignes, Lorraine Research Laboratory in Computer Science and its Applications, University of Lorraine
- Stepan Petrosyan, Bioinformatics Institute student (2019/2020)
- Maria Mamayeva, Bioinformatics Institute student (2019/2020)
