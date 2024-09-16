# GEODES 
This repository contains **GEODES** Python package. This tool is based upon the collection of Python modules for calculating various geometrical descriptors for proteins given their structure as a PDB-file and outputting them in a form of a dataframe or CSV-table.


## Installation and pre-requisites
The libraries that are necessary for using GEODES can be installed by downloading this repository and running:

`pip install -r requirements.txt `

Several descriptors also require DSSP module, this can be done using:

`conda install salilab::dssp`

The alternative installation using:

`sudo apt-get install dssp`

should be chosen with caution as currently the clean installation of DSSP 4.0.4 (intended for Ubuntu 22.04) does not work properly. It might be working on other versions though.

Additionally one may need KPAX software for providing user settings for protein helices borders: http://kpax.loria.fr/

GEODES package can be installed by running this command from the root of the repository:

`pip install .`

## Usage

To see an example of usage of this package, see [`usage_example.ipynb`](/examples/usage_example.ipynb)

To see an example of result analysis, see[`analysis_example.ipynb`](/examples/geodes-analysis/analysis_example.ipynb)


## Estimated calculation time
To run calculations in non-parallel mode, it takes ~3 min for 100 PDB files using Ubuntu 22.04 on virtual machine with 2048 MB RAM and 4 CPUs (cores).

## Authors
- Karina Pats, Visiting Research Assistant, Nazarbayev University; PhD candidate, ITMO University;
- Igor Glukhov, MSc student (2022-2023), ITMO University


## Acknowledgements
- Dr. Ferdinand Molnár, Associate Professor, Nazarbayev University
- Dr. Marie-Dominique Devignes, Lorraine Research Laboratory in Computer Science and its Applications, University of Lorraine
- Stepan Petrosyan, Bioinformatics Institute student (2019/2020)
- Maria Mamaeva, Bioinformatics Institute student (2019/2020)
