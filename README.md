# PNcsp

This repository represents the efforts of the [TCCDEM](https://github.com/tccdem/) which is an implementation of Phase Prediction via Crystal Structure Similarity in the Periodic Number Pepresentation. 

This work introduces a novel strategy for crystal structure prediction founded upon the principle of chemical similarity.  Our method uses Mendeleev's Periodic Number (PN) as a quantitative measure of substitutability to identify potential crystal structures for unexplored chemical systems. Representation of the workflow for predicting stable materials based on PN similarity is shown below. 

![Graphical_Abstract_3](https://github.com/user-attachments/assets/cf590168-ab66-4dc7-8954-de794dfbf780)

PNcsp is a generaliezed version of this workflow. It scans the [OQMD](https://www.oqmd.org/) for a desired chemical system for given order of neighbors and proposes initial crystal structures that are ready for further analysis.

### Installation & Usage
This program is based on Python 3 under Anaconda. 

1) Clone the repository.
2) Open terminal and locate the directory.
3) Install reuqirements or run:
```bash
   pip install -r requirements.txt
```
4) Help page:
```bash
   python PNcsp.py -h
```
5) Run the Python code:
```bash
  python PNcsp.py <formula> -n <neighbor_order>  -f <energy_filter>
```
\<formula\>: Pretty formula of chemical system for query (Ex. Cu2Mn1Al1)

\<neighbor_order\>: Order of neighbors to be considered in the similarity search.. (Default: 1 (first order neighbors))

\<energy_filter\>: Selected neighbors are limited to those below the energy filter value. (default: 0) unit: [eV/atom]. Use "none" to disable filter.

Created prototypes are shown in "output" folder in current directory.

### Example usage
```bash
python Similarity.py Na2Cl1 -n 3 -f 0.1
```
