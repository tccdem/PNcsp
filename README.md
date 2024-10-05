# PNcsp

This repository represents the efforts of the [TCCDEM](https://github.com/tccdem/) which is an implementation of Phase Prediction via Crystal Structure Similarity in the Periodic Number Pepresentation. 

This work introduces a novel strategy for crystal structure prediction founded upon the principle of chemical similarity.  Our method leverages Mendeleev's Periodic Number (PN) as a quantitative measure of substitutability to identify potential crystal structures for unexplored chemical systems. We apply this concept to predict stable equiatomic binary phases, demonstrating its efficacy in uncovering both known and novel structures.

![Graphical_Abstract_3](https://github.com/user-attachments/assets/cf590168-ab66-4dc7-8954-de794dfbf780)

Representation of the workflow for predicting stable materials based on PN similarity. Selected prototypes are subjected to further analyses for the investigation of electronic properties.



### Installation & Usage
This program use Python 3.12.3 under Anaconda3. 

1) Clone the repository.
2) Open the command line in the directory.
3) Install reuqirements or run:
```bash
   pip install -r requirements.txt
```
5) Run the Python code:
```bash
  python PNcsp.py <formula> -n <#_of_neighbor>  -f <energy_filter>
```
\<formula\>: Pretty formula of chemical system for query (Ex. Cu2Mn1Al1)

\<#_of_neighbor\>: How many neighbor will be taken into account in similarity search. (default: 1)

\<energy_filter\>: Selected neighbors can be limited with Energy Threshold parameter. (default: "<0") unit: [eV/atom]

Created prototypes are shown in "output" folder in current directory.

### Example usage
```bash
python Similarity.py Na2Cl1 -n 3 -f 0.1
```
