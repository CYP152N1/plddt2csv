# plddt2csv
Extract the plddt value of Alphafold2 from .pkl to .csv file
This program was only checked my condition. (centos7, python3)

--------------
How to Use
--------------
$ cd /path/to/storage/

$ git clone https://github.com/CYP152N1/plddt2csv

$ cd /path/to/output/of/Alphafold2

$Python3 /path/to/storage/plddt2csv.py

--------------
This program generates [plddt.csv] and [relaxed_"current folder name"_?.pdb] 
form [features.pkl], [result_model_?.pkl] and [relaxed_model_?.pdb]

plddt.csv show plddt value of each residue of all relaxed models.
output pdb file contain plddt value in B factor region.
High plddt value means high condidence of the output structure.
You can check plddt put the colored by b-factor on the structure using model viewer such as pymol.  

arguments can be chack by $ Python3 /path/to/storage/plddt2csv.py -h
