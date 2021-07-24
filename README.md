# plddt2csv
Extract the plddt value of Alphafold2 from .pkl to .csv file
This program was only checked my condition. (centos7, python3)

--------------
How to Use
--------------
$ cd /path/to/storage/

$ git clone https://github.com/CYP152N1/plddt2csv

$ cd /path/to/output/of/Alphafold2

$ python3 /path/to/storage/plddt2csv.py


![image](https://user-images.githubusercontent.com/87903303/126872882-a89776f9-50aa-434d-b773-f96bf775f2fb.png)

--------------
This program generates [plddt.csv] and [relaxed_"current folder name"_?.pdb] 

form [features.pkl], [result_model_?.pkl] and [relaxed_model_?.pdb]

plddt.csv show plddt value of each residue of all relaxed models.

![image](https://user-images.githubusercontent.com/87903303/126872927-1b69da3d-35a6-4978-ac5b-22ddf87b9f69.png)

output pdb file contain plddt value in B factor region.
High plddt value means high condidence of the output structure.
You can check plddt put the colored by b-factor on the structure using model viewer such as pymol.  

![ppp](https://user-images.githubusercontent.com/87903303/126873215-8ed0f006-525d-45ae-abb4-ac8d0bab71f3.png)

arguments can be chack by $ Python3 /path/to/storage/plddt2csv.py -h

