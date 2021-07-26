# plddt2csv
Extract the plddt value of Alphafold2 from .pkl to .csv & .pdb file.

This program was only checked by my environment. (centos7, python3)

Introduction movie was upload on Youtube.

https://youtu.be/0ozOPnC6IWk

Update information
--------------

plddt_align.py was released

    $ cd /path/to/output/of/Alphafold2 

    $ python3 /path/to/storage/plddt2csv/plddt_align.py

This program automatically generates [align_xxxx_?.pdb]s.

These .pdb files will be aligned each other.

These .pdb including pLDDT instead of b-factor.

![image](https://user-images.githubusercontent.com/87903303/127006396-5ec09e1c-cbb6-4f8a-b914-d48c1e5c73c7.png)

--------------
-d arguments was introduced.

    python3 plddt2csv.py -d 1-30_53-

               or

    python3 plddt_align.py -d 1-30_53-

avarage pLDDT will be caliculate without residue number from 1 to 30 and from 53 to C-terminus. 

["current folder name"_?_pLDDT??.pdb] will be also generated. default value is ""

![image](https://user-images.githubusercontent.com/87903303/126935487-5393b821-93ae-4460-a4dc-9a1034ea904e.png)

--------------
How to Use
--------------
    $ cd /path/to/storage/

    $ git clone https://github.com/CYP152N1/plddt2csv

    $ cd /path/to/output/of/Alphafold2

    $ python3 /path/to/storage/plddt2csv/plddt2csv.py
    
                         or

    $ python3 /path/to/storage//plddt2csv/plddt_align.py

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

