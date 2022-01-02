# plddt2csv
Extract the PAE & pLDDT value of Alphafold ver. 2.0.0 from .pkl to .csv & .pdb file.

Alinment program of .pdb files were updated for Alphafold ver. 2.1.1

https://github.com/CYP152N1/AF2align

--------------

Introduction movie was upload on Youtube.

Japanese var. (plddt_align.py)

https://youtu.be/qiw5DbzZcUE

No sound var. (plddt2csv.py)

https://youtu.be/0ozOPnC6IWk

for Alphafold 2.1.1 (align_mono.py, align_ptm.py, align_multi.py)

https://youtu.be/sDsV7pr6T_A?t=358

Update information
--------------

"align_multi.py", "align_mono.py" and "align_ptm.py" were released.

This version was created for alphafold 2.1.1. to align .pdb files  

This script can be parformed without any argument.

    $ cd /path/to/output/directory
    
    $ python align_mono.py

It is omitted that the inthersion of pLDDT at b-factor on ".pdb format", 
because pLDDT was already writtened in out put file of AF 2.1.    

additional argument can be observed by 
    
    $ python align_mono.py -h

multimer version can be select alinment chain by

    $ python align_multi.py -alc 2

The default value of -arc argument is 1 for the alignment of A chain. 

These alignment program will be separated to another programs, because program content and name don't match.

--------------

"pae2csv.py", "pae2png.py", "plddt_align_ptm.py", and "ptm_docker.py" were released

[pae2csv.py] generates [p2c_pae_?.csv] & [p2x_plddt.csv].

[p2c_pae_?.csv] contains Predicted Aligned Error (PAE). 

You can visualize using Excel like the figure shown below.

![image](https://user-images.githubusercontent.com/87903303/127757080-41d251db-6b8f-4499-bf2a-f0ab71b953c9.png)

Thanks Dr. Yoshitaka Moriwaki (@Ag_smith) for showing how to export predicted aligned error in alphafold2.


To extract "PAE" into [result_model_?.pkl], we need to change [docker/run_docker.py]

[ptm_docker.py] creates run_docker_ptm.py from your run_docker.py.

model_names of [run_docker_ptm.py] was changed to "model_?_ptm" from "model_?" 

    $ cd /path/to/storage/

    $ git clone https://github.com/CYP152N1/plddt2csv

    $ cd /path/to/Alphafold2
   
    $ python3 /path/to/Alphafold2/ptm_docker.py
    
    $ python3 docker/run_docker_ptm.py --fasta_paths=/path/to/.fasta --max_template_date=2020-05-14
    
If you performe "run_docker_ptm.py" instead of "run_docker.py", 

result_model_?_ptm.pkl posesses predicted_aligned_error value.
    
    $ cd /path/to/output/of/Alphafold2/
    
    $ python3 /path/to/storage/plddt2csv/pae2csv.py

               or

    $ python3 /path/to/storage/plddt2csv/plddt_align_ptm.py
    
If you perform "run_docker_ptm.py" instead of "run_docker.py"_ptm, 

you cannot use plddt_align.py and plddt2csv.py

![image](https://user-images.githubusercontent.com/87903303/127762914-224b0ddb-b785-4263-aa51-0fb606e1176f.png)


[pae2png.py] generates [p2c_pae_?.png] using matplotlib & seaborn.

You need to perform below comands.

    $ pip3 install matplotlib
    
    $ pip3 install seaborn

![image](https://user-images.githubusercontent.com/87903303/127762852-eb4047fd-8df9-40d7-b1a5-ad44b2e5aae1.png)

--------------

-m arguments was introduced in plddt_align.py & plddt_align_ptm.py.

If you don't perform amber_minimization,

or amber_minimization was clashed on first model

    python3 plddt_align.py -m unrelaxed_model

               or

    python3 plddt_align.py -m unrelaxed_model -n 1

--------------
plddt_align.py was released

    $ cd /path/to/output/of/Alphafold2 

    $ python3 /path/to/storage/plddt2csv/plddt_align.py

This program automatically generates [align_xxxx_?.pdb]s.

These .pdb files will be aligned each other.

These .pdb including pLDDT instead of b-factor.

![image](https://user-images.githubusercontent.com/87903303/127757060-8d7081f6-6b15-4c92-b133-6988f1c50d66.png)

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

    $ python3 /path/to/storage/plddt2csv/plddt_align.py

--------------
This program generates [plddt.csv] and [relaxed_"current folder name"_?.pdb] or [align_"current folder name"_?.pdb]

form [features.pkl], [result_model_?.pkl] and [relaxed_model_?.pdb]

plddt.csv show plddt value of each residue of all relaxed models.

![image](https://user-images.githubusercontent.com/87903303/126872927-1b69da3d-35a6-4978-ac5b-22ddf87b9f69.png)

output pdb file contain plddt value in B factor region.
High plddt value means high condidence of the output structure.
You can check plddt put the colored by b-factor on the structure using model viewer such as pymol.  

![ppp](https://user-images.githubusercontent.com/87903303/126873215-8ed0f006-525d-45ae-abb4-ac8d0bab71f3.png)

arguments can be chack by $ Python3 /path/to/storage/plddt2csv.py -h

