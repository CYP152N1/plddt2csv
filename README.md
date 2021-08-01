# plddt2csv
Extract the plddt value of Alphafold2 from .pkl to .csv & .pdb file.

This program was only checked by my environment. (centos7, python3)

Introduction movie was upload on Youtube.

Japanese var. (plddt_align.py)

https://youtu.be/qiw5DbzZcUE

No sound var. (plddt2csv.py)

https://youtu.be/0ozOPnC6IWk

Update information
--------------

"pae2csv.py", "plddt_align_ptm.py", and "ptm_docker.py" were released

[pae2csv.py] generates [p2c_pae_?.csv] & [p2x_plddt.csv].

[p2c_pae_?.csv] contains Predicted Aligned Error (PAE). 

You can visualize using Excel like the figure shown below.

![image](https://user-images.githubusercontent.com/87903303/127755745-abd6a406-1df9-4a08-8f80-75e1e5a979bb.png)

Thanks Dr. Yoshitaka Moriwaki (@Ag_smith) for initially showing how to export predicted aligned error in alphafold2.


To extract "PAE" into [result_model_?.pkl], we need to change [docker/run_docker.py]

[ptm_docker.py] creates run_docker_ptm.py from your run_docker.py.

model_names of [run_docker_ptm.py] was changed to "model_?_ptm" from "model_?" 

    $ cd /path/to/storage/

    $ git clone https://github.com/CYP152N1/plddt2csv

    $ cd /path/to/Alphafold2
   
    $ python3 /path/to/Alphafold2/ptm_docker.py
    
    $ python3 docker/run_docker_ptm.py --fasta_paths=/path/to/.fasta --max_template_date=2020-05-14
    
If you purforme "run_docker_ptm.py" instead of "run_docker.py", 

result_model_?_ptm.pkl posesses predicted_aligned_error value.
    
    $ cd /path/to/output/of/Alphafold2/
    
    $ python3 /path/to/storage/plddt2csv/pae2csv.py

               or

    $ python3 /path/to/storage/plddt2csv/plddt_align_ptm.py
    
If you perform "run_docker_ptm.py" instead of "run_docker.py"_ptm, 

you cannot use plddt_align.py and plddt2csv.py


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

![image](https://user-images.githubusercontent.com/87903303/127021043-ba8f4d1e-c168-49d1-8ece-482e8a917399.png)

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

