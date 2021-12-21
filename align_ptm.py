# Import
import os
import sys
import pickle
import csv
import sys
import time
import numpy as np
import argparse
from scipy.stats import rankdata
# Import

# Header
now = time.ctime()
cnvtime = time.strptime(now)
print("------------")
print(time.strftime("%Y/%m/%d %H:%M", cnvtime))
print("Start:"+str(sys.argv))
print("Ver.2.1_(2021Dec.)")
print("by Hiroki Onoda, YCU")
print("------------")
print("The program align output pdb file  and extract pae plddt values from [result_model_ptm.pkl] to [.csv] ")
print("Thanks Dr. Yoshitaka Moriwaki (@Ag_smith) for initially showing how to export predicted aligned error in alphafold2.")
print("You need to change model_names to model_?_ptm from model_? at /Docker/run_docker.py.")
print("for Alphafold2.1.1")
print("------------")
# Header

# Argument & Help
parser = argparse.ArgumentParser(            description="The program extract pae from [result_model.pkl] to [.csv]")
parser.add_argument('-i', default='result_model',   help='input XXXXX if [XXXXX_Y_ptm.pkl]; ex) result_model (def.)')
parser.add_argument('-m', default='relaxed_model',  help='input ZZZZZ if [ZZZZZ_Y_ptm.pdb]; ex) result_model (def.)')
parser.add_argument('-n', default='5',              help='input Y if you want open [XXXXX_1_ptm.pkl,XXXXX_2_ptm.pkl ... ,XXXXX_5_ptm.pkl]; ex) 5 (def.)')
parser.add_argument('-d', default='',               help='ignored residue number for average pLDDTs caliculation to determined the ranking of models [1-5_56-75_105-135]; ex) "" (def.)')
parser.add_argument('-al1', default=0,              help='Center the mid point of al1 and al3, the point (al1) shift on x-axis')
parser.add_argument('-al2', default=0,              help='The point (al2) shift on xy-plane')
parser.add_argument('-al3', default=0,              help='Center the mid point of al1 and al3, the point (al3) shift on x-asis and opozit site of al1')
args = parser.parse_args()
# Argument & Help

#l----------------------------------------------l
#l           Valiable information               l
#l----------------------------------------------l

# Input file name
iname=str(args.i)
print("Input file name:"+str(iname)+"_?_ptm.pkl")
print("It can be changed by [-i "+iname+"]")
print("")

mname=str(args.m)
print("Input file name:"+str(mname)+"_?_ptm.pdb")
print("It can be changed by [-m "+iname+"]")
print("")

# Number of models
namen=int(str(args.n))

# Ignore residue
dres=str(args.d)
ures=dres.split("_")
print("Ignore residues:"+str(ures))
sres=[]
for i in range(len(ures)):
    hres=[]
    hres=ures[i].split("-")
    sres.extend(hres)
fres=list(filter(None,sres))
mres=np.empty_like(fres)
mres[::1]=0
mres[::2]=1
fres_f = np.array(fres).astype(np.float32)
mres_f = mres.astype(np.float32)
ares=fres_f-mres_f
print("ON/OFF pisition:"+str(ares))

#l-----------------------------------------------l
#l           Valiable information                l
#l-----------------------------------------------l


def read():    
    cpass=os.getcwd()
    lpass=cpass.split("/")
    
    # Get sequene from features.pkl 
    with open(cpass+"/features.pkl", 'rb') as cn:
        cnn=pickle.load(cn)
        sequ=cnn['residue_index']
        rein=np.array(list(str(cnn['sequence'][0])[2:-1]))
        print(str(rein))
        twoda=np.column_stack((sequ+1,rein))
        cn.close()
        pass
    # Get sequence from features.pkl
    
    # Make array to determine the caliculation range of pLDDT average
    eres=np.empty_like(rein)
    for i in range(len(eres)):
        cococo=0
        for j in range(len(ares)):
            if i > int(ares[j])-1:
                cococo+=1
                pass
            pass
        if cococo % 2 == 0:
            eres[i]=1
        else:
            eres[i]=0
        pass
    # Make array to determine the caliculation range of pLDDT average
    print(eres)
    eres_f = eres.astype(np.float32)

    np.set_printoptions(formatter={'float': '{:.1f}'.format})

    # Model alignment
    import math
    o4=int(len(eres)/4)
    w4=int(o4*2)
    t4=int(o4*3)
    
    if int(args.al1)>0:
        if int(args.al1) < len(eres):
            o4=int(args.al1)
    if int(args.al1)>0:
        if int(args.al2) < len(eres):
            w4=int(args.al2)
    if int(args.al1)>0:
        if int(args.al3) < len(eres):
            t4=int(args.al3)
    
    xyz1=np.zeros((4,3),dtype=float)
    for i in range(namen):
        www=str(int(float(i)+1))
        with open(cpass+"/"+mname+"_"+www+"_ptm.pdb", mode='r') as g1:
            for line in g1:
                if line[0:4]=="ATOM":
                    if o4==int(float(line[22:26].strip())):
                        if "CA"==line[12:16].strip():
                            print(line)
                            xyz1[0,0] = float(line[30:38].strip())
                            xyz1[0,1] = float(line[38:46].strip())
                            xyz1[0,2] = float(line[46:54].strip())
                            pass
                        pass
                    elif w4==int(float(line[22:26].strip())):
                        if "CA"==line[12:16].strip():
                            print(line)
                            xyz1[1,0] = float(line[30:38].strip())
                            xyz1[1,1] = float(line[38:46].strip())
                            xyz1[1,2] = float(line[46:54].strip())
                            pass
                        pass
                    elif t4==int(float(line[22:26].strip())):
                        if "CA"==line[12:16].strip():
                            print(line)
                            xyz1[2,0] = float(line[30:38].strip())
                            xyz1[2,1] = float(line[38:46].strip())
                            xyz1[2,2] = float(line[46:54].strip())
                            pass
                        pass
                    else:
                        pass
                    pass
                else:
                    pass
                pass
            pass
        print("------------")
        xyz1[3,0] = round((xyz1[0,0]+xyz1[2,0])/2,3)
        xyz1[3,1] = round((xyz1[0,1]+xyz1[2,1])/2,3)
        xyz1[3,2] = round((xyz1[0,2]+xyz1[2,2])/2,3)
        xyzm=np.zeros((1,3),dtype=float)
        xyzm[0,0] = xyz1[3,0]
        xyzm[0,1] = xyz1[3,1]
        xyzm[0,2] = xyz1[3,2]
        print("Initial xyz")
        print(xyz1)
        xyz2=xyz1-xyzm
        print("Center xyz")
        #print(xyz2)
        #print("------------")
        # rotation z fit oy and ty to zero
        if xyz2[0,1]>0:
            pmy=-1
            pass
        else:
            pmy=1
            pass
        distaxy=math.sqrt(math.pow(xyz2[0,0], 2)+math.pow(xyz2[0,1], 2))
        radax=pmy*math.acos(xyz2[0,0]/distaxy)
        rot1=np.zeros((3,3),dtype=float)
        rot1[0,0]=math.cos(radax)
        rot1[1,0]=-math.sin(radax)
        rot1[0,1]=math.sin(radax)
        rot1[1,1]=math.cos(radax)
        rot1[2,2]=1
        xyz3=np.dot(xyz2,rot1)
        print("rotation z fit 1y and 3y to zero")
        #print(xyz3)
        #print("------------")
        # rotation y fit ox to zero
        if xyz3[0,2]>0:
            pmz=-1
            pass
        else:
            pmz=1
            pass
        dista1=math.sqrt(math.pow(xyz3[0,0], 2)+math.pow(xyz3[0,1], 2)+math.pow(xyz3[0,2], 2))
        radaxy=pmz*math.acos(xyz3[0,0]/dista1)
        rot2=np.zeros((3,3),dtype=float)
        rot2[0,0]=math.cos(radaxy)
        rot2[2,0]=-math.sin(radaxy)
        rot2[0,2]=math.sin(radaxy)
        rot2[2,2]=math.cos(radaxy)
        rot2[1,1]=1
        xyz4=np.dot(xyz3,rot2)
        print("rotation y fit 1z and 3z to zero")
        #print(xyz4)
        #print("------------")
        # rotation x fit wz to zero
        if xyz4[1,2]>0:
            pmyz=-1
            pass
        else:
            pmyz=1
            pass
        distbbyz=math.sqrt(math.pow(xyz4[1,1], 2)+math.pow(xyz4[1,2], 2))
        radbbyz=pmyz*math.acos(xyz4[1,1]/distbbyz)
        rot3=np.zeros((3,3),dtype=float)
        rot3[1,1]=math.cos(radbbyz)
        rot3[2,1]=-math.sin(radbbyz)
        rot3[1,2]=math.sin(radbbyz)
        rot3[2,2]=math.cos(radbbyz)
        rot3[0,0]=1
        xyz5=np.dot(xyz4,rot3)
        print("rotation x fit 2z to zero")
        #print(xyz5)
        #print("------------")
        # Check all shift & roation 
        xyz6=np.dot(np.dot(np.dot(xyz1-xyzm,rot1),rot2),rot3)
        print("Check all shift & roation")
        print(xyz6)
        print("------------")
        xyz=np.zeros((1,3),dtype=float)
        liii=""
        with open(cpass+"/"+mname+"_"+www+"_ptm.pdb", mode='r') as g1:
            for line in g1:
                if line[0:4]=="ATOM":
                    xyz[::1]=0
                    xyz[0,0] = float(line[30:38].strip())
                    xyz[0,1] = float(line[38:46].strip())
                    xyz[0,2] = float(line[46:54].strip())
                    xyzn=np.dot(np.dot(np.dot(xyz-xyzm,rot1),rot2),rot3)
                    qxxl=" "*(8-len(str(round(xyzn[0,0],4))))+str(round(xyzn[0,0],4))
                    qyyl=" "*(8-len(str(round(xyzn[0,1],4))))+str(round(xyzn[0,1],4))
                    qzzl=" "*(8-len(str(round(xyzn[0,2],4))))+str(round(xyzn[0,2],4))
                    liii+=line[0:30]+qxxl+qyyl+qzzl+line[54:]
                    pass
                else:
                    liii+=line
                pass
            pass
        with open(cpass+"/align_"+str(lpass[-1])+"_m"+www+"_ptm.pdb", mode='w') as i:
            i.write(liii)
            i.close()
            pass
        pass


    
    plddts=np.empty(namen)
    # Repeat for the caliculation of all models
    for i in range(namen):
        www=str(int(float(i)+1))
        # Read pLDDT from result_model_?.pkl
        with open(cpass+"/"+iname+"_"+www+"_ptm.pkl", 'rb') as c1:
            c1n=pickle.load(c1)
            # Extract PAE
            pae = c1n['predicted_aligned_error'].astype(np.float32)
            np.savetxt(cpass+'/'+str(lpass[-1])+'_pae_m'+www+'_ptm.csv',pae,delimiter=',',fmt="%s")
            print("Save:"+cpass+'/'+str(lpass[-1])+'_pae_m'+www+'_ptm.csv')
            # Extract PAE
            c1n_f = c1n['plddt'].astype(np.float32)
            # Caliculation of average pLDDT
            na_mul=c1n_f*eres_f
            plddts[i]=np.sum(na_mul)/np.sum(eres_f)
            # Caliculation of average pLDDT
            pld1=np.round(c1n['plddt'], decimals=2)
            twoda=np.column_stack((twoda,pld1))
            pd1=['{:.2f}'.format(n) for n in pld1.tolist()]
            c1.close()
            pass
        # Read pLDDT from result_model_?.pkl
        
        # Write pLDDT in the PDB_Bfactor
        print("Save:"+cpass+"/align_"+str(lpass[-1])+"_m"+www+"_ptm.pdb")
        pass
    print("------------")
    
    

    
    # Caliculation of pLDDT rank without specificic range (-d argument)
    print("number of calculated residues: "+str(int(np.sum(eres_f)))+" (Ignore: "+str(ures)+")")
    np.set_printoptions(formatter={'float': '{:.0f}'.format})
    print("Rank  :"+str(rankdata(-plddts)))
    np.set_printoptions(formatter={'float': '{:.1f}'.format})
    twoda=np.column_stack((twoda,eres))
    print("plddts:"+str(plddts))
    print("------------")
    # Caliculation of pLDDT rank without specificic range (-d argument)
    
    # Write pLDDT in csv file  
    np.savetxt(cpass+'/'+str(lpass[-1])+'_plddt.csv',twoda,delimiter=',',fmt="%s")
    print("Save:"+cpass+'/'+str(lpass[-1])+'_plddt.csv')
    pass



#l-----------------------------------------------l
#l                   Execution                   l
#l-----------------------------------------------l

print("------------")
read()
print("------------")
now = time.ctime()
cnvtime = time.strptime(now)
print(time.strftime("%Y/%m/%d %H:%M", cnvtime))
print("Finish:"+str(sys.argv))
print("------------")
