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
print("Ver.2_(2021Jul.)")
print("by Hiroki Onoda, YCU")
print("------------")
print("The program extract pae and plddt values from [result_model_?_ptm.pkl] to [.csv] ")
print("You need to change model name to model_?_ptm from model_? at Docker/run_docker.py")
print("for Alphafold2")
print("------------")
# Header

# Argument & Help
parser = argparse.ArgumentParser(description="The program extract plddt from [result_model.pkl] to [.csv]")
parser.add_argument('-i', default='result_model', help='input XXXXX if [XXXXX_Y.pkl]; ex) result_model (def.)')
parser.add_argument('-n', default='5', help='input Y if you want open [XXXXX_1.pkl,XXXXX_2.pkl ... ,XXXXX_5.pkl]; ex) 5 (def.)')
parser.add_argument('-o', default='p2c', help='input XXXXX if you want to name the output [XXXXX.csv]; ex) plddt (def.)')
parser.add_argument('-d', default='', help='ignored residue number for average pLDDTs caliculation to determined the ranking of models [1-5_56-75_105-135]; ex) "" (def.)')
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

# Number of models
namen=int(str(args.n))

# Output file name
oname=str(args.o)
print("Output file name:"+str(oname)+".csv")
print("It can be changed by [-o "+oname+"]")
print("")

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

    
    plddts=np.empty(namen)
    # Repeat for the caliculation of all models
    for i in range(namen):
        www=str(int(float(i)+1))
        # Read pLDDT from result_model_?.pkl
        with open(cpass+"/"+iname+"_"+www+"_ptm.pkl", 'rb') as c1:
            c1n=pickle.load(c1)
            # Extract PAE
            pae = c1n['predicted_aligned_error'].astype(np.float32)
            np.savetxt(cpass+'/'+oname+'_pae_'+www+'.csv',pae,delimiter=',',fmt="%s")
            print("Save:"+cpass+'/'+oname+'_pae_'+www+'.csv')
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
    np.savetxt(cpass+'/'+oname+'_plddt.csv',twoda,delimiter=',',fmt="%s")
    print("Save:"+cpass+'/'+oname+'_plddt.csv')
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
