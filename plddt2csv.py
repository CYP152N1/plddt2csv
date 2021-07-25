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
print("The program extract plddt values from [result_model.pkl] to [.csv] ")
print("for Alphafold")
print("------------")
# Header

# Argument & Help
parser = argparse.ArgumentParser(description="The program extract plddt from [result_model.pkl] to [.csv]")
parser.add_argument('-i', default='result_model', help='input XXXXX if [XXXXX_Y.pkl]; ex) result_model (def.)')
parser.add_argument('-n', default='5', help='input Y if you want open [XXXXX_1.pkl,XXXXX_2.pkl ... ,XXXXX_5.pkl]; ex) 5 (def.)')
parser.add_argument('-o', default='plddt', help='input XXXXX if you want to name the output [XXXXX.csv]; ex) plddt (def.)')
parser.add_argument('-d', default='', help='ignoreed residue number for average pLDDTs caliculation to determined the ranking of models [1-5_56-75_105-135]; ex) "" (def.)')
args = parser.parse_args()
# Argument & Help

#l----------------------------------------------l
#l           Valiable information               l
#l----------------------------------------------l

# Input file name
iname=str(args.i)
print("Input file name:"+str(iname)+"_?.pkl")
print("It can be changed by [-i "+iname+"]")
print("")

# Number of model
namen=int(str(args.n))

# Output file name
oname=str(args.o)
print("Output file name:"+str(oname)+".csv")
print("It can be changed by [-o "+oname+"]")
print("")

# Igunore residue
dres=str(args.d)
ures=dres.split("_")
print("Igunore residues:"+str(ures))
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
    with open(cpass+"/features.pkl", 'rb') as cn:
        cnn=pickle.load(cn)
        sequ=cnn['residue_index']
        rein=np.array(list(str(cnn['sequence'][0])[2:-1]))
        twoda=np.column_stack((sequ+1,rein))
        cn.close()
        pass
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
    print(eres)
    print()
    eres_f = eres.astype(np.float32)
    plddts=np.empty(namen)
    for i in range(namen):
        www=str(int(float(i)+1))
        with open(cpass+"/"+iname+"_"+www+".pkl", 'rb') as c1:
            c1n=pickle.load(c1)
            c1n_f = c1n['plddt'].astype(np.float32)
            na_mul=c1n_f*eres_f
            plddts[i]=np.sum(na_mul)/np.sum(eres_f)
            print(na_mul)
            pld1=np.round(c1n['plddt'], decimals=2)
            twoda=np.column_stack((twoda,pld1))
            pd1=['{:.2f}'.format(n) for n in pld1.tolist()]
            c1.close()
            pass
        ll1=""
        with open(cpass+"/relaxed_model_"+www+".pdb", mode='r') as g1:
            for line in g1:
                if  line[0:4]=="ATOM":
                    pdd1=pd1[int(float(line[22:26].strip())-1)]
                    ll1+=line[0:60]+" "*(6-len(pdd1))+pdd1+line[66:]
                    pass
                else:
                    ll1+=line
                    pass
                pass
            g1.close()
            pass
        with open(cpass+"/relaxed_"+str(lpass[-1])+"_"+www+".pdb", mode='w') as i:
            i.write(ll1)
            i.close()
            pass
        print("Save:"+cpass+"/relaxed_"+str(lpass[-1])+"_"+www+".pdb")
        pass
    print("------------")
    print()
    print("numbers of residues: "+str(int(np.sum(eres_f)))+" (Ignore: "+str(ures)+")")
    print("Rank  :"+str(rankdata(list(plddts))))
    np.set_printoptions(formatter={'float': '{:.1f}'.format})
    twoda=np.column_stack((twoda,eres))
    print("plddts:"+str(plddts))
    np.savetxt(cpass+'/'+oname+'.csv',twoda,delimiter=',',fmt="%s")

print("------------")
read()
print("------------")
now = time.ctime()
cnvtime = time.strptime(now)
print(time.strftime("%Y/%m/%d %H:%M", cnvtime))
print("Finish:"+str(sys.argv))
print("------------")
