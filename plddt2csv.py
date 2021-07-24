# Import
import os
import sys
import pickle
import csv
import sys
import time
import numpy as np
import argparse
# Import

# Header
now = time.ctime()
cnvtime = time.strptime(now)
print("------------")
print(time.strftime("%Y/%m/%d %H:%M", cnvtime))
print("Start:"+str(sys.argv))
print("Ver.1_(2021Jul.)")
print("by Hiroki Onoda, YCU")
print("------------")
print("The program extract plddt from [result_model.pkl] to [.csv] ")
print("for Alphafold")
print("------------")
# Header

# Argument & Help
parser = argparse.ArgumentParser(description="The program extract plddt from [result_model.pkl] to [.csv]")
parser.add_argument('-i', default='result_model', help='input XXXXX if [XXXXX_Y.pkl]; ex) result_model (def.)')
parser.add_argument('-n', default='5', help='input Y if you want open [XXXXX_1.pkl,XXXXX_2.pkl ... ,XXXXX_5.pkl]; ex) 5 (def.)')
parser.add_argument('-o', default='plddt', help='input XXXXX if you want to name the output [XXXXX.csv]; ex) plddt (def.)')
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
        pass
    for i in range(namen):
        www=str(int(float(i)+1))
        with open(cpass+"/"+iname+"_"+www+".pkl", 'rb') as c1:
            c1n=pickle.load(c1)
            pld1=np.round(c1n['plddt'], decimals=2)
            twoda=np.column_stack((twoda,pld1))
            pd1=['{:.2f}'.format(n) for n in pld1.tolist()]
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
            pass
        with open(cpass+"/relaxed_"+str(lpass[-1])+"_"+www+".pdb", mode='w') as i:
            i.write(ll1)
            i.close()
            pass
        print("Save:"+cpass+"/relaxed_"+str(lpass[-1])+"_"+www+".pdb")
        pass
    print(twoda)
    np.savetxt(cpass+'/'+oname+'.csv',twoda,delimiter=',',fmt="%s")

print("------------")
read()
print("------------")
now = time.ctime()
cnvtime = time.strptime(now)
print(time.strftime("%Y/%m/%d %H:%M", cnvtime))
print("Finish:"+str(sys.argv))
print("------------")
