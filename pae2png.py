# Import
import os
import sys
import pickle
import csv
import sys
import time
import numpy as np
import argparse
import matplotlib.pyplot as plt
import seaborn as sns
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
print("The program extract pae values from [result_model_?_ptm.pkl] to [.png] ")
print("Thanks Dr. Yoshitaka Moriwaki (@Ag_smith) for initially showing how to export predicted aligned error in alphafold2.")
print("You need to change model_names to model_?_ptm from model_? at Docker/run_docker.py")
print("for Alphafold2")
print("------------")
# Header

# Argument & Help
parser = argparse.ArgumentParser(description="The program extract plddt from [result_model.pkl] to [.csv]")
parser.add_argument('-i', default='result_model', help='input XXXXX if [XXXXX_Y.pkl]; ex) result_model (def.)')
parser.add_argument('-n', default='5', help='input Y if you want open [XXXXX_1.pkl,XXXXX_2.pkl ... ,XXXXX_5.pkl]; ex) 5 (def.)')
parser.add_argument('-o', default='p2c', help='input XXXXX if you want to name the output [XXXXX.csv]; ex) plddt (def.)')
parser.add_argument('-c', default='Greens_r', help='color of heat map (PAE; ex) Greans_r (def.)')
parser.add_argument('-s', default='None', help="plt.show was skiped when -s 'Skip'")
parser.add_argument('-m', default=30, help="vmax of sns.heatmap ex) -m 30 (def.)")
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


# Color
cname=str(args.c)
print("Output file color:"+str(cname))
print("It can be changed by [-c "+cname+"]")
print("")


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
        eres=np.array(list(str(cnn['sequence'][0])[2:-1]))
        cn.close()
        pass
    np.set_printoptions(formatter={'float': '{:.1f}'.format})

    
    # Repeat for the caliculation of all models
    for i in range(namen):
        www=str(int(float(i)+1))
        # Read pLDDT from result_model_?_ptm.pkl
        with open(cpass+"/"+iname+"_"+www+"_ptm.pkl", 'rb') as c1:
            c1n=pickle.load(c1)
            # Extract PAE
            pae = c1n['predicted_aligned_error'].astype(np.float32)
            np.savetxt(cpass+'/'+oname+'_pae_'+www+'.csv',pae,delimiter=',',fmt="%s")
            print("Display:"+cpass+'/'+oname+'_pae_'+www+'.png')
            plt.figure(figsize=(15, 12))
            ax = sns.heatmap(pae, linewidth=0, cmap=cname, square=True, xticklabels=50, yticklabels=50, vmin=0, vmax=str(args.m))
            plt.savefig(cpass+'/'+oname+'_pae_'+www+'.png')
            if str(args.s)=='None':
                plt.show()
                pass
            # Extract PAE
            pass
        # Read pLDDT from result_model_?_ptm.pkl
        pass
    pass


#l-----------------------------------------------l
#l                   Execution                   l
#l-----------------------------------------------l

print("------------")
read()
print("------------")
print("Finish:"+str(sys.argv))
print("------------")
