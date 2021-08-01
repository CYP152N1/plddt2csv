# Import
import os
import sys
import sys
import time
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
print("The program make run_docker_ptm.py to extract pae value")
print("for Alphafold2")
print("------------")
# Header

# Argument & Help
parser = argparse.ArgumentParser(description="The program make run_docker_ptm for caliculate pae value in alphafold2")
parser.add_argument('-i', default='run_docker', help='input .py file name ; ex) run_docker (def.)')
parser.add_argument('-o', default='run_docker_ptm', help='output .py file name ; ex) run_docker_ptm (def.)')
args = parser.parse_args()
# Argument & Help

#l----------------------------------------------l
#l           Valiable information               l
#l----------------------------------------------l

# Input file name
iname=str(args.i)
print("Input file name:"+str(iname)+"py")
print("It can be changed by [-i "+iname+"]")
print("")

# Output file name
oname=str(args.o)
print("Output file name:"+str(oname)+"py")
print("It can be changed by [-o "+oname+"]")
print("")

#l-----------------------------------------------l
#l           Valiable information                l
#l-----------------------------------------------l


def read():    
    cpass=os.getcwd()
    lpass=cpass.split("/")
    cowco=6
    liii=""
    with open(cpass+"/docker/run_docker.py", mode='r') as g:
        print("Read:"+cpass+"/docker/run_docker.py")
        for line in g:
            if cowco < 5.5:
                print(cowco)
                liii+=line[:12]+"_ptm"+line[12:]
                cowco+=1
                print(line)
                print(line[:12]+"_ptm"+line[12:])
                pass
            else:
                liii+=line
                if line[0:13]=="model_names =":
                    cowco+=-5
                    print("model_name was found")
                pass
            pass
        g.close()
        pass
    with open(cpass+"/docker/run_docker_ptm.py", mode='w') as i:
        i.write(liii)
        i.close()
        pass
    # Write pLDDT in the PDB_Bfactor
    print("Save:"+cpass+"/docker/run_docker_ptm.py")
    print("------------")
    
 
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
