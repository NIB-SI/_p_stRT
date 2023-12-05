#!/usr/bin/env python
# zagor
from __future__ import print_function
from subprocess import check_call
import os
import treelib
import errno
import glob
from colorama import Fore, Back, Style, init
from termcolor import colored, cprint
import numbers

# ref: https://www.geeksforgeeks.org/print-colors-python-terminal/
def prRed(skk):
    print(f"\033[91m {skk}\033[00m") 
def prGreen(skk):
    print(f"\033[92m {skk}\033[00m") 
def prYellow(skk):
    print(f"\033[93m {skk}\033[00m") 
def prLightPurple(skk):
    print(f"\033[94m {skk}\033[00m") 
def prPurple(skk):
    print(f"\033[95m {skk}\033[00m") 
def prCyan(skk):
    print(f"\033[96m {skk}\033[00m") 
def prLightGray(skk):
    print(f"\033[97m {skk}\033[00m") 
def prBlack(skk):
    print(f"\033[98m {skk}\033[00m")
def prBlue(skk):
    print(f"\033[94m {skk}\033[00m")

# from BioPython
def read_fasta(fp):
    name, seq = None, []
    for line in fp:
        line = line.rstrip()
        if line.startswith(">"):
            if name: yield (name, ''.join(seq))
            name, seq = line, []
        else:
            seq.append(line)
    if name: yield (name, ''.join(seq))

prYellow("\nHello!\n")
# print(Style.RESET_ALL)


prRed ("Did you generate the input data files?\n")
prGreen ("1. Yes")
prPurple ("2. No")
prBlue ("3. I do not know")

is_valid=0
while not is_valid :
    try :
            choice = int ( raw_input('\nEnter your choice [1-3] : ') )
            is_valid = 1
    except ValueError, e :
            print ("'%s' is not a valid choice (integer)." % e.args[0].split(": ")[1])
if choice == 3:
    print("\nSee:")
    prBlue("####  " * 5)
    os.system('tree ./testSet/')
    prBlue("####  " * 5)
    print("\n")
elif choice == 2:
    prCyan("\n Well it's all right, riding around in the breeze\n Well it's all right, if you live the life you please\n Well it's all right, doing the best you can\n Well it's all right, as long as you lend a hand\n")
elif choice == 1:
    # Open a file
    print("\n")
    is_valid=False
    while not is_valid:
        fpath1 = str(raw_input('Path to your reference .fasta file: '))
        cnt = 0
        with open(fpath1) as fp:
            for name, seq in read_fasta(fp):
                # print(name)
                cnt = cnt + 1
        if cnt >=1:
            is_valid=True
    try:
        open(fpath1)
    except IOError as e:
        prPurple(os.strerror(e.errno))
    else:
        is_valid=False
        while not is_valid: 
            fpath2 = str(raw_input('Path to your folder containing files with written gene IDs of interest, line by line: '))
            is_valid = os.path.isdir(fpath2)
        res = str(raw_input('Where would you like to store your results (directory name, e.g. here/and/there)?: '))
        is_valid=False
        while not is_valid: 
            woa = str(raw_input('Preferred width of the alignment is (full or integer): '))
            if (isinstance(woa.isdigit(),numbers.Integral) or (woa == "full")):
                is_valid=True
            else:
                print("Enter integer number or full for unlimited!")
        is_valid=False
        while not is_valid: 
            cons = str(raw_input('Consensus on or off: '))
            if cons in ("on", "off"):
                is_valid=True
            else:
                print("Type on or off")
        print (15 * '-')
        print ("Choose aligner:")
        print ("1. ClustalOmega")
        print ("2. Muscle")
        print ("3. MAFFT")
        print (15 * '-')
        is_valid=0
        while not is_valid :
                try :
                        choice2 = int ( raw_input('Enter your choice [1-3] : ') )
                        is_valid = 1
                except ValueError, e :
                        print ("'%s' is not a valid choice (integer)." % e.args[0].split(": ")[1])
        if choice2 == 1:
            is_valid=False
            while not is_valid: 
                nHMM = str(raw_input('Number of iterations: '))
                if (isinstance(nHMM.isdigit(),numbers.Integral)):
                    is_valid=True
                else:
                    print("Enter integer number!")
            is_valid=False
            while not is_valid: 
                nCPU = str(raw_input('Number of CPU (up to 30): '))
                if (isinstance(nCPU.isdigit(),numbers.Integral)):
                    is_valid=True
                else:
                    print("Enter integer number!")
            # is_valid=False
            # while not is_valid: 
                # nRAM = str(raw_input('RAM MB (e.g. 8000): '))
                # if (isinstance(nRAM.isdigit(),numbers.Integral)):
                    # is_valid=True
                # else:
                    # print("Enter integer number!")
        elif choice2 == 2:
            is_valid=False
            while not is_valid: 
                niter = str(raw_input('Number of iterations: '))
                if (isinstance(niter.isdigit(),numbers.Integral)):
                    is_valid=True
                else:
                    print("Enter integer number!")
            is_valid=False
            while not is_valid: 
                nH = str(raw_input('Maximum time to iterate in hours: '))
                if (isinstance(nH.isdigit(),numbers.Integral)):
                    is_valid=True
                else:
                    print("Enter integer number!")
        elif choice2 == 3:
            is_valid=False
            while not is_valid: 
                niter = str(raw_input('Number of iterations: '))
                if (isinstance(niter.isdigit(),numbers.Integral)):
                    is_valid=True
                else:
                    print("Enter integer number!")
            is_valid=False
            while not is_valid: 
                nCPU = str(raw_input('Number of CPU (up to 30): '))
                if (isinstance(nCPU.isdigit(),numbers.Integral)):
                    is_valid=True
                else:
                    print("Enter integer number!")
        else:
                print ("Invalid number. Try again...")
        fpath1 = fpath1.strip()
        fpath2 = fpath2.strip()
        res = res.strip()
        if res[-1] == "/":
            res = res[0:-1]
        fname1 = fpath1.split('/')
        if (fname1[(len(fname1)-1)]):
            temp = fname1[(len(fname1)-1)]
        else:
            temp = fname1[(len(fname1)-2)]
        fname1 = temp
        # fname1 = temp.split('.')[0]
        # print(fname1)
        fname2 = fpath2.split('/')
        if (fname2[(len(fname2)-1)]):
            temp = fname2[(len(fname2)-1)]
        else:
            temp = fname2[(len(fname2)-2)]
        fname2 = temp
        # print(fname2)
        # print(fname2)
        # print(res)
        tmp = (os.listdir(fpath2))
        # print(tmp)
        # os.makedirs(res)
        print("\n")
        fnamej = fpath1.split('/')
        temp = fnamej[(len(fnamej)-1)]
        fnamej=temp
        fpathN = os.path.splitext(os.path.basename(fpath2))[0]
        for i in tmp:
            fnamei = i.split('/')
            temp = fnamei[(len(fnamei)-1)]
            fnamei=temp
            # fnamei = temp.split('.txt')[0]
            # print(fnamei)
            # check_call(["./print.sh", res+"_"+woa+"_"+cons+"_"+nHMM+"_"+fnamej, fnamei, fpath1, fpath2+"/", fnamej, woa, cons, nHMM, nCPU, nRAM ])
            if choice2 == 1:
                check_call(["./print.sh", res+"_ClustalOmega"+"_"+fpathN+"_width-"+woa+"_consensus-"+cons+"_"+nHMM+"-iterations_"+nCPU+"-threads", fnamei, fpath1, fpath2+"/", fnamej, woa, cons, nHMM, nCPU, "clustalo" ])
            elif choice2 == 2:
                check_call(["./print.sh", res+"_Muscle"+"_"+fpathN+"_width-"+woa+"_consensus-"+cons+"_"+niter+"-iterations_"+nH+"-hours", fnamei, fpath1, fpath2+"/", fnamej, woa, cons, niter, nH, "muscle" ])
            elif choice2 == 3:
                check_call(["./print.sh", res+"_MAFFT"+"_"+fpathN+"_width-"+woa+"_consensus-"+cons+"_"+niter+"-iterations__"+nCPU+"-threads", fnamei, fpath1, fpath2+"/", fnamej, woa, cons, niter, nCPU, "mafft" ])
            else:
                    print ("Invalid number. Try again...")
    prYellow("\nFinito :)\n")
    # print(Style.RESET_ALL) 
else:
        print ("Invalid number. Reenter your choice")



