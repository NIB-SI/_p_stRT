#!/usr/bin/env python2
#~*~ coding:utf-8 ~*~
import glob, re, os
# assuming this script is within scripts subdirectory
pathTo1 = '../input' # change appropriately - where assembler files are
pathTo2 = '../output' # change appropriately - where assembler renamed files and log file will be
try:
    os.makedirs(pathTo2)
except OSError:
    pass
f0 = open(pathTo2+"/Genotype_fasta-headers-rename.log", "a+") # change appropriately how do you want to call log table
fname=glob.glob(pathTo1+"/assemblerNameGenotype*.fasta") # change appropriately pattern in your file naming
for i in range (len(fname)):
    vname=fname[i].split("/")
    vname=vname[len(vname)-1]
    idname=re.sub('\.fasta$', '', vname)
    f4 = open(fname[i], "r")
    f5 = open(pathTo2+"/"+str(idname) + "_renamed.fasta", "w")
    line4=f4.readline()
    count=1
    while line4:
        line4=str.rstrip(line4)
        if '>' in line4:
            line44=">"+str(idname)+"_"+str(count)
            f5.write("%s\n" % (line44))
            f0.write("%s\t" % (line4))
            f0.write("%s\n" % (line44))
            count=count+1
        else:
            f5.write("%s\n" % (line4))
        line4 = f4.readline()
    f4.close()
    f5.close()
f0.close()
