#zagor
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
with open(f"{pathTo2}/Genotype_fasta-headers-rename.log", "a+") as f0:
    fname = glob.glob(f"{pathTo1}/assemblerNameGenotype*.fasta")
    for item in fname:
        vname = item.split("/")
        vname=vname[len(vname)-1]
        idname=re.sub('\.fasta$', '', vname)
        with open(item, "r") as f4:
            f5 = open(f"{pathTo2}/{str(idname)}_renamed.fasta", "w")
            line4=f4.readline()
            count=1
            while line4:
                line4=str.rstrip(line4)
                if '>' in line4:
                    line44 = f">{str(idname)}_{str(count)}"
                    f5.write("%s\n" % (line44))
                    f0.write("%s\t" % (line4))
                    f0.write("%s\n" % (line44))
                    count=count+1
                else:
                    f5.write("%s\n" % (line4))
                line4 = f4.readline()
        f5.close()
