#!/bin/bash
# zagor


# pip2 install pyfaidx

# Num  Colour    #define         R G B
# 0    black     COLOR_BLACK     0,0,0
# 1    red       COLOR_RED       1,0,0
# 2    green     COLOR_GREEN     0,1,0
# 3    yellow    COLOR_YELLOW    1,1,0
# 4    blue      COLOR_BLUE      0,0,1
# 5    magenta   COLOR_MAGENTA   1,0,1
# 6    cyan      COLOR_CYAN      0,1,1
# 7    white     COLOR_WHITE     1,1,1

# write proper paths
myMuscle=/home/administrator/Software/muscle3.8.31_i86linux64
export PERL5LIB=/home/administrator/Software/mview-1.66/lib:$PERL5LIB
myMview=/home/administrator/Software/mview-1.66/bin/mview
myClustal=clustalw # CLUSTAL 2.1 Multiple Sequence Alignments
myPhylo=simple_phylogeny.pl # clustalw
myClustalo=clustalo # Clustal Omega - 1.2.1 (AndreaGiacomo)
myMAFFT=mafft # MAFFT v7.271 (2016/1/6) http://mafft.cbrc.jp/alignment/software/ MBE 30:772-780 (2013), NAR 30:3059-3066 (2002)

chmod 775 ../scripts/
chmod 755 ../input/

echo "output folder: ../output/"$1
echo "File with IDs: "$2
# echo $3
# echo $4
echo "Fasta file: "$5
# echo $6
# echo $7
# echo $8
# echo $9
# echo "${10}"
# echo $1_"${10}"


dos2unix $4$2

mkdir -p ../output/
mkdir -pv $1/

outfile="./"$1"/"$2"_myLog.txt"
date > $outfile
echo  2>&1 | tee -a $outfile

# echo $1
echo "$(tput setaf 3)####  ####  ####  ####  ####  ####  ####  ####  ####  ####  ####  ####  $(tput sgr 0)"
echo "$(tput setaf 3)Using: "$2" and "$5"$(tput sgr 0)"
echo "Using: "$4$2" and "$3 >> $outfile
echo  2>&1 | tee -a $outfile


f1="$2_fasta"
f2="$2_.err.txt"
f3="$2.fasta"
f4="$2_aligned"
f5="$2.aligned"
f6="$2.aligned.mview.html"


echo "$(tput setaf 2)get sequences$(tput sgr 0)"
echo "####  ########  ########  ########  ####" >> $outfile
echo "get sequences" >> $outfile
echo  2>&1 | tee -a $outfile

xargs faidx -d ' ' $3 < $4$2 > "./"$1"/"$f1 2> "./"$1"/"$f2;
fasta_formatter -i "./"$1"/"$f1 -o "./"$1"/"$f3 -w 0;
# fastx_collapser
rm "./"$1"/"$f1;

# [ -s  "./"$1"/"$f3 ] && echo "File not empty" || echo "File empty" 



if [ -s "./"$1"/"$f3 ]
then
    # By default, MUSCLE re-arranges sequences so that similar sequences are adjacent in the output file. (This is done by ordering sequences according to a prefix traversal of the guide tree). This makes the alignment easier to evaluate by eye. If you want to the sequences to be output in the same order as the input file, you can use the –stable option.
    # echo "$(tput setaf 2)muscle$(tput sgr 0)"
    # echo "$(tput setaf 2)clustalo$(tput sgr 0)"
    echo "####  ########  ########  ########  ####" >> $outfile
    # echo "muscle" >> $outfile
    # echo "clustalo" >> $outfile
    if [ "${10}" == "clustalo" ]
    then
        echo "$(tput setaf 2)ClustalOmega$(tput sgr 0)"
        $myClustalo --in "./"$1"/"$f3 --out "./"$1"/"$f5 --resno --threads $9 --iterations $8 --max-hmm-iterations $8 --percent-id --full --full-iter --clustering-out="./"$1"/"$f5"clucstalo.cluster" --force 2>&1 | tee -a $outfile;
        echo "ClustalOmega" >> $outfile
    elif [ "${10}" == "muscle" ]
    then
        echo "$(tput setaf 2)Muscle$(tput sgr 0)"
        $myMuscle -in "./"$1"/"$f3 -out "./"$1"/"$f5 -maxiters $8 -maxhours $9 2>&1 | tee -a $outfile;
        echo "Muscle" >> $outfile
    elif [ "${10}" == "mafft" ]
    then
        echo "$(tput setaf 2)MAFFT$(tput sgr 0)"
        $myMAFFT --maxiterate $8 --thread $9 "./"$1"/"$f3 > "./"$1"/"$f5
        echo "MAFFT" >> $outfile
    else
       echo "Error - none of the condition met."
    fi
    # $myMuscle -in "./"$1"/"$f3 -out "./"$1"/"$f5 -maxiters $8 -maxhours $9 2>&1 | tee -a $outfile;
    # $myClustalo --in "./"$1"/"$f3 --out "./"$1"/"$f5 --resno --threads $9 --iterations $8 --max-hmm-iterations $8 --pileup --percent-id --full --full-iter --clustering-out="./"$1"/"$f5"clucstalo.cluster" --force 2>&1 | tee -a $outfile;
    # $myMAFFT --maxiterate $8 --thread $9 "./"$1"/"$f3 > "./"$1"/"$f5
    xargs faidx -d ' ' "./"$1"/"$f5 < $4$2 > "./"$1"/"$f4; # reorder
    fasta_formatter -i "./"$1"/"$f4 -o "./"$1"/"$f5 -w 0;
    rm "./"$1"/"$f4;
    
    echo  2>&1 | tee -a $outfile
    echo "$(tput setaf 2)mview$(tput sgr 0)"
    echo "####  ########  ########  ########  ####" >> $outfile
    echo "mview" >> $outfile
    echo  2>&1 | tee -a $outfile
    
    $myMview -in fasta "./"$1"/"$f5 -consensus $7 \
    -reference 1 -pcid reference \
    -coloring identity -colormap CLUSTAL \
    -threshold 50 -html full -css on \
    -title "STRT2019" -width $6 -conservation on \
    -out mview > "./"$1"/"$f6 2>&1 | tee -a $outfile;
    
    echo "$(tput setaf 2)tree$(tput sgr 0)"
    echo "####  ########  ########  ########  ####" >> $outfile
    echo "tree" >> $outfile
    
    $myClustal $myPhylo -INFILE="./"$1"/"$f5 -clustering=nj -TREE -OUTORDER=INPUT -BOOTSTRAP=1000 -KIMURA -BOOTLABELS=branch 2>&1 | tee -a $outfile;
    
    # echo  2>&1 | tee -a $outfile
    
    rm "./"$1"/$2.aligned.fai"
    
    echo "$(tput setaf 1)Results at ./"$1"$(tput sgr 0)"
    echo "Results at ./"$1 >> $outfile
    echo  2>&1 | tee -a $outfile
else 
    echo $4$2 "not in" $3
fi



date >> $outfile

chmod 775 ../output/
mkdir -p ../output/$1
mv $1/* ../output/$1
rm -rf $1
