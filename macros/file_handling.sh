#!/bin/bash

input_list=$1
sys_list=$2

for i in `cat ${input_list}`; do echo $i; tar -xzf $i; done

#fix some files that may have corrupted names
rename '_mc15_13TeV:' '_' flav_Akt4EMTo*mc15_13*root

echo "moving files in corresponding dirs..."
for i in `cat ${sys_list}`; do echo $i; mkdir $i; mv flav_Akt4EMToTR*${i}* $i/; done; 
mv flav_Akt4EMTo_* Nominal

echo "making lists..."
for i in `cat ${sys_list}`; do echo $i; ls $i/* > ${i}.list; done; 
#remove the first dir name
for i in `cat ${sys_list}`; do echo $i; sed -i "s/${i}\///g" ${i}.list; done; 
#remove syst name from the syst files 
for i in `cat ${sys_list}`; do echo $i; sed -i "s/TRK_${i}//g" ${i}.list; echo "nFiles "`wc -l ${i}.list`; done; 

echo "list comparison..."
for i in `cat ${sys_list}`; do 
    echo $i
    echo "nFiles "`wc -l ${i}.list`
    echo "Nominal diff"
    diff Nominal.list ${i}.list
done