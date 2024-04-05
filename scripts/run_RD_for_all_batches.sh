#!/bin/bash
batches_string="{"
for i in {598..598} #all:567-571 and 575-584
do
    for f in /scratch1/gccb/data/Jan2023Beam/results/${i}_sifi_b*.root; 
    do 
#         batches+=($f)
        batches_string+="\""
        batches_string+=$f
        batches_string+="\","
#         echo "adding element to batches"
    done
    batches_string=${batches_string::-1}
    batches_string+="}"
#     echo "${#batches[@]}"
    echo "$batches_string"
#     for i in ${#batches[@]}; 
#     for ((i=0;i<${#batches[@]};i++));
#     do 
#         echo ${batches[$i]}
#     done
#         root -b -l -q "/home/magda/project/macros/rawDisplayFibers_HeatMapsOnly.C({\"/scratch1/gccb/data/Jan2023Beam/results/00579_sifi_testing.root\"})"
        root -b -l -q "/home/magda/project/macros/rawDisplayFibers_HeatMapsOnly.C($batches_string)"
        batches_string="{"

done
