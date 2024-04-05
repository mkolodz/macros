#!/bin/bash
path="/scratch1/gccb/data/Jan2023Beam/root/"
cd path
for i in {491..616}
do
    echo $path"run00${i} being converted..."
    
#     ./convert_raw_to_singles --config /home/lab/Desktop/DAQ/test_v_att_diff_bias_ig/config.ini -i /scratch1/gccb/data/TOFPET2/run00446 -o /scratch1/gccb/data/TOFPET2/root/raw00446_single.root --writeRoot
#     ./convert_raw_to_singles --config /home/lab/Desktop/DAQ/test_v_att_diff_bias_ig/config.ini -i $path"run00$i" -o $path"run00$i\_single_nlc.root" --writeRoot
done
