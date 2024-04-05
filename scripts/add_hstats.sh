#!/bin/bash
for i in {567..571} #all:567-571 and 575-584
do
    root -b -l -q "/home/magda/project/macros/postprocessHeatMaps.C(\"/scratch1/gccb/data/Jan2023Beam/results/${i}_sifi_b0_HEATMAPS_ONLY.root\", \"/scratch1/gccb/data/Jan2023Beam/root/preselectionTimes/run00${i}_single_nlc_spills.root\")"
done

# TString path = "/scratch1/gccb/data/Jan2023Beam/results/00579_sifi_testing_HEATMAPS_ONLY.root", TString stats = "/scratch1/gccb/data/Jan2023Beam/root/preselectionTimes/run00579_single_nlc_spills.root"
