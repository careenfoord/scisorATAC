#!/bin/bash

#local_path=$1

Rscript -e 'cat(system.file("bash", "v1.1a_exonInclusion_CTspecific_case_control.sh", package = "scisorATAC"))' > path_to_casesVcontrols

casesVcontrols_path=$(<path_to_casesVcontrols)

caseList=$1
controlList=$2
chrom_file=$3
numThreads=$4
annotation_path=$5

#find "$(cd ..; pwd)" -name "other-scripts" | awk '{print $1"/"}' > path_to_scripts
Rscript -e 'cat(system.file("bash", "other-scripts", package = "scisorATAC"))' > path_to_scripts
scripts_path=$(<path_to_scripts)

ci_low=$6
ci_upper=$7
min_reads=$8
zipping_function=${9}
OL_fraction=${10}
CellTypeFile=${11}
OutputDir=${12}

echo "$caseList" > caseList
#caseList_path=$(<caseList)
echo "$controlList" > controlList
#controlList_path=$(<controlList)


mkdir "$OutputDir"
#cd "$OutputDir"

echo "$(pwd)" > tmpdirs; time bash "$casesVcontrols_path" \
caseList controlList \
tmpdirs "$chrom_file" "$numThreads" "$annotation_path" \
"$scripts_path" "$ci_low" "$ci_upper" "$min_reads" "$zipping_function" "$OL_fraction" "$CellTypeFile" &> report

#Sig DPSI
zcat cases_vs_controls.counts.pVal.FDR.LOR.tab.gz | awk '$7 < .05' |awk '{split($1, a, "_"); dpsi = $2 / ($2 + $3) - ($4 / ($4 + $5)); \
print a[1]"_"a[2]"_"a[3]"_"a[5]"\t"a[4]"\t"dpsi"\t"$7}' > cases_vs_controls.counts.SIG_DPSI_FDR.tab

#All DPSI
zcat cases_vs_controls.counts.pVal.FDR.LOR.tab.gz | awk '{split($1, a, "_"); dpsi = $2 / ($2 + $3) - ($4 / ($4 + $5)); \
print a[1]"_"a[2]"_"a[3]"_"a[5]"\t"a[4]"\t"dpsi"\t"$7}' > cases_vs_controls.counts.All_DPSI_FDR.tab

#Large DPSI exons (0.5)
cat cases_vs_controls.counts.All_DPSI_FDR.tab | awk '($3 >.5 || $3 < -0.5){ print $1"\t"$2"\t"$3}' > LargeDPSI_exons
#Medium DPSI exons (0.25)
cat cases_vs_controls.counts.All_DPSI_FDR.tab | awk '($3 >.25 || $3 < -0.25){ print $1"\t"$2"\t"$3}' > MediumDPSI_exons

mv cases_vs_controls* "$OutputDir"
mv all.* "$OutputDir"
mv sampleBoth* "$OutputDir" 
mv *_exons "$OutputDir"
mv report "$OutputDir"
mv tmpdirs "$OutputDir"
mv *inc.exc.tab.gz "$OutputDir"
mv universe* "$OutputDir"

rm caseList
rm controlList
rm path_to_scripts
rm path_to_casesVcontrols
