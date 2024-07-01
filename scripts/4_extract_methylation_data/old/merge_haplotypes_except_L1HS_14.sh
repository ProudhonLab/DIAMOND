#!/bin/bash

rootDir=${1}
primers=(2 15 4 8 9 7 6)

for primer in ${primers[@]}; do
	head -n 1 ${rootDir}/merged.L1HS_${primer}.methylation.csv > ${rootDir}/sorted.merged.L1HS_${primer}.methylation.csv
	sed 1d ${rootDir}/merged.L1HS_${primer}.methylation.csv | sort -V >> ${rootDir}/sorted.merged.L1HS_${primer}.methylation.csv
	cut -d ',' -f 1 ${rootDir}/sorted.merged.L1HS_${primer}.methylation.csv > ${rootDir}/diff_L1HS_${primer}
done

for primer1 in ${primers[@]}; do
	for primer2 in ${primers[@]}; do
		diff ${rootDir}/diff_L1HS_${primer1} ${rootDir}/diff_L1HS_${primer2} || exit 1
	done
done

cut -d ',' -f 1,2 ${rootDir}/sorted.merged.L1HS_${primers[0]}.methylation.csv > ${rootDir}/tmp_sample_biological-name

for primer in ${primers[@]}; do
	cut -d ',' -f 3- ${rootDir}/sorted.merged.L1HS_${primer}.methylation.csv > ${rootDir}/sorted.merged.L1HS_${primer}.only_haplotypes.methylation.csv
done

filesToConcatenate=()
for primer in ${primers[@]}; do
	filesToConcatenate+="${rootDir}/sorted.merged.L1HS_${primer}.only_haplotypes.methylation.csv "
done

paste -d ',' ${rootDir}/tmp_sample_biological-name ${filesToConcatenate} > ${rootDir}/sorted.merged.all_primers_except_L1HS_14.methylation.csv

# Cleanup
for primer in ${primers[@]}; do
	rm -f ${rootDir}/sorted.merged.L1HS_${primer}.only_haplotypes.methylation.csv
	rm -f ${rootDir}/diff_L1HS_${primer}
	rm -f ${rootDir}/sorted.merged.L1HS_${primer}.methylation.csv
done

rm -f ${rootDir}/tmp_sample_biological-name
