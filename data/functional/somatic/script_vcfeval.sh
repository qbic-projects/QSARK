#!/usr/bin/bash

results_folder="./results_vcfeval/"

# for mapper in bwa bwamem2 dragmap
#        do

#               for sample in EA FD NV
#               do

#               bcftools view -f 'PASS' \
#                      ../pipeline_results/somatic/sarek_30/${mapper}/results_trimmed/variant_calling/strelka/HCC1395T_${sample}_vs_HCC1395N_${sample}/HCC1395T_${sample}_vs_HCC1395N_${sample}.strelka.somatic_snvs.vcf.gz \
#                      -o ${results_folder}${mapper}_PASS_HCC1395T_${sample}_vs_HCC1395N_${sample}.strelka.somatic_snvs.vcf

#               bgzip ${results_folder}${mapper}_PASS_HCC1395T_${sample}_vs_HCC1395N_${sample}.strelka.somatic_snvs.vcf
#               tabix ${results_folder}${mapper}_PASS_HCC1395T_${sample}_vs_HCC1395N_${sample}.strelka.somatic_snvs.vcf.gz

#               echo ""

#               bcftools view -f 'PASS' \
#                      ../pipeline_results/somatic/sarek_30/${mapper}/results_trimmed/variant_calling/strelka/HCC1395T_${sample}_vs_HCC1395N_${sample}/HCC1395T_${sample}_vs_HCC1395N_${sample}.strelka.somatic_indels.vcf.gz \
#                      -o ${results_folder}${mapper}_PASS_HCC1395T_${sample}_vs_HCC1395N_${sample}.strelka.somatic_indels.vcf

#               bgzip ${results_folder}${mapper}_PASS_HCC1395T_${sample}_vs_HCC1395N_${sample}.strelka.somatic_indels.vcf
#               tabix ${results_folder}${mapper}_PASS_HCC1395T_${sample}_vs_HCC1395N_${sample}.strelka.somatic_indels.vcf.gz

#               echo ""

#               bcftools view -f 'PASS' \
#                      ../pipeline_results/somatic/sarek_30/${mapper}/results_trimmed/variant_calling/mutect2/HCC1395T_${sample}_vs_HCC1395N_${sample}/HCC1395T_${sample}_vs_HCC1395N_${sample}.mutect2.filtered.vcf.gz \
#                      -o ${results_folder}${mapper}_PASS_HCC1395T_${sample}_vs_HCC1395N_${sample}.mutect2.filtered.vcf

#               bgzip ${results_folder}${mapper}_PASS_HCC1395T_${sample}_vs_HCC1395N_${sample}.mutect2.filtered.vcf
#               tabix ${results_folder}${mapper}_PASS_HCC1395T_${sample}_vs_HCC1395N_${sample}.mutect2.filtered.vcf.gz

#               echo ""

#               bcftools view -f 'PASS,.' \
#                      ../pipeline_results/somatic/sarek_30/${mapper}/results_trimmed/variant_calling/mutect2/HCC1395T_${sample}_vs_HCC1395N_${sample}/HCC1395T_${sample}_vs_HCC1395N_${sample}.mutect2.vcf.gz \
#                      -o ${results_folder}${mapper}_PASS_HCC1395T_${sample}_vs_HCC1395N_${sample}.mutect2.vcf

#               bgzip ${results_folder}${mapper}_PASS_HCC1395T_${sample}_vs_HCC1395N_${sample}.mutect2.vcf
#               tabix ${results_folder}${mapper}_PASS_HCC1395T_${sample}_vs_HCC1395N_${sample}.mutect2.vcf.gz

#               echo ""

#               bcftools view -f 'PASS,.' \
#                      ../pipeline_results/somatic/sarek_30/${mapper}/results_trimmed/variant_calling/freebayes/HCC1395T_${sample}_vs_HCC1395N_${sample}/HCC1395T_${sample}_vs_HCC1395N_${sample}.freebayes.vcf.gz \
#                      -o ${results_folder}${mapper}_PASS_HCC1395T_${sample}_vs_HCC1395N_${sample}.freebayes.vcf

#               bgzip ${results_folder}${mapper}_PASS_HCC1395T_${sample}_vs_HCC1395N_${sample}.freebayes.vcf
#               tabix ${results_folder}${mapper}_PASS_HCC1395T_${sample}_vs_HCC1395N_${sample}.freebayes.vcf.gz

#               echo ""

#        done
# done

for type in SNV INDEL
do
       for mapper in bwa bwamem2 dragmap
       do

              for sample in EA FD NV
              do

              echo ${type}"_"${mapper}"_"${sample}

              rtg vcfeval -c ${results_folder}${mapper}_PASS_HCC1395T_${sample}_vs_HCC1395N_${sample}.strelka.somatic_${type,,}s.vcf.gz \
                     -b ./reference_files/somatic_ref/high-confidence_s${type}_in_HC_regions_v1.2.vcf.gz \
                     -o ${results_folder}${mapper}_${sample}_strelka_${type,,}/ \
                     -t /sfs/9/ws/iizha01-QSARK/functional/happy_results/grch38.sdf \
                     -e ./reference_files/somatic_ref/High-Confidence_Regions_v1.2.bed \
                     --squash-ploidy --all-records --sample=ALT \
                     --bed-regions ./reference_files/S07604624_Padded_Agilent_SureSelectXT_allexons_V6_UTR.bed

              awk -v d="strelka_"$type"_"$mapper"_"$sample"" 'FNR==1{a="run"} FNR>1{a=d} {print $0"\t"a}' ${results_folder}${mapper}_${sample}_strelka_${type,,}/summary.txt > ${results_folder}/formatted_summary/strelka_${type}_${mapper}_${sample}.txt

              rtg vcfeval -c ${results_folder}${mapper}_PASS_HCC1395T_${sample}_vs_HCC1395N_${sample}.mutect2.filtered.vcf.gz \
                     -b ./reference_files/somatic_ref/high-confidence_s${type}_in_HC_regions_v1.2.vcf.gz \
                     -o ${results_folder}${mapper}_${sample}_mutect2_filtered_${type}/ \
                     -t /sfs/9/ws/iizha01-QSARK/functional/happy_results/grch38.sdf \
                     -e ./reference_files/somatic_ref/High-Confidence_Regions_v1.2.bed \
                     --squash-ploidy --all-records --sample=ALT \
                     --bed-regions ./reference_files/S07604624_Padded_Agilent_SureSelectXT_allexons_V6_UTR.bed

              awk -v d="mutect2_filtered_"$type"_"$mapper"_"$sample"" 'FNR==1{a="run"} FNR>1{a=d} {print $0"\t"a}' ${results_folder}${mapper}_${sample}_mutect2_filtered_${type}/summary.txt > ${results_folder}/formatted_summary/mutect2_filtered_${type}_${mapper}_${sample}.txt

              rtg vcfeval -c ${results_folder}${mapper}_PASS_HCC1395T_${sample}_vs_HCC1395N_${sample}.mutect2.vcf.gz \
                     -b ./reference_files/somatic_ref/high-confidence_s${type}_in_HC_regions_v1.2.vcf.gz \
                     -o ${results_folder}${mapper}_${sample}_mutect2_${type}/ \
                     -t /sfs/9/ws/iizha01-QSARK/functional/happy_results/grch38.sdf \
                     -e ./reference_files/somatic_ref/High-Confidence_Regions_v1.2.bed \
                     --squash-ploidy --all-records --sample=ALT \
                     --bed-regions ./reference_files/S07604624_Padded_Agilent_SureSelectXT_allexons_V6_UTR.bed

              awk -v d="mutect2_"$type"_"$mapper"_"$sample"" 'FNR==1{a="run"} FNR>1{a=d} {print $0"\t"a}' ${results_folder}${mapper}_${sample}_mutect2_${type}/summary.txt > ${results_folder}/formatted_summary/mutect2_${type}_${mapper}_${sample}.txt

              rtg vcfeval -c ${results_folder}${mapper}_PASS_HCC1395T_${sample}_vs_HCC1395N_${sample}.freebayes.vcf.gz \
                     -b ./reference_files/somatic_ref/high-confidence_s${type}_in_HC_regions_v1.2.vcf.gz \
                     -o ${results_folder}${mapper}_${sample}_freebayes_${type}/ \
                     -t /sfs/9/ws/iizha01-QSARK/functional/happy_results/grch38.sdf \
                     -e ./reference_files/somatic_ref/High-Confidence_Regions_v1.2.bed \
                     --squash-ploidy --all-records --sample=ALT \
                     --bed-regions ./reference_files/S07604624_Padded_Agilent_SureSelectXT_allexons_V6_UTR.bed

              awk -v d="freebayes_"$type"_"$mapper"_"$sample"" 'FNR==1{a="run"} FNR>1{a=d} {print $0"\t"a}' ${results_folder}${mapper}_${sample}_freebayes_${type}/summary.txt > ${results_folder}/formatted_summary/freebayes_${type}_${mapper}_${sample}.txt

              done
       done

done

awk '(NR == 1) || (FNR > 2)' ${results_folder}formatted_summary/*.txt > merged.txt