#!/usr/bin/bash

export HGREF=/nfsmounts/igenomes/Homo_sapiens/GATK/GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta
## Parameters from this preprint: https://www.biorxiv.org/content/biorxiv/early/2019/05/02/625624/DC2/embed/media-2.pdf?download=true & https://sites.google.com/view/seqc2/home/benchmarking-examples

for type in SNV INDEL
do
       for mapper in bwa bwamem2 dragmap
       do

       som.py reference_files/somatic_ref/high-confidence_s${type}_in_HC_regions_v1.2.vcf.gz \
              ../pipeline_results/sarek_30/${mapper}/results/variant_calling/strelka/HCC1395T_vs_HCC1395N/HCC1395T_vs_HCC1395N.strelka.somatic_${type,,}s.vcf.gz \
              --output results/strelka_${type}_${mapper} \
              --reference /nfsmounts/igenomes/Homo_sapiens/GATK/GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta \
              --normalize-all \
              -f reference_files/somatic_ref/High-Confidence_Regions_v1.2.bed \
              --target-regions reference_files/S07604624_Padded_Agilent_SureSelectXT_allexons_V6_UTR.bed \
              --logfile results/strelka_${type}_${mapper}.log
       awk -v d="strelka_"$type"_"$mapper"" -F"," 'FNR==1{a="run"} FNR>1{a=d} {print $0",\t"a}' results/strelka_${type}_${mapper}.stats.csv > results/strelka_${type}_${mapper}_formatted.stats.csv

       som.py reference_files/somatic_ref/high-confidence_s${type}_in_HC_regions_v1.2.vcf.gz \
              ../pipeline_results/sarek_30/${mapper}/results/variant_calling/mutect2/HCC1395T_vs_HCC1395N/HCC1395T_vs_HCC1395N.mutect2.filtered.vcf.gz \
              --output results/mutect2_filtered_${type}_${mapper} \
              --reference /nfsmounts/igenomes/Homo_sapiens/GATK/GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta \
              --normalize-all \
              -f reference_files/somatic_ref/High-Confidence_Regions_v1.2.bed \
              --target-regions reference_files/S07604624_Padded_Agilent_SureSelectXT_allexons_V6_UTR.bed \
              --logfile results/mutect2_filtered_${type}_${mapper}.log
       awk -v d="mutect2_filtered_"$type"_"$mapper"" -F"," 'FNR==1{a="run"} FNR>1{a=d} {print $0",\t"a}' results/mutect2_filtered_${type}_${mapper}.stats.csv > results/mutect2_filtered_${type}_${mapper}_formatted.stats.csv

       som.py reference_files/somatic_ref/high-confidence_s${type}_in_HC_regions_v1.2.vcf.gz \
              ../pipeline_results/sarek_30/${mapper}/results/variant_calling/mutect2/HCC1395T_vs_HCC1395N/HCC1395T_vs_HCC1395N.mutect2.vcf.gz \
              --output results/mutect2_${type}_${mapper} \
              --reference /nfsmounts/igenomes/Homo_sapiens/GATK/GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta \
              --normalize-all \
              -f reference_files/somatic_ref/High-Confidence_Regions_v1.2.bed \
              --target-regions reference_files/S07604624_Padded_Agilent_SureSelectXT_allexons_V6_UTR.bed \
              --logfile results/mutect2_${type}_${mapper}.log
       awk -v d="mutect2_"$type"_"$mapper"" -F"," 'FNR==1{a="run"} FNR>1{a=d} {print $0",\t"a}' results/mutect2_${type}_${mapper}.stats.csv > results/mutect2_${type}_${mapper}_formatted.stats.csv

       som.py reference_files/somatic_ref/high-confidence_s${type}_in_HC_regions_v1.2.vcf.gz \
              ../pipeline_results/sarek_30/bwa/results/variant_calling/freebayes/HCC1395T_vs_HCC1395N/HCC1395T_vs_HCC1395N.freebayes.vcf.gz \
              --output results/freebayes_${type}_${mapper} \
              --reference /nfsmounts/igenomes/Homo_sapiens/GATK/GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta \
              --normalize-all \
              -f reference_files/somatic_ref/High-Confidence_Regions_v1.2.bed \
              --target-regions reference_files/S07604624_Padded_Agilent_SureSelectXT_allexons_V6_UTR.bed \
              --logfile results/freebayes_${type}_${mapper}.log
       awk -v d="freebayes_"$type"_"$mapper"" -F"," 'FNR==1{a="run"} FNR>1{a=d} {print $0",\t"a}' results/freebayes_${type}_${mapper}.stats.csv > results/freebayes_${type}_${mapper}_formatted.stats.csv

       done

done

awk '(NR == 1) || (FNR > 1)' results/*formatted.stats.csv > merged_30.csv