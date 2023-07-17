nextflow pull nf-core/sarek -r 3.1.1

nextflow run nf-core/sarek -r 3.1.1 -profile test_full,cfc \
--input ../../input.csv \
--intervals /sfs/9/ws/iizha01-QSARK/functional/pipeline_results/somatic/data/S07604624_Padded_Agilent_SureSelectXT_allexons_V6_UTR.bed \
--trim_fastq \
--outdir ./results_trimmed \
--tools "strelka,manta,freebayes,mutect2" \
-c ~/trace.config \
-resume
