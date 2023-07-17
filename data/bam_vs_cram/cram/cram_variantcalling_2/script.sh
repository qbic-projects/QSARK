nextflow pull nf-core/sarek -r 3.1.1

nextflow run nf-core/sarek -r 3.1.1 -profile cfc \
--input /sfs/9/ws/iizha01-QSARK/bam_vs_cram/cram/non_spark/result_thirdtry/csv/recalibrated.csv \
-c ~/trace.config -c config.config \
--tools deepvariant,haplotypecaller,mutect2,strelka,freebayes,ascat,controlfreec,cnvkit,manta,tiddit,msisensorpro \
--outdir ./results_vc_secondtry \
--step variant_calling \
-resume
