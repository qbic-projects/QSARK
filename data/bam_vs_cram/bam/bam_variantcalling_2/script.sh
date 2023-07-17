nextflow pull FriederikeHanssen/sarek -r bam_31

nextflow run FriederikeHanssen/sarek -r bam_31 -profile cfc \
--step 'variant_calling' \
--outdir './results_vc_secondtry' \
--tools 'deepvariant,haplotypecaller,mutect2,strelka,freebayes,ascat,controlfreec,cnvkit,manta,tiddit,msisensorpro' \
-c ~/trace.config -c config.config \
--input './results_thirdtry/csv/recalibrated.csv' \
-resume