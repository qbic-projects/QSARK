nextflow pull nf-core/sarek -r 3.1.1

NXF_VER=22.10.2 nextflow run nf-core/sarek -r 3.1.1 -profile cfc \
--step 'prepare_recalibration' \
--input /sfs/9/ws/iizha01-QSARK/intervals/fastp_12_intervals_21/results/csv/markduplicates_no_table.csv \
-c ~/trace.config -c config.config \
--tools deepvariant,haplotypecaller,mutect2,strelka,freebayes,ascat,controlfreec,cnvkit,manta,tiddit,msisensorpro  \
--nucleotides_per_second 400000 \
-resume
