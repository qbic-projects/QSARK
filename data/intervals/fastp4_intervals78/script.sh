nextflow pull nf-core/sarek -r 3.1.1

nextflow run nf-core/sarek -r 3.1.1 -profile cfc \
--input /sfs/9/ws/iizha01-QSARK/input.csv \
-c ~/trace.config -c config.config \
--tools deepvariant,haplotypecaller,mutect2,strelka,freebayes,ascat,controlfreec,cnvkit,manta,tiddit,msisensorpro \
--nucleotides_per_second 10001 \
--split_fastq 10000000000 \
-resume
