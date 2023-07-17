nextflow pull nf-core/sarek -r 3.1.1

nextflow run nf-core/sarek -r 3.1.1 -profile cfc \
--input /sfs/9/ws/iizha01-QSARK/input.csv \
-c ~/trace.config -c config.config \
--tools deepvariant,haplotypecaller,mutect2,strelka,freebayes,ascat,controlfreec,cnvkit,manta,tiddit,msisensorpro \
--nucleotides_per_second 5000000 \
--split_fastq 100000000 \
-resume
