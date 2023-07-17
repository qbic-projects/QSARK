nextflow pull nf-core/sarek -r 3.1.2

nextflow run nf-core/sarek -r 3.1.2 -profile cfc \
--input ./input.csv \
--outdir /sfs/9/ws/iizha01-QSARK/functional/pipeline_results/germline/sarek_30/bwamem2 \
-c ~/trace.config -c config.config \
--tools "deepvariant,freebayes,mpileup,strelka,manta,haplotypecaller" \
--aligner "bwa-mem2" \
--nucleotides_per_second 200000 \
-resume
