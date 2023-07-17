nextflow pull nf-core/sarek -r 3.1.2

nextflow run nf-core/sarek -r 3.1.2 -profile cfc \
--input ./input.csv \
-c ~/trace.config -c config.config \
--tools "deepvariant,freebayes,mpileup,strelka,manta,haplotypecaller" \
--nucleotides_per_second 200000 \
-resume
