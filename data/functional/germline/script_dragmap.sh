nextflow pull nf-core/sarek -r 3.1.2

nextflow run nf-core/sarek -r 3.1.2 -profile cfc \
--input /sfs/9/ws/iizha01-QSARK/functional/input.csv \
-c ~/trace.config \
--tools "deepvariant,freebayes,mpileup,strelka,manta,haplotypecaller" \
--aligner "dragmap" \
--skip_tools "baserecalibrator" \
--nucleotides_per_second 200000 \
-resume
