nextflow pull nf-core/sarek -r 3.1.1

nextflow run nf-core/sarek -r 3.1.1 -profile cfc \
--input /sfs/9/ws/iizha01-QSARK/input.csv \
-c ~/trace.config -c config.config \
--split_fastq 100000000 \
--skip_tools "baserecalibrator" \
-resume
