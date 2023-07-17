nextflow pull nf-core/sarek -r 3.1.1

nextflow run nf-core/sarek -r 3.1.1 -profile cfc \
--input /sfs/9/ws/iizha01-QSARK/input.csv \
--outdir ./results_thirdtry \
-c ~/trace.config -c config.config \
--use_gatk_spark "markduplicates" \
--skip_tools "fastqc,samtools,mosdepth,baserecalibrator" \
--aligner bwa-mem2 \
-resume
