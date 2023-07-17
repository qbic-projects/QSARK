nextflow pull FriederikeHanssen/sarek -r bam_31

nextflow run FriederikeHanssen/sarek -r bam_31 -profile cfc \
--input /sfs/9/ws/iizha01-QSARK/input.csv \
--outdir ./results_secondtry \
-c ~/trace.config -c config.config \
--use_gatk_spark "markduplicates" \
--skip_tools "fastqc,samtools,mosdepth,baserecalibrator" \
--aligner bwa-mem2 \
-resume
