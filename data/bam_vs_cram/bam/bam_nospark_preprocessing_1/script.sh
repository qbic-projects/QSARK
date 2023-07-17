nextflow pull FriederikeHanssen/sarek -r bam_31

nextflow run FriederikeHanssen/sarek -r bam_31 -profile cfc \
--input /sfs/9/ws/iizha01-QSARK/input.csv \
--outdir ./results_firsttry \
-c ~/trace.config -c config.config \
-resume
