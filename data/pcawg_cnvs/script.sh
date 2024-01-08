nextflow pull nf-core/sarek -r 3.4.0

nextflow run nf-core/sarek -r 3.4.0 -profile cfc \
--input input.csv \
--outdir results \
--tools ascat,controlfreec,cnvkit \
--only_paired_variant_calling \
-resume -c ~/trace.config -c ../custom.config -with-tower


