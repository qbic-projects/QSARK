# CNV comparison

The WGS BAM files were downloaded from the ICGC portal as described in their documentation using the PDC client. Each patient was then processed with `nf-core/sarek` with the `script.sh` file.

In addition the CNV files for each patient were downloaded. The VCF files were converted to a table format with:

```bash
gatk VariantsToTable -V .somatic.cnv.vcf -O output.table
```

The CNV calls of all patients and all callers (nf-core/sarek (ASCAT, ControLFREEC, CNVKit), ICGC (DKFZ, Sanger))were then compared and plotted.