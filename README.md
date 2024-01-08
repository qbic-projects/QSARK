# Overview

This repository contains the run scripts, relevant result files, and plotting scripts for the benchmarking of nf-ore/sarek 3.x.

There are two folders: `data` and `code`.

## Data

The `data` folder contains subfolders for each experiment type, scripts to run the experiment, additional configuration files, and input files needed to plot the data:

```bash

├── bam_vs_cram
│   ├── bam
│   └── cram
├── functional
│   ├── germline
│   └── somatic
├── intervals
│   ├── fastp0_intervals20
│   ├── fastp12
│   ├── fastp12_intervals10
│   ├── fastp16_intervals0
│   ├── fastp4_intervals78
│   ├── fastp8_intervals40
│   └── intervals10
├── pcawg_cnvs
│   ├── DO44888
│   ├── DO44889
│   ├── DO44890
│   ├── DO44919
│   └── DO44930
└── readme.md
```

> _Disclaimer: The CNV calls for each patient are not uploaded here. The original PCAWG results can be obtained from [https://dcc.icgc.org/](https://dcc.icgc.org/)_
.


## Code

The `code` folder contains Rmarkdown scripts, auxiliary files, a yaml file to create the conda environment to run the scripts, and a `results` folder, which contains the actual output plots. The `functions.R` file contains various functions used in the scripts.:

```bash
├── QBiCLogo.png
├── Rmarkdown_BAM_vs_CRAM.Rmd
├── Rmarkdown_BAM_vs_CRAM.html
├── Rmarkdown_cnvvalidation.Rmd
├── Rmarkdown_cnvvalidation.html
├── Rmarkdown_dataflow.Rmd
├── Rmarkdown_dataflow.html
├── Rmarkdown_functional_benchmark.Rmd
├── Rmarkdown_functional_benchmark.html
├── environment.yml
├── functions.R
├── qbic-style.css
└── results
    ├── bam_vs_cram
    ├── cnvs
    ├── dataflow
    ├── functional
    └── images
```




