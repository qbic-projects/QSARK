export HGREF=/nfsmounts/igenomes/Homo_sapiens/GATK/GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta

for mapper in dragmap bwa bwamem2
do
    for ref in HG003 HG004 HG002
    do

        hap.py /sfs/9/ws/iizha01-QSARK/functional/happy_results/reference_files/${ref}_GRCh38_1_22_v4.2.1_benchmark.vcf.gz \
        /sfs/9/ws/iizha01-QSARK/functional/pipeline_results/germline/sarek_30/${mapper}/results/variant_calling/deepvariant/${ref}/${ref}.deepvariant.vcf.gz \
        -o results/${ref}_${mapper}_deepvariant \
        -V --engine=vcfeval --threads 3 --engine-vcfeval-template grch38.sdf \
        -f /sfs/9/ws/iizha01-QSARK/functional/happy_results/reference_files/${ref}_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed \
        --logfile results/${ref}_${mapper}_deepvariant.log \
        --scratch-prefix .
        awk -v d="deepvariant_"$ref"_"$mapper"" -F"," 'FNR==1{a="run"} FNR>1{a=d} {print $0",\t"a}' results/${ref}_${mapper}_deepvariant.summary.csv > results/${ref}_${mapper}_deepvariant.formatted.summary.csv

        hap.py /sfs/9/ws/iizha01-QSARK/functional/happy_results/reference_files/${ref}_GRCh38_1_22_v4.2.1_benchmark.vcf.gz \
        /sfs/9/ws/iizha01-QSARK/functional/pipeline_results/germline/sarek_30/${mapper}/results/variant_calling/freebayes/${ref}/${ref}.freebayes.vcf.gz \
        -o results/${ref}_${mapper}_freebayes \
        -V --engine=vcfeval --threads 3 --engine-vcfeval-template grch38.sdf \
        -f /sfs/9/ws/iizha01-QSARK/functional/happy_results/reference_files/${ref}_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed \
        --logfile results/${ref}_${mapper}_freebayes.log \
        --scratch-prefix .
        awk -v d="freebayes_"$ref"_"$mapper"" -F"," 'FNR==1{a="run"} FNR>1{a=d} {print $0",\t"a}' results/${ref}_${mapper}_freebayes.summary.csv > results/${ref}_${mapper}_freebayes.formatted.summary.csv

        hap.py /sfs/9/ws/iizha01-QSARK/functional/happy_results/reference_files/${ref}_GRCh38_1_22_v4.2.1_benchmark.vcf.gz \
        /sfs/9/ws/iizha01-QSARK/functional/pipeline_results/germline/sarek_30/${mapper}/results/variant_calling/haplotypecaller/${ref}/${ref}.haplotypecaller.vcf.gz \
        -o results/${ref}_${mapper}_haplotypecaller \
        -V --engine=vcfeval --threads 3 --engine-vcfeval-template grch38.sdf \
        -f /sfs/9/ws/iizha01-QSARK/functional/happy_results/reference_files/${ref}_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed \
        --logfile results/${ref}_${mapper}_haplotypecaller.log \
        --scratch-prefix .
        awk -v d="haplotypecaller_"$ref"_"$mapper"" -F"," 'FNR==1{a="run"} FNR>1{a=d} {print $0",\t"a}' results/${ref}_${mapper}_haplotypecaller.summary.csv > results/${ref}_${mapper}_haplotypecaller.formatted.summary.csv

        hap.py /sfs/9/ws/iizha01-QSARK/functional/happy_results/reference_files/${ref}_GRCh38_1_22_v4.2.1_benchmark.vcf.gz \
        /sfs/9/ws/iizha01-QSARK/functional/pipeline_results/germline/sarek_30/${mapper}/results/variant_calling/haplotypecaller/${ref}/${ref}.haplotypecaller.filtered.vcf.gz \
        -o results/${ref}_${mapper}_haplotypecaller_filtered \
        -V --engine=vcfeval --threads 3 --engine-vcfeval-template grch38.sdf \
        -f /sfs/9/ws/iizha01-QSARK/functional/happy_results/reference_files/${ref}_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed \
        --logfile results/${ref}_${mapper}_haplotypecaller_filtered.log \
        --scratch-prefix .
        awk -v d="haplotypecaller_filtered_"$ref"_"$mapper"" -F"," 'FNR==1{a="run"} FNR>1{a=d} {print $0",\t"a}' results/${ref}_${mapper}_haplotypecaller_filtered.summary.csv > results/${ref}_${mapper}_haplotypecaller.filtered.formatted.summary.csv

        hap.py /sfs/9/ws/iizha01-QSARK/functional/happy_results/reference_files/${ref}_GRCh38_1_22_v4.2.1_benchmark.vcf.gz \
        /sfs/9/ws/iizha01-QSARK/functional/pipeline_results/germline/sarek_30/${mapper}/results/variant_calling/strelka/${ref}/${ref}.strelka.variants.vcf.gz \
        -o results/${ref}_${mapper}_strelka \
        -V --engine=vcfeval --threads 3 --engine-vcfeval-template grch38.sdf \
        -f /sfs/9/ws/iizha01-QSARK/functional/happy_results/reference_files/${ref}_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed \
        --logfile results/${ref}_${mapper}_strelka.log \
        --scratch-prefix .
        awk -v d="strelka_"$ref"_"$mapper"" -F"," 'FNR==1{a="run"} FNR>1{a=d} {print $0",\t"a}' results/${ref}_${mapper}_strelka.summary.csv > results/${ref}_${mapper}_strelka.formatted.summary.csv
    done
done

awk '(NR == 1) || (FNR > 1)' results/*.formatted.summary.csv > merged_30_germline.csv