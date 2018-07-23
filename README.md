# ProteomeGenerator 2.0
DNA sequencing integration to proteomegenerator

WIP - workflow outlined in CBSP poster

Addition of 4 snakemake files to original proteomegenerator, all following GATK Best Practice protocol:

1) Preprocessing
*** Convert unaligned fastq reads to aligned+unmapped reads
2) Base Quality Score Recalibration
*** Base recalibration and sort sam to convert files into analysis-ready BAM files
3) Mutect2
*** Somatic variant call to form VCF files
4) Custom Reference
*** Create custom .fasta and .gtf files by incorporating SNPs and Indels VCF files


![alt text](https://github.com/kentsisresearchgroup/proteomegenerator2.0/blob/master/workflow.png)