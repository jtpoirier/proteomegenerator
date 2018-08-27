# ProteomeGenerator 2.0
DNA sequencing integration to proteomegenerator

WIP - workflow outlined in CBSP poster

Addition of 2 snakemake files to original proteomegenerator, all following GATK Best Practice protocol:

1) Preprocessing

* Convert unaligned fastq reads to aligned+unmapped reads
* Base recalibration and sort sam to convert files into analysis-ready BAM files

2) GATK

* Somatic variant call to form VCF files
* Create custom .fasta and .gtf files by incorporating SNPs and Indels VCF files

Current Problems:

* R markdown not working properly
* How to detect translocations? (https://academic.oup.com/gigascience/article/6/10/gix091/4160384) seems promising
* How to patch reference genome?
  * Can use USCS tool if given chain file (this is difficult to generate)
  * Otherwise, wait for g2gtools or GATK


![alt text](https://github.com/jtpoirier/proteomegenerator/blob/cr/workflow.png)