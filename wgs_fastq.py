#WD="/data/kentsis/testing_scratchpad/RPE-GFP-PGBD5"
WD=config['project_directory']
workdir: WD
BWA='bwa'
BWA_VERSION='0.7.15-r1140'
SAMTOOLS='samtools'
JAVA='java'
PICARD='picard'
PICARDJAR='~/picard.jar'
TMP='/scratch/kwokn/'
GATK37JAR='~/miniconda2/envs/wgs/opt/gatk-3.7/GenomeAnalysisTK.jar'
GATK4JAR='~/miniconda2/envs/wgs/share/gatk4-4.0b2-0/gatk-package-4.beta.2-local.jar'
SAMBAMBA='~/sambamba_v0.6.6'
BEDTOOLS='bedtools'

snakemake.utils.makedirs({WD})
snakemake.utils.makedirs('out/benchmarks')
snakemake.utils.makedirs('out/all-merge')
snakemake.utils.makedirs('out/logs')

#SAMPLES=list(config['samples'].keys())
#BAMSAMPLES=list(config['bam_samples'].keys())
#BAMSAMPLE='RPE_GFP_control'
#BAMSAMPLE_INPUT='/data/kentsis/WGS/K052_HemePACT/readgroups/s_OA_amlcl_003_CL_bc03_Proj_5932_L000.bam'
#BAMSAMPLE_INPUT='/data/kentsis/WGS/K052_HemePACT/s_OA_amlcl_003_CL_bc03_Proj_5932_L000_mrg_cl_aln_srt_MD_IR_BR.bam'
#FASTQSAMPLE='RPE-GFP-PGBD5'
#UNITS=list(config['units'].keys())
#ALIGNMENTS=['aligned','realigned.alignments_merged']
NUM_TARGETS_RANGE=list(range(1,len(open(config['targets'],'r').readlines())+1))
NUM_VARIANT_INTERVALS=50 #Broad standard for GRCh38

PAIRED_SAMPLE='/data/kentsis/testing_scratchpad/K052_3/out/K052-normal.realigned_ubam-merged_RGmerged_dedup_fixedtags_recal.bam'

include: "wgs_common.py"

#rule all:
#    input: expand("out/RPE-GFP_control.realigned_ubam-merged_RGmerged_dedup_fixedtags_recal_HTC.vcf")
 #expand("out/recal/K052.realigned_ubam-merged_RGmerged_dedup_fixedtags_{group}_recal.bam", group=NUM_TARGETS_RANGE), \
	   #expand("out/intervals/{interval}-scattered.intervals", interval=[str(x).zfill(4) for x in range(NUM_VARIANT_INTERVALS)]), \
	   #expand("out/RPE-GFP_control.realigned_ubam-merged_RGmerged_dedup_fixedtags_recal_HTC.vcf")
	   #expand("out/K052.realigned_ubam-merged_RGmerged_dedup_fixedtags_recal_m2.vcf"), \
	   #expand("out/RPE-GFP-PGBD5.aligned-fromfq_RGmerged_dedup_sorted_recal_HTC.vcf")

 	   #expand("out/intervals/K052.realigned_ubam_merged.RGmerged.dedup.fixedtags.recal.{interval}.vcf.gz", interval=[str(x).zfill(4) for x in range(NUM_VARIANT_INTERVALS)])

snakemake.utils.makedirs('out/readgroups')

rule AlignFastqReadsByRG:
    input: fq_files=lambda wildcards: config["samples"][wildcards.sample]['read_groups'][wildcards.readgroup]['files'], ref_idx=config['ref']+".fai", ref_dict=config['dict']
    output: temp("out/readgroups/{sample}/{readgroup}.aligned.bam")
    benchmark: "out/benchmarks/{sample}.{readgroup}.align.txt"
    params: n="48", R="'span[hosts=1] rusage[mem=8]'", o="out/logs/bwa_fastq.out", eo="out/logs/bwa_fastq.err", J="bwa_align_{wildcards.unit}", rg_header= lambda wildcards: config["samples"][wildcards.sample]['read_groups'][wildcards.readgroup]['header']
    shell: "source activate wgs; {BWA} mem -M \
              -t {params.n} \
              -R '{params.rg_header}' \
              {config[ref]} \
              {input.fq_files} | \
            {SAMTOOLS} view -Shu -@ {params.n} - | \
            {SAMBAMBA} sort -n -m 370G --tmpdir {TMP} -t {params.n} -o {output} /dev/stdin"

rule GenerateRGWiseUbamsFromFq:
    input: lambda wildcards: config["samples"][wildcards.sample]['read_groups'][wildcards.readgroup]['files']
    output: temp("out/readgroups/{sample}/{readgroup}.unmapped.bam")
    params: n="64", R="'span[hosts=1] rusage[mem=4]'", o="out/logs/fastq2ubam.out", eo="out/logs/fastq2ubam.err", read_group_header= lambda wildcards: config['samples'][wildcards.sample]['read_groups'][wildcards.readgroup]['header'], J="fq2ubams"
    run:
        rg_list=params.read_group_header.split('\\t')
        print(rg_list)
        rg_dict={k:v for k,v in (x.split(':') for x in [y for y in rg_list if ':' in y])}
        tags=['ID','SM','LB','PU','PL','CN','PI','PG','PM','DS','DT']
        rg_dict={x:rg_dict[x] if x in rg_dict else 'null' for x in tags}
        shell("source activate wgs;" + \
            "{PICARD} FastqToSam -Xmx64g TMP_DIR={TMP} " + \
            "FASTQ={input[0]} FASTQ2={input[1]} OUTPUT={output} " + \
            "READ_GROUP_NAME={rg_dict[ID]} SAMPLE_NAME={rg_dict[SM]} " + \
            "LIBRARY_NAME={rg_dict[LB]} PLATFORM_UNIT={rg_dict[PU]} " + \
            "PLATFORM={rg_dict[PL]} SEQUENCING_CENTER={rg_dict[CN]} " + \
            "PREDICTED_INSERT_SIZE={rg_dict[PI]} PROGRAM_GROUP={rg_dict[PG]} " + \
            "PLATFORM_MODEL={rg_dict[PM]} DESCRIPTION={rg_dict[DS]} " + \
            "RUN_DATE={rg_dict[DT]}")

rule MergeBamsByRG:
    input: aligned=temp("out/readgroups/{sample}/{readgroup}.aligned.bam"),unmapped=temp("out/readgroups/{sample}/{readgroup}.unmapped.bam")
    output: temp("out/readgroups/{sample}/{readgroup}.aligned_ubam-merged.bam")
    params: n="36", R="'span[hosts=1] rusage[mem=11]'", o="out/logs/merge_byRG.out", eo="out/logs/merge_byRG.err", J="merge_byRG"
    shell: "source activate wgs; \
            {PICARD} MergeBamAlignment -Xmx360g ALIGNED_BAM={input.aligned} UNMAPPED_BAM={input.unmapped} \
             OUTPUT={output} R={config[ref]} CREATE_INDEX=true ADD_MATE_CIGAR=true \
             CLIP_ADAPTERS=false CLIP_OVERLAPPING_READS=true \
             INCLUDE_SECONDARY_ALIGNMENTS=true MAX_INSERTIONS_OR_DELETIONS=-1 \
             PRIMARY_ALIGNMENT_STRATEGY=MostDistant ATTRIBUTES_TO_RETAIN=XS \
             MAX_RECORDS_IN_RAM=90000000 TMP_DIR={TMP}"

