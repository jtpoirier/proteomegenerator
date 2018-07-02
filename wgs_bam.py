WD="/data/kentsis/testing_scratchpad/rpe_p5"
workdir: WD
BWA='bwa'
BWA_VERSION='0.7.15-r1140'
SAMBLASTER='samblaster'
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
#FASTQSAMPLE='RPE-GFP-PGBD5'
UNITS=list(config['units'].keys())
#ALIGNMENTS=['aligned','realigned.alignments_merged']
NUM_TARGETS_RANGE=list(range(1,len(open(config['targets'],'r').readlines())+1))
NUM_VARIANT_INTERVALS=50 #Broad standard for GRCh38

PAIRED_SAMPLE='/data/kentsis/testing_scratchpad/K052_3/out/K052-normal.realigned_ubam-merged_RGmerged_dedup_fixedtags_recal.bam'

include: "wgs_common"

rule all:
    input: "out/RPE_GFP_P5.aligned_ubam-merged_RG-merged_dedup_fixedtags_recal_HTC.vcf"

if not config['bamContainsRGs']:
    rule AddRGInfoIfMissing:
        input: lambda wildcards: config['bam_sample'][wildcards.sample]
        output: temp("out/{sample}_readyToRevert.bam")
        params: n="4", R="'span[hosts=1] rusage[mem=9]'", o="out/logs/add_RGinfo.out", eo="out/logs/add_RGinfo.err", J="add_RGinfo"
        shell: '{PICARD} AddOrReplaceReadGroups -Xmx32g I={input} O={output} TMP_DIR={TMP} RGLB=lib1 \
                  RGPL=illumina RGPU=unit1 RGSM=xeno_RPE_GFP_P5 VALIDATION_STRINGENCY=LENIENT'

rule RevertToUnmappedBAM:
    #input: lambda wildcards: config['samples'][wildcards.sample]['read_groups'][wildcards.readgroup]['files']
    input: "out/{sample}_readyToRevert.bam"
    output: "out/readgroups/{sample}/reverted/{readgroup}.bam"
    params: n="48", R="'span[hosts=1] rusage[mem=8]'", o="out/logs/revert_bam.out", eo="out/logs/revert_bam.err", J="revert_bam"
    shell: "{PICARD} RevertSam -Xmx360g I={input} O={WD}/out/readgroups/{wildcards.sample}/reverted TMP_DIR={TMP} \
              MAX_RECORDS_IN_RAM=90000000 SANITIZE=true ATTRIBUTE_TO_CLEAR=XT \
              ATTRIBUTE_TO_CLEAR=XN ATTRIBUTE_TO_CLEAR=OC ATTRIBUTE_TO_CLEAR=OP \
              OUTPUT_BY_READGROUP=true VALIDATION_STRINGENCY=LENIENT"

#rule RevertToUnmappedBAM:
#    input: lambda wildcards: config['samples'][wildcards.sample]['read_groups'][wildcards.readgroup]['files']
#    output: "out/readgroups/{sample}/reverted/{readgroup}.bam"
#    params: n="48", R="'span[hosts=1] rusage[mem=8]'", o="out/logs/revert_bam.out", eo="out/logs/revert_bam.err", J="revert_bam"
#    shell: "{PICARD} RevertSam -Xmx360g I={input} O={WD}/out/readgroups/{wildcards.sample}/reverted TMP_DIR={TMP} \
#              MAX_RECORDS_IN_RAM=90000000 SANITIZE=true ATTRIBUTE_TO_CLEAR=XT \
#              ATTRIBUTE_TO_CLEAR=XN ATTRIBUTE_TO_CLEAR=OC ATTRIBUTE_TO_CLEAR=OP \
#              OUTPUT_BY_READGROUP=true"


rule RealignRevertedBAM:
    input: reverted="out/readgroups/{sample}/reverted/{readgroup}.bam", ref_idx=config['ref']+".fai", ref_dict=config['dict']
    #output: fq1=temp("out/readgroups/{sample}/{readgroup}.end1.fq"),fq2=temp("out/readgroups/{sample}/{readgroup}.end2.fq")
    output: "out/readgroups/{sample}/{readgroup}.aligned_ubam-merged.bam"
    params: n="48", R="'span[hosts=1] rusage[mem=8]'", o="out/logs/recreate_fq.out", eo="out/logs/recreate_fq.err", J="recreate_fq"
    benchmark: "out/benchmarks/realignment.txt"
    shell: "source activate wgs; \
            {PICARD} SamToFastq -Xmx360g I={input.reverted} FASTQ=/dev/stdout \
              CLIPPING_ATTRIBUTE=XT CLIPPING_ACTION=2 INTERLEAVE=true NON_PF=true \
              MAX_RECORDS_IN_RAM=90000000 TMP_DIR={TMP} | \
            {BWA} mem -M -t {params.n} -p {config[ref]} /dev/stdin | \
            {PICARD} MergeBamAlignment -Xmx360g ALIGNED_BAM=/dev/stdin UNMAPPED_BAM={input.reverted} \
              OUTPUT={output} R={config[ref]} CREATE_INDEX=true ADD_MATE_CIGAR=true \
              CLIP_ADAPTERS=false CLIP_OVERLAPPING_READS=true \
              INCLUDE_SECONDARY_ALIGNMENTS=true MAX_INSERTIONS_OR_DELETIONS=-1 \
              PRIMARY_ALIGNMENT_STRATEGY=MostDistant ATTRIBUTES_TO_RETAIN=XS \
              MAX_RECORDS_IN_RAM=90000000 TMP_DIR={TMP}"
