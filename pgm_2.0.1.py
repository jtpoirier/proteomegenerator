import glob,os,subprocess

SAMPLES= "3-27-18_test1".split() # should only ever list 1 sample

WD="/data/kentsis/testing_scratchpad/" + SAMPLES[0]
workdir: WD

SAMPLE_VCF = "/data/kentsis/testing_scratchpad/K052_3/out/s_OA_amlcl_003_CL.realigned_ubam_merged.RGmerged.dedup.fixedtags.recal.finished.vcf"
VCF_INDIVIDUAL = subprocess.check_output(['bcftools','query','-l',SAMPLE_VCF]).decode('ASCII').strip()

REF_FASTA='/data/kentsis/indexes/GRCh38/GRCh38.genome.fa'
REF_GTF='/data/kentsis/indexes/GRCh38/gencode.v20.annotation.gtf'
REF_FASTA_NAME=os.path.split(config['ref'])[1].strip('.fa')
print(REF_FASTA_NAME)
REF_DIR=os.path.split(config['ref'])[0]


# variables
PY = '/home/kwokn/miniconda2/bin/python'
UNIPROT="/home/kwokn/proteomegenerator/indexes/uniprot/UP000005640.fasta"
PGX="/home/kwokn/PGx"
RSCRIPT='/opt/common/CentOS_7/R/R-3.3.3/bin/Rscript'
SCRIPTS='/home/kwokn/scripts'
TMP = '/scratch/kwokn'
EMBOSS = "/home/kwokn/miniconda2/pkgs/emboss-6.5.7-2/bin"
MSFRAG = "/home/kwokn/MSFragger_20170103"
TD_DIR='~/miniconda2/opt/transdecoder'

#SAMPLE = os.path.split(WD)[1]
#SAMPLES = "Sample_AH_RNAseq_GFP-P5-Tumor_01".split()
#SPECTRA_SAMPLES = "b1906_293T_proteinID_01A".split()
SPECTRA_SAMPLES=[os.path.abspath(a) for a in glob.glob("/data/kentsis/proteomics/K052/*.mzXML")]
MS_DIR="/data/kentsis/proteomics/K052"
MS_RELATIVE_DIR=os.path.relpath(MS_DIR,WD)
MS_SAMPLES=['160116_K052_OffLRP_RP_f','20160205_PGM_K052_SCX_RP_']
MS_RUNS=[str(x).zfill(2) for x in list(range(1,24+1))]

#SAMPLE_RNASAMPLE_RUNS = expand("/data/kentsis/RNAseq/K0562/FCH9EFLADXX-HUMbghEAACRAAPEI-225_L{lane}_{run}.fq.gz",lane=[1,2],run=[1,2])
RNASAMPLE_DIR="/data/kentsis/RNAseq/K0562"
RNASAMPLE_SAMPLE="FCH9EFLADXX-HUMbghEAACRAAPEI-225"
RNASAMPLE_RUNS = [1,2]
RNASAMPLE_LANES = [1,2]

#MODELS = 'merged reference'.split()
MODELS = 'K052'.split()
ADBS = expand('adb{n}',n=range(1,8))
DBS = ADBS + 'merged reference uniprot Ecoli loki'.split()

REF_GENOME_STAR='/data/kentsis/indexes/GRCh38/STAR'
INDEX=REF_GENOME_STAR

if 'REF_GTF' in locals():
    ruleorder: STAR_GTF > STAR_denovo
    ruleorder: StringTie_GTF > StringTie_denovo
    ruleorder: UCSC_GTF > UCSC_denovo
    MODELS = 'sample reference'.split()
    #MODELS = MODELS + ['reference']
else:
    ruleorder: STAR_denovo > STAR_GTF
    ruleorder: StringTie_denovo > StringTie_GTF
    ruleorder: UCSC_denovo > UCSC_GTF
    MODELS = 'sample'

snakemake.utils.makedirs(WD)
snakemake.utils.makedirs('out/benchmarks')
snakemake.utils.makedirs('out/all-merge')
snakemake.utils.makedirs('out/logs')
snakemake.utils.makedirs('out/rnasample_runs')

rule all:
    input: expand("out/{model}.proteome/proteome.pickle",model=['sample','reference'])
	   #expand("out/{SAMPLE}.proteome/proteome.pickle",SAMPLE=SAMPLES), \
           #expand(MS_RELATIVE_DIR+'/{ms_sample}{run}.{SAMPLE}.db.pepXML',ms_sample=MS_SAMPLES,run=MS_RUNS,SAMPLE=SAMPLES)
        #expand("out/all-merge/{db}.db.fasta.pepdigest",db=DBS),

rule STAR_index:
    input: fasta = REF_FASTA
    output: INDEX
    benchmark: "out/benchmarks/index.txt"
    log: "out/logs/index.txt"
    params: n="12", R="'span[hosts=1] rusage[mem=15]'", J="index", o="out/logs/index.out", eo="out/logs/index.err"
    shell: "mkdir {output} ; \
            {config[STAR]} \
            --runThreadN {params.n} \
            --runMode genomeGenerate \
            --genomeDir {output} \
            --genomeFastaFiles {input.fasta} 2> {log}"

rule STAR_GTF:
    #input: expand(RNASAMPLE_DIR+"/"+RNASAMPLE_SAMPLE+"_L{lane}_{{run}}.fq.gz",lane=RNASAMPLE_LANES)
    input: sample=expand(RNASAMPLE_DIR+"/"+RNASAMPLE_SAMPLE+"_L{lane}_{{run}}.fq.gz",lane=RNASAMPLE_LANES), index=INDEX
    output: temp("out/rnasample_runs/run{run}.Aligned.sortedByCoord.out.bam")
    benchmark: "out/benchmarks/{run}.align.json"
    log: "out/logs/STAR_GTF.txt"
    params: n="12", R="'span[hosts=1] rusage[mem=20]'", J="align", o="out/logs/align.out", eo="out/logs/align.err"
    shell: "{config[STAR]} \
        --genomeDir {input.index} \
        --readFilesIn {input.sample} \
        --outFileNamePrefix out/rnasample_runs/run{wildcards.run}. \
        --outSAMattributes NH HI XS \
        --outSAMattrRGline ID:{wildcards.run} LB:1 PL:illumina PU:1 SM:{wildcards.run} \
        --runThreadN {params.n} \
        --outSAMtype BAM SortedByCoordinate \
        --clip3pAdapterSeq AGATCGGAAGAG \
        --readFilesCommand zcat \
        --twopassMode Basic \
        --outSAMstrandField intronMotif \
        --outFilterIntronMotifs None \
        --outReadsUnmapped None \
        --chimSegmentMin 15 \
        --chimJunctionOverhangMin 15 \
        --alignMatesGapMax 1000000 \
        --alignIntronMax 1000000 \
        --outFilterType Normal \
        --alignSJDBoverhangMin 1 \
        --alignSJoverhangMin 8 \
        --outFilterMismatchNmax 1 \
        --outSJfilterReads Unique \
        --outFilterMultimapNmax 10 \
        --sjdbOverhang 100 \
        --sjdbGTFfile {REF_GTF} \2 > {log}"

rule STAR_denovo:
    #input: expand(RNASAMPLE_DIR+"/"+RNASAMPLE_SAMPLE+"_L{lane}_{{run}}.fq.gz",lane=RNASAMPLE_LANES)
    input: sample=expand(RNASAMPLE_DIR+"/"+RNASAMPLE_SAMPLE+"_L{lane}_{{run}}.fq.gz",lane=RNASAMPLE_LANES), index=INDEX
    output: temp("out/rnasample_runs/run{run}.Aligned.sortedByCoord.out.bam")
    benchmark: "out/benchmarks/{run}.align.json"
    log: "out/logs/STAR_denovo.txt"
    params: n="12", R="'span[hosts=1] rusage[mem=20]'", J="align", o="out/logs/align.out", eo="out/logs/align.err"
    shell: "{config[STAR]} \
        --genomeDir {input.index} \
        --readFilesIn {input.sample} \
        --outFileNamePrefix out/rnasample_runs/run{wildcards.run}. \
        --outSAMattributes NH HI XS \
        --outSAMattrRGline ID:{wildcards.run} LB:1 PL:illumina PU:1 SM:{wildcards.run} \
        --runThreadN {params.n} \
        --outSAMtype BAM SortedByCoordinate \
        --clip3pAdapterSeq AGATCGGAAGAG \
        --readFilesCommand zcat \
        --twopassMode Basic \
        --outSAMstrandField intronMotif \
        --outFilterIntronMotifs None \
        --outReadsUnmapped None \
        --chimSegmentMin 15 \
        --chimJunctionOverhangMin 15 \
        --alignMatesGapMax 1000000 \
        --alignIntronMax 1000000 \
        --outFilterType Normal \
        --alignSJDBoverhangMin 1 \
        --alignSJoverhangMin 8 \
        --outFilterMismatchNmax 1 \
        --outSJfilterReads Unique \
        --outFilterMultimapNmax 10 \
        --sjdbOverhang 100 \2 > {log}"

rule filter:
    input: bam="out/rnasample_runs/run{run}.Aligned.sortedByCoord.out.bam"
    output: "out/rnasample_runs/run{run}.Aligned.trimmed.out.bam"
    log: "out/logs/rna_sample/run{run}.filter.txt"
    benchmark: "out/benchmarks/run{run}.filter.txt"
    params: n="1", R="'span[hosts=1] rusage[mem=10]'", J="filter", o="out/logs/filter.out", eo="out/logs/filter.err"
    shell: "{config[SAMTOOLS]} view -b -h -F 4 -F 8 -F 256 -F 512 -F 2048 -q 30 {input.bam} > {output} 2> {log}"

rule BuildBamIndex:
    input: "out/rnasample_runs/run{run}.Aligned.trimmed.out.bam"
    output: "out/rnasample_runs/run{run}.Aligned.trimmed.out.bai"
    log: "out/logs/rna_sample/run{run}.bamIndex.txt"
    benchmark: "out/benchmarks/run{run}.BuildBamIndex"
    params: n="1", R="'span[hosts=1] rusage[mem=10]'", \
            o="out/logs/out.buildbamindex", eo="out/logs/error.buildbamindex", \
            J="BuildBamIndex"
    shell: "java -Djava.io.tmpdir={TMP} -Xmx6g -jar {config[PICARD]} \
            BuildBamIndex \
            INPUT={input} 2> {log}"

rule StringTie_GTF:
    input: bam="out/rnasample_runs/run{run}.Aligned.trimmed.out.bam", bai="out/rnasample_runs/run{run}.Aligned.trimmed.out.bai"
    output: "out/rnasample_runs/run{run}.stringtie.gtf"
    log: "out/logs/rnasample_runs/run{run}.filterAndTrimBed.txt"
    benchmark: "out/benchmarks/run{run}.StringTie_GTF"
    params: n="6", R="'span[hosts=1] rusage[mem=20]'", \
            o="out/logs/out.stringtie", eo="out/logs/error.stringtie", \
	    J="StringTie"
    shell: "{config[STRINGTIE]} \
            -G {REF_GTF} \
            {input.bam} \
            -p {params.n} \
            -o {output} \
            -c 2.5 \
            -m 300 \
            -f .01 2> {log}"

rule StringTie_denovo:
    input: bam="out/rnasample_runs/run{run}.Aligned.trimmed.out.bam", bai="out/rnasample_runs/run{run}.Aligned.trimmed.out.bai"
    output: "out/rnasample_runs/run{run}.stringtie.gtf"
    log: "out/logs/rnasample_runs/run{run}.filterAndTrimBed.txt"
    benchmark: "out/benchmarks/run{run}.StringTie_denovo"
    params: n="6", R="'span[hosts=1] rusage[mem=20]'", \
            o="out/logs/out.stringtie", eo="out/logs/error.stringtie", \
        J="StringTie"
    shell: "{config[STRINGTIE]} \
            {input.bam} \
            -p {params.n} \
            -o {output} \
            -c 2.5 \
            -m 300 \
            -f .01 2> {log}"

snakemake.utils.makedirs('out/all-merge')

rule merge:
    input: expand("out/rnasample_runs/run{run}.stringtie.gtf",run=RNASAMPLE_RUNS)
    output: "out/all-merge/stringtie.merged.gtf"
    log: "out/logs/stringtie.merge.txt"    
    benchmark: "out/benchmarks/merge.txt"
    params: n="12", R="'span[ptile=72] rusage[mem=4]'", \
            o="out/logs/out.merge", eo="out/logs/error.merge", \
            J="merge"
    shell: "{config[STRINGTIE]} \
            --merge \
	    -o {output} \
            -p {params.n} \
            -c 2.5 \
            -m 300 \
            -T 1 \
            -f .01 \
            -i \
	    {input} 2> {log}"        


## both outputs for this rule are identical. done this way for the sake of semantics downstream ##
# rule UCSC_GTF:
#     input: "out/all-merge/stringtie.merged.gtf"
#     output: sample="out/all-merge/sample-UCSC.gtf", reference="out/all-merge/reference.gtf"
#     benchmark: "out/benchmarks/UCSC.txt"
#     params: n="1", R="'span[hosts=1] rusage[mem=10]'", \
#             o="out/logs/UCSC.out", eo="out/logs/UCSC.err", \
#             J="UCSC"
#     shell: "cat {input} | grep chr > {output.reference}; \
#             cat {input} | grep chr > {output.sample}"

rule UCSC_GTF:
    input: "out/all-merge/stringtie.merged.gtf"
    output: sample="out/all-merge/sample-UCSC.gtf", reference="out/all-merge/reference.gtf"
    benchmark: "out/benchmarks/UCSC.txt"
    params: n="1", R="'span[hosts=1] rusage[mem=10]'", J="UCSC", o="out/logs/UCSC.out", eo="out/logs/UCSC.err"
    shell: "cat {REF_GTF} | grep chr > {output.reference}; \
            cat {input} | grep chr > {output.sample} 2> {log}"

rule UCSC_denovo:
        input: "out/all-merge/stringtie.merged.gtf"
        output: sample="out/all-merge/sample-UCSC.gtf"
        benchmark: "out/benchmarks/UCSC.txt"
        log: "out/logs/UCSC.txt"
        params: n="1", R="'span[hosts=1] rusage[mem=10]'", J="UCSC", o="out/logs/UCSC.out", eo="out/logs/UCSC.err"
        shell: "cat {input} | grep chr > {output.sample} 2> {log}"

snakemake.utils.makedirs('out/all-merge/custom_ref')

rule CreateIndelVcf:
    input: SAMPLE_VCF
    output: "out/custom_ref/indels_only.recode.vcf"
    params: n="1", R="'span[hosts=1] rusage[mem=16]'", \
	    o="out/logs/indel_vcf.out", eo="out/logs/indel_vcf.err", \
	    J="indel_vcf", output_prefix="out/custom_ref/indels_only"
    shell: "source activate g2gtools; vcftools --vcf {input} --out {params.output_prefix} --keep-only-indels --recode --recode-INFO-all"
   
rule CreateSnpVcf:
    input: SAMPLE_VCF
    output: "out/custom_ref/snps_only.recode.vcf"
    params: n="1", R="'span[hosts=1] rusage[mem=16]'", \
	    o="out/logs/snp_vcf.out", eo="out/logs/snp_vcf.err", \
	    J="snp_vcf",output_prefix="out/custom_ref/snps_only"
    shell: "source activate g2gtools; vcftools --vcf {input} --out {params.output_prefix} --remove-indels --recode --recode-INFO-all"

rule IndexSeparatedVcfs:
    input: "out/custom_ref/{mutation}_only.recode.vcf"
    output: gz="out/custom_ref/{mutation}_only.recode.vcf.gz",idx="out/custom_ref/{mutation}_only.recode.vcf.gz.tbi"
    benchmark: "out/benchmarks/index_separated_vcfs.txt"
    params: n="1", R="'span[hosts=1] rusage[mem=16]'", \
	    o="out/logs/index_{mutation}_vcf.out", eo="out/logs/index_{mutation}_vcf.err", \
	    J="index_{mutation}_vcf"
    shell: "source activate g2gtools; bgzip {input}; tabix -p vcf {output.gz}"

### NOTE: g2gtools only works on index files with non-ALT contigs,
###       hence the index file maneuvering in the next few rules.

rule LiftIndelsOverRef:
    input: indels_gz="out/custom_ref/indels_only.recode.vcf.gz", ref={REF_FASTA}, idx="out/custom_ref/indels_only.recode.vcf.gz.tbi"
    output: "out/custom_ref/REF-to-" + VCF_INDIVIDUAL + ".chain"
    benchmark: "out/benchmarks/lift_indels.txt"
    params: n="4", R="'span[hosts=1] rusage[mem=8]'", \
	    o="out/logs/lift_indels.out", eo="out/logs/lift_indels.err", \
	    J="lift_indels_over_ref"
    shell: "source activate g2gtools; mv '{input.ref}.fai' '{input.ref}_tmp.fai'; cat '{input.ref}_tmp.fai' | grep 'chr' > '{input.ref}.fai'; g2gtools vcf2chain -f {input.ref} -i {input.indels_gz} -s {VCF_INDIVIDUAL} -o {output}; mv {input.ref}_tmp.fai {input.ref}.fai"

rule PatchRefGenome:
    input: snps_gz="out/custom_ref/snps_only.recode.vcf.gz", ref={REF_FASTA}, idx="out/custom_ref/snps_only.recode.vcf.gz.tbi"
    output: "out/custom_ref/snps_patched.fa"
    params: n="4", R="'span[hosts=1] rusage[mem=8]'", \
	    o="out/logs/patch_ref.out", eo="out/logs/patch_ref.err", \
	    J="patch_ref"
    shell: "source activate g2gtools; g2gtools patch -i {input.ref} -s {VCF_INDIVIDUAL} -v {input.snps_gz} -o {output}; cat '{input.ref}.fai' | grep 'chr' > '{output}.fai'"

rule CreateSampleRef:
    input: patched="out/custom_ref/snps_patched.fa",chain="out/custom_ref/REF-to-"+VCF_INDIVIDUAL+".chain"
    output: "out/custom_ref/sample.fa"
    params: n="4", R="'span[hosts=1] rusage[mem=16]'", \
	    o="out/logs/create_sample_ref.out", eo="out/logs/create_sample_ref.err", \
	    J="create_sample_ref"
    shell: "source activate g2gtools; g2gtools transform -i {input.patched} -c {input.chain} -o {output}" #; cp '{input.patched}.fai' '{output}.fai'"

rule CreateSampleGTF:
    input: custom_ref="out/custom_ref/sample.fa", gtf="out/all-merge/sample-UCSC.gtf", chain="out/custom_ref/REF-to-" + VCF_INDIVIDUAL + ".chain"
    output: "out/all-merge/sample.gtf"
    benchmark: "out/benchmarks/create_sample_gtf.txt"
    params: n="4", R="'span[hosts=1] rusage[mem=8]'", \
	    o="out/logs/create_custom_gtf.out", eo= "out/logs/create_custom_gtf.err", \
	    J="create_custom_gtf"
    shell: "source activate g2gtools; g2gtools convert -c {input.chain} -i {input.gtf} -f gtf -o {output}"

rule CreateSampleDB:
    input: "out/all-merge/sample.gtf"
    output: "out/all-merge/sample.gtf.db"
    benchmark: "out/benchmarks/create_sample_db.txt"
    params: n="1", R="'span[hosts=1] rusage[mem=8]'", \
	    o="out/logs/create_sample_gtf.out", eo= "out/logs/create_sample_gtf.err", \
	    J="create_sample_gtf"
    shell: "source activate g2gtools; g2gtools gtf2db -i {input} -o {output}"

rule TempCopyRefGenome:
    input: REF_DIR+"/"+REF_FASTA_NAME+".fa"
    output: fa=temp("out/custom_ref/reference.fa"),idx="out/custom_ref/reference.fa.fai"
    params: n="1", R="'span[hosts=1] rusage[mem=8]'", \
	    o="out/logs/copy_ref_genome.out", eo= "out/logs/copy_ref_genome.err", \
	    J="copy_ref_genome"
    shell: "cp {input} {WD}/out/custom_ref/reference.fa; cp {input}.fai {WD}/out/custom_ref/reference.fa.fai"

rule gtf_file_to_cDNA_seqs: 
    input: gtf="out/all-merge/{model}.gtf", ref_genome="out/custom_ref/{model}.fa"
    output: fasta="out/all-merge/{model}.transcripts.fasta",
            gtf="out/all-merge/{model}.transcripts.gtf"
    benchmark: "out/benchmarks/gtf_file_to_cDNA_seqs.txt"
    params: n="1", R="'span[hosts=1] rusage[mem=16]'", \
            o="out/logs/{model}.gtf_file_to_cDNA_seqs.out", eo="out/logs/{model}.gtf_file_to_cDNA_seqs.err", \
            J="gtf2cDNA_{model}"
    shell: "gffread {input.gtf} -T -o {output.gtf} \
        --no-pseudo \
        --force-exons \
        -M -Q; \
	gffread -w {output.fasta} -g {input.ref_genome} {output.gtf}"

rule LongOrfs:
    input: "out/all-merge/{model}.transcripts.fasta"#,"out/all-merge/{SAMPLE}.transcripts.fasta.idx"
    output: "out/all-merge/{model}.transcripts.fasta.transdecoder_dir/longest_orfs.pep"
    benchmark: "out/benchmarks/{model}.LongOrfs.json"
    params: n="16", R="'span[ptile=16] rusage[mem=8]'", J="LongOrfs_{model}", \
	    o="out/logs/{model}.LongOrfs.out", eo="out/logs/{model}.LongOrfs.err"
    shell: "cd {WD}/out/all-merge; \
            TransDecoder.LongOrfs \
	    -t {wildcards.model}.transcripts.fasta \
            -p 0"

rule blastp:
    input: "out/all-merge/{model}.transcripts.fasta.transdecoder_dir/longest_orfs.pep"
    output: "out/all-merge/{model}.blastp.outfmt6"
    benchmark: "out/benchmarks/{model}.blastp.json"
    params: n="36", R="'span[hosts=1] rusage[mem=4]'", J="blastp", \
	    o="out/logs/{model}.blastp.out", eo="out/logs/{model}.blastp.err"
    shell: "blastp \
            -num_threads {params.n} \
            -query {input} \
            -db {UNIPROT} \
            -max_target_seqs 1 \
            -outfmt 6 \
            -evalue 1e-5 \
            > {output}"
    
rule Predict:
    input: orfs="out/all-merge/{model}.transcripts.fasta.transdecoder_dir/longest_orfs.pep",
        fasta="out/all-merge/{model}.transcripts.fasta",
        blastp="out/all-merge/{model}.blastp.outfmt6"
    output: "out/all-merge/{model}.transcripts.fasta.transdecoder.pep",
        gff3="out/all-merge/{model}.transcripts.fasta.transdecoder.gff3"
    benchmark: "out/benchmarks/{model}.Predict.json"
    params: n="36", R="'span[hosts=1] rusage[mem=4]'", J="Predict_{model}", \
	    o="out/logs/{model}.Predict.out", eo="out/logs/{model}.Predict.err"
    shell: "cd {WD}/out/all-merge; {TD_DIR}/TransDecoder.Predict \
        -t {wildcards.model}.transcripts.fasta \
        --cpu {params.n} \
	--single_best_orf \
	--retain_blastp_hits {wildcards.model}.blastp.outfmt6"
    
rule gtf_to_alignment_gff3:
    input: "out/all-merge/{model}.transcripts.gtf"
    output: "out/all-merge/{model}.transcripts.gff3"
    benchmark: "out/benchmarks/{model}.gtf_to_alignment_gff3.txt"
    params: n="8", R="'span[ptile=8] rusage[mem=16]'", J="gtf2align_gff3_{model}", \
	    o="out/logs/{model}.gtf_to_alignment_gff3.out", eo="out/logs/{model}.gtf_to_alignment_gff3.err"
    shell: "cufflinks_gtf_to_alignment_gff3.pl {input} > {output}"

rule cdna_alignment_orf_to_genome_orf:
    input: gff3="out/all-merge/{model}.transcripts.gff3",
        fasta_td="out/all-merge/{model}.transcripts.fasta",
        gff3_td="out/all-merge/{model}.transcripts.fasta.transdecoder.gff3"
    output: "out/all-merge/{model}.transcripts.genome.gff3"
    benchmark: "out/benchmarks/{model}.cdna_alignment_orf_to_genome_orf.txt"
    params: n="8", R="'span[ptile=8] rusage[mem=16]'", J="cdna_orf2genome_orf_{model}", \
	    o="out/logs/{model}.cdna_alignment_orf_to_genome_orf.out", eo="out/logs/{model}.cdna_alignment_orf_to_genome_orf.err"
    shell: "cdna_alignment_orf_to_genome_orf.pl {input.gff3_td} {input.gff3} {input.fasta_td} > {output}"

rule gff3_file_to_bed:
    input: "out/all-merge/{model}.transcripts.genome.gff3"
    output: "out/{model}.proteome/proteome.bed"
    #output: "out/all-merge/{model}.proteome.bed"
    benchmark: "out/benchmarks/{model}.gff3_file_to_bed.txt"
    params: n="8", R="'span[ptile=8] rusage[mem=16]'", err="~/error/error.gff3_file_to_bed", J="gff2bed_{model}", \
	    o="out/logs/{model}.gff3_file_to_bed.out", eo="out/logs/{model}.gff3_file_to_bed.err"
    shell: "cat {input} | grep -P \"\tCDS\t\" | gffread --force-exons - -o- | gff3_file_to_bed.pl /dev/stdin | tail -n +2 > {output}"#"| awk '{{ print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12; OFS=\"\t\" }}' > {output}"

rule gff3_file_to_proteins:
    input: gff3="out/all-merge/{model}.transcripts.genome.gff3", ref="out/custom_ref/{model}.fa"
    output: "out/{model}.proteome/proteome.fasta"
    benchmark: "out/benchmarks/{model}.gff3_file_to_proteins.txt"
    params: n="8", R="'span[ptile=8] rusage[mem=16]'", err="~/error/error.gff3_file_to_proteins", J="gff2proteins_{model}", \
	    o="out/logs/{model}.gff3_file_to_proteins.out", eo="out/logs/{model}.gff3_file_to_proteins.err"
    shell: "cat {input.gff3} | grep -P \"\tCDS\t\" | gffread --force-exons - -o- | {TD_DIR}/util/gff3_file_to_proteins_jtp.pl /dev/stdin {input.ref} | egrep -o '^[^*]+' > {output}"

rule reorderFASTA:
    input: "out/{model}.proteome/proteome.fasta"
    output: "out/{model}.proteome/proteome.unique.fasta"
    benchmark: "out/benchmarks/{model}.reorderFASTA.txt"
    params: n="4", R="'span[hosts=1] rusage[mem=4]'", J="reorderFASTA_{model}", \
	    o="out/logs/{model}.reorderFASTA.out", eo="out/logs/{model}.reorderFASTA.err"
    shell: "{RSCRIPT} {SCRIPTS}/reorderFASTA.R {input} {output}"

# rule getfasta:
#     input: "out/all-merge/{SAMPLE}/proteome.bed"
#     output: "out/all-merge/{SAMPLE}/proteome.fasta"
#     benchmark: "out/benchmarks/{SAMPLE}.getfasta.txt"
#     params: pe="1", mem="h_vmem=10G,virtual_free=8G", err="~/error/error.getfasta", name="getfasta"
#     shell: "{BEDTOOLS} getfasta -name -s -fi {REF_FASTA} -bed {input} -fo /dev/stdout | {EMBOSS}/transeq -sequence /dev/stdin -outseq {output}"


# rule PGx_prepare:
#     input: "out/all-merge/{SAMPLE}/transcripts.genome.gff3"
#     output: fasta="out/all-merge/{SAMPLE}/proteome.fasta",
#         bed="out/all-merge/{SAMPLE}/proteome.bed"
#     benchmark: "out/benchmarks/{SAMPLE}.PGx_prepare.txt"
#     params: pe="1", mem="h_vmem=10G,virtual_free=8G", err="~/error/error.PGx_prepare", name="PGx_prepare"
#     shell: "{RSCRIPT} {SCRIPTS}/reorderFASTA.R {WD}/{input} {WD}/{output.fasta} {WD}/{output.bed}"


rule PGx_index:
    input: fasta="out/{model}.proteome/proteome.fasta",
        bed="out/{model}.proteome/proteome.bed"
    output: "out/{model}.proteome/proteome.pickle"
    benchmark: "out/benchmarks/{model}.PGx_index.txt"
    params: n="4", R="'span[hosts=1] rusage[mem=4]'", J="PGx_index_{model}", \
	    o="out/logs/{model}.PGx_index.out", eo="out/logs/{model}.PGx_index.err"
    shell: "{PY} {PGX}/pgx_index.py \
        out/{wildcards.model}.proteome/"

### NEED SEPARATE PEPTIDES QUERY FILE FOR PGx_query AND PGx_bed ###

#rule PGx_query:
#    input: pep="out/all-merge/{SAMPLE}.peptides.txt",
#        index="out/all-merge/{SAMPLE}.proteome.pickle"
#    output: hits="out/all-merge/{SAMPLE}.hits.txt"
#    benchmark: "out/benchmarks/{SAMPLE}.PGx_index.txt"
#    params: n="4", R="'span[hosts=1] rusage[mem=4]'", J="PGx_query_{SAMPLE}", \
#	    o="out/logs/{SAMPLE}.PGx_query.out", eo="out/logs/{SAMPLE}.PGx_query.err"
#    shell: "{PY} {PGX}/pgx_query.py \
#        out/all-merge/ \
#        {input.index} > {output}"


#rule PGx_bed:
#    input: hits="out/all-merge/{SAMPLE}.hits.txt",
#        index="out/all-merge/{SAMPLE}.proteome.pickle"
#    output: bed="out/all-merge/{SAMPLE}.hits.bed"
#    benchmark: "out/benchmarks/{SAMPLE}.PGx_index.txt"
#    params: n="4", R="'span[hosts=1] rusage[mem=4]'", J="PGx_bed_{SAMPLE}", \
#	    o="out/logs/{SAMPLE}.PGx_bed.out", eo="out/logs/{SAMPLE}.PGx_bed.err"
#    shell: "{PY} {PGX}/pgx_bed.py \
#        out/all-merge \
#        {input.hits} \
#        > {output}"

#rule pepdigest:
#    input: "out/all-merge/{db}.db.fasta"
#    output: "out/all-merge/{db}.db.fasta.pepdigest"
#    benchmark: "out/benchmarks/{db}.pepdigest.txt"
#    params: pe="1", mem="h_vmem=10G,virtual_free=8G", err="~/error/error.pepdigest", name="pepdigest"
#    shell: "{EMBOSS}/pepdigest \
#            -seqall {input} \
#            -menu 1 \
#            -outfile {output} \
#            -auto"

DB_PATTERN = "\(= \).*\(.fasta\)"
EXT_PATTERN = "\(= \).*\(pepXML\)"

snakemake.utils.makedirs('out/logs/fragger')

rule MSFragger:
    #input: spectra=expand('{spectra_sample}',spectra_sample=SPECTRA_SAMPLES),\
    input: spectra=MS_RELATIVE_DIR+'/{ms_sample}{run}.mzXML', \
        proteome_db='out/{SAMPLE}.proteome/proteome.unique.fasta', \
	   fragger_params='/home/kwokn/MSFragger_20170103/fragger.params'
    output: MS_RELATIVE_DIR+'/{ms_sample}{run}.{SAMPLE}.db.pepXML'
    benchmark: 'out/benchmarks/{SAMPLE}.MSFragger.txt'
    params: n="16", R="'span[hosts=1] rusage[mem=4]'", J="msfrag_{ms_sample}{run}_{SAMPLE}", \
	o="out/logs/fragger/{ms_sample}{run}.out", eo="out/logs/{ms_sample}{run}.err"
    shell: "cd {MSFRAG}; \
            cat {input.fragger_params} | sed '1s@{DB_PATTERN}@\\1{WD}/{input.proteome_db}@' \
	       | sed '34s@{EXT_PATTERN}@\\1{wildcards.SAMPLE}\.db\.\\2@' > temp.params.{wildcards.SAMPLE}; \
	    java -Xmx32G -jar MSFragger.jar temp.params.{wildcards.SAMPLE} {WD}/{input.spectra}"
	    #rm temp.params.{wildcards.SAMPLE}"

#snakemake --snakefile pgm-lilac --cluster "bsub -J {params.J} -n {params.n} -R {params.R} -W 3:00 -o {params.o} -eo {params.eo}" -j 100 -k --ri --latency-wait 30 -np
#