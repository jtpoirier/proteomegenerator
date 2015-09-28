# variables
WD = "/ifs/e63data/poirierj/AK"
workdir: WD
PY = '/ifs/opt/common/CentOS_6/python/python-2.7.8/bin/python'
REF = '/ifs/e63data/poirierj/indexes/GRCh38/STAR'
DICT = '/ifs/e63data/poirierj/indexes/GRCh38/BWAIndex/GRCh38.genome.dict'
STAR = '/ifs/e63data/poirierj/bin/STAR/bin/Linux_x86_64_static/STAR'
SF='/ifs/e63data/poirierj/bin/STAR-Fusion'
PICARD = '/ifs/e63data/poirierj/bin/picard-tools-1.134/picard.jar'
STRINGTIE='/ifs/e63data/poirierj/bin/stringtie-1.0.4.Linux_x86_64/stringtie'
GTF='/ifs/e63data/poirierj/indexes/GRCh38/gencode.v20.annotation.gtf'
CUFFLINKS='/ifs/e63data/poirierj/bin/cufflinks-2.2.1.Linux_x86_64'
FASTA='/ifs/e63data/poirierj/indexes/GRCh38/GRCh38.genome.fa'
TD="/ifs/e63data/poirierj/bin//TransDecoder"
BLASTP="/ifs/e63data/poirierj/bin/ncbi-blast-2.2.31+/bin/blastp"
HMMSCAN="/ifs/e63data/poirierj/bin/hmmer-3.1b2-linux-intel-x86_64/binaries/hmmscan"
UNIPROT="/ifs/e63data/poirierj/indexes/uniprot/UP000005640.fasta"
PFAM="/ifs/e63data/poirierj/indexes/pfam/Pfam-A.hmm"
RSCRIPT='/ifs/opt/common/CentOS_6/R/R-3.2.0/bin/Rscript'
SCRIPTS='/home/poirierj/scripts'
TMP = '/scratch/poirierj'

GFP = 'Sample_AH_RNAseq_RPE-GFP_01 Sample_AH_RNAseq_RPE-GFP_02 Sample_AH_RNAseq_RPE-GFP_03'.split()
GFPP5 = 'Sample_AH_RNAseq_RPE-GFP-P5_01 Sample_AH_RNAseq_RPE-GFP-P5_02 Sample_AH_RNAseq_RPE-GFP-P5_03'.split()
PGBD5P5 = 'Sample_AH_RNAseq_GFP-P5-Tumor_01 Sample_AH_RNAseq_GFP-P5-Tumor_02 Sample_AH_RNAseq_GFP-P5-Tumor_03'.split()
K0562 = 'FCH9EFLADXX-HUMbghEAACRAAPEI-225'.split()

SAMPLES = GFP + GFPP5 + PGBD5P5
GROUPS="all GFP GFPP5 PGBD5P5".split()

snakemake.utils.makedirs('/scratch/poirierj')
snakemake.utils.makedirs('out/benchmarks')
snakemake.utils.makedirs('out/all-cuffmerge')

#benchmarking figures
onsuccess:
    print("Workflow finished, no error.")
    shell("{RSCRIPT} {SCRIPTS}/stats.R {WD}/out/benchmarks")
onerror:
    print("Workflow finished with errors.")
    shell("{RSCRIPT} {SCRIPTS}/stats.R {WD}/out/benchmarks")

# rule all:
#     input: "out/all-cuffmerge/transcripts.fasta.transdecoder.pep"#,expand("out/{sample}.fusion_candidates.final",sample=SAMPLES)

rule all:
    input: expand("out/{sample}.fusion_candidates.final",sample=SAMPLES), "out/all-cuffmerge/transcripts.fasta.transdecoder.pep"

rule STAR:
    input: "/ifs/e63data/kentsislab/RKoche/UpdateSept2015_RNAseq_RPE_P5/{sample}/{sample}_cat_1.fastq.gz", "/ifs/e63data/kentsislab/RKoche/UpdateSept2015_RNAseq_RPE_P5/{sample}/{sample}_cat_2.fastq.gz"
    output: "out/{sample}.Aligned.sortedByCoord.out.bam", "out/{sample}.Chimeric.out.junction"
    benchmark: "out/benchmarks/{sample}.align.json"
    params: pe="2", mem="h_vmem=20G,virtual_free=18G", err="~/tmp/error.snake", name="align"
    shell: "{STAR} \
        --genomeDir {REF} \
        --readFilesIn {input} \
        --outFileNamePrefix out/{wildcards.sample}. \
        --outSAMattributes NH HI \
        --outSAMattrRGline ID:{wildcards.sample} LB:1 PL:illumina PU:1 SM:{wildcards.sample} \
        --runThreadN {params.pe} \
        --outSAMtype BAM SortedByCoordinate \
        --clip3pAdapterSeq AGATCGGAAGAG \
        --readFilesCommand zcat \
        --outSAMmapqUnique 60 \
        --twopassMode Basic \
        --chimOutType WithinBAM \
        --outSAMstrandField intronMotif \
        --outFilterIntronMotifs RemoveNoncanonical \
        --outReadsUnmapped None \
        --chimSegmentMin 15 \
        --chimJunctionOverhangMin 15 \
        --alignMatesGapMax 200000 \
        --alignIntronMax 200000"

rule BuildBamIndex:
    input: "out/{sample}.Aligned.sortedByCoord.out.bam"
    output: "out/{sample}.Aligned.sortedByCoord.out.bai"
    benchmark: "out/benchmarks/{sample}.BuildBamIndex.json"
    params: pe="1", mem="h_vmem=10G,virtual_free=8G", err="~/tmp/error.snake", name="BuildBamIndex"
    shell: "java -Djava.io.tmpdir={TMP} -Xmx6g -jar {PICARD} \
            BuildBamIndex \
            INPUT={input}"

rule StringTie:
    input: bam="out/{sample}.Aligned.sortedByCoord.out.bam", bai="out/{sample}.Aligned.sortedByCoord.out.bai"
    output: "out/{sample}-stringtie.gtf"
    benchmark: "out/benchmarks/{sample}.StringTie.json"
    params: pe="1", mem="h_vmem=20G,virtual_free=18G", err="~/tmp/error.snake", name="StringTie"
    shell: "{STRINGTIE} \
            {input.bam} \
            -p {params.pe} \
            -G {GTF} \
            -o {output}"

snakemake.utils.makedirs('out/all-cuffmerge')
            
rule compose_merge:
    input: expand("out/{sample}-stringtie.gtf",sample=SAMPLES)
    output: txt="out/all-cuffmerge/conditions.txt"
    params: pe="1", mem="h_vmem=10G,virtual_free=8G", err="~/tmp/error.snake", name="cuffmerge"
    run:
        with open(output.txt, "w") as out:
            out.write("\n".join(input))
            
rule cuffmerge:
    input: "out/all-cuffmerge/conditions.txt"
    output: "out/all-cuffmerge/merged.gtf"
    benchmark: "out/benchmarks/cuffmerge.json"
    params: pe="12", mem="h_vmem=4G,virtual_free=3G", err="~/tmp/error.cuffmerge", name="cuffmerge"
    shell: "{PY} {CUFFLINKS}/cuffmerge \
    	-o out/all-cuffmerge \
    	--ref-gtf {GTF} \
    	--ref-sequence {FASTA} \
    	--num-threads {params.pe} \
    	{input}"

rule UCSC:
    input: "out/all-cuffmerge/merged.gtf"
    output: "out/all-cuffmerge/merged-UCSC.gtf"
    benchmark: "out/benchmarks/UCSC.json"
    params: pe="1", mem="h_vmem=10G,virtual_free=8G", err="~/tmp/error.snake", name="UCSC"
    shell: "cat {input} | grep chr > {output}"

rule gtf_file_to_cDNA_seqs:
    input: "out/all-cuffmerge/merged-UCSC.gtf"
    output: "out/all-cuffmerge/transcripts.fasta"
    benchmark: "out/benchmarks/gtf_file_to_cDNA_seqs.json"
    params: pe="1", mem="h_vmem=10G,virtual_free=8G", err="~/tmp/error.gtf_file_to_cDNA_seqs", name="gtf_file_to_cDNA_seqs"
    shell: "{SF}/util/gtf_file_to_cDNA_seqs.pl {input} {FASTA} > {output}"

rule index_cdna_seqs:
    input: "out/all-cuffmerge/transcripts.fasta"
    output: "out/all-cuffmerge/transcripts.fasta.idx"
    benchmark: "out/benchmarks/index_cdna_seqs.json"
    params: pe="1", mem="h_vmem=20G,virtual_free=18G", err="~/tmp/error.idx", name="idx"
    shell: "{SF}/util/index_cdna_seqs.pl {input}"

# rule cdna_fasta:
#     input: "out/all-cuffmerge/merged-UCSC.gtf"
#     output: "out/all-cuffmerge/transcripts.fasta"
#     benchmark: "out/benchmarks/cdna_fasta.json"
#     params: pe="1", mem="h_vmem=10G,virtual_free=8G", err="~/tmp/error.snake", name="cdna_fasta"
#     shell: "perl {TD}/util/cufflinks_gtf_genome_to_cdna_fasta.pl \
#         {input} \
#         {FASTA} \
#         > {output}"

rule gff3:
    input: "out/all-cuffmerge/merged-UCSC.gtf"
    output: "out/all-cuffmerge/merged-UCSC.gff3"
    benchmark: "out/benchmarks/gff3.json"
    params: pe="1", mem="h_vmem=10G,virtual_free=8G", err="~/tmp/error.snake", name="gff3"
    shell: "perl {TD}/util/cufflinks_gtf_to_alignment_gff3.pl \
        {input} \
        > {output}"

rule LongOrfs:
    input: "out/all-cuffmerge/transcripts.fasta"
    output: "out/all-cuffmerge/transcripts.fasta.transdecoder_dir/longest_orfs.pep"
    benchmark: "out/benchmarks/LongOrfs.json"
    params: pe="1", mem="h_vmem=10G,virtual_free=8G", err="~/tmp/error.snake", name="LongOrfs"
    shell: "cd out/all-cuffmerge ; {TD}/TransDecoder.LongOrfs -t transcripts.fasta"

rule blastp:
    input: "out/all-cuffmerge/transcripts.fasta.transdecoder_dir/longest_orfs.pep"
    output: "out/all-cuffmerge/blastp.outfmt6"
    benchmark: "out/benchmarks/blastp.json"
    params: pe="12", mem="h_vmem=4G,virtual_free=3G", err="~/tmp/error.snake", name="blastp"
    shell: "{BLASTP} \
        -num_threads {params.pe} \
        -query {input}  \
        -db {UNIPROT}  \
        -max_target_seqs 1 \
        -outfmt 6 \
        -evalue 1e-5 \
        > {output}"

rule pfam:
    input: "out/all-cuffmerge/transcripts.fasta.transdecoder_dir/longest_orfs.pep"
    output: "out/all-cuffmerge/pfam.domtblout"
    benchmark: "out/benchmarks/pfam.json"
    params: pe="12", mem="h_vmem=10G,virtual_free=8G", err="~/tmp/error.snake", name="pfam"
    shell: "{HMMSCAN} \
        --cpu {params.pe} \
        --domtblout {output} \
        {PFAM} \
        {input}"

rule Predict:
    input: orfs="out/all-cuffmerge/transcripts.fasta.transdecoder_dir/longest_orfs.pep", fasta="out/all-cuffmerge/transcripts.fasta", blastp="out/all-cuffmerge/blastp.outfmt6"#, pfam="out/all-cuffmerge/pfam.domtblout"
    output: "out/all-cuffmerge/transcripts.fasta.transdecoder.pep"
    benchmark: "out/benchmarks/Predict.json"
    params: pe="1", mem="h_vmem=225G,virtual_free=200G", err="~/tmp/error.snake", name="Predict"
    shell: "cd out/all-cuffmerge ; {TD}/TransDecoder.Predict \
        -t transcripts.fasta \
        --retain_blastp_hits blastp.outfmt6"# \
 #       --retain_pfam_hits pfam.domtblout"

rule orf:
    input: TD="out/all-cuffmerge/transcripts.fasta.transdecoder.gff3", UCSC="out/all-cuffmerge/merged-UCSC.gff3", fasta="out/all-cuffmerge/transcripts.fasta"
    output: "out/all-cuffmerge/transcripts.fasta.transdecoder.genome.gff3"
    benchmark: "out/benchmarks/{sample}.orf.json"
    params: pe="1", mem="h_vmem=10G,virtual_free=8G", err="~/tmp/error.snake", name="orf"
    shell: "perl {TD}/util/cdna_alignment_orf_to_genome_orf.pl \
        {input.TD} \
        {input.UCSC} \
        {input.fasta} \
        > {output}"

rule fusion:
    input: GTF="out/all-cuffmerge/merged-UCSC.gtf", fasta="out/all-cuffmerge/transcripts.fasta", idx="out/all-cuffmerge/transcripts.fasta.idx", junc="out/{sample}.Chimeric.out.junction"
    output: "out/{sample}.fusion_candidates.final"
    benchmark: "out/benchmarks/{sample}.fusion.json"
    params: pe="1", mem="h_vmem=10G,virtual_free=8G", err="~/tmp/error.fusion", name="fusion"
    shell: "{SF}/STAR-Fusion \
        -G {input.GTF} \
        -C {input.fasta} \
        -J {input.junc} \
        --out_prefix out/{wildcards.sample} \
        --tmpdir {TMP}"
        
# R script to fetch peptide sequences and write to fasta

# rule translate_fusions:
#     input: "out/{sample}.fusion_candidates.final"
#     output: "out/{sample}.fusions.pep"
#     benchmark: "out/benchmarks/{sample}.translateFusion.json"
#     params: pe="1", mem="h_vmem=10G,virtual_free=8G", err="~/tmp/error.translateFusion", name="translateFusion"
#     shell: "{RSCRIPT} {SCRIPTS}/translateFusion.R {WD}/{input}"

# cat fusion fasta to end of *.pep file or create a new file (not that big)

#snakemake --cluster "qsub -V -pe smp {params.pe} -l {params.mem} -e {params.err} -o ~/tmp/fusion.out" --jn {rulename}.{jobid}.sj --stats ~/Projects/AK/stats.json -j 100 -k -n
# --rulegraph | dot -Tpdf > dag.pdf