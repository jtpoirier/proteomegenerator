rule CreateRefSequenceDict:
    input: config['ref']
    output: config['dict']
    params: n="1", R="'span[hosts=1] rusage[mem=18]'", o="out/logs/create_refDict.out", eo="out/logs/create_refDict.err", J="create_refDict"
    shell: "{PICARD} -Xmx16g CreateSequenceDictionary R={input} O={output}"

rule CreateRefSequenceIndex:
    input: config['ref']
    output: config['ref']+".fai"
    params: n="1", R="'span[hosts=1] rusage[mem=18]'", o="out/logs/create_refIdx.out", eo="out/logs/create_refDict.err", J="create_refIdx"
    shell: "{SAMTOOLS} faidx {input}"

rule SambambaMergeAllRGs:
    input: lambda wildcards: expand("out/readgroups/{sample}/{readgroup}.aligned_ubam-merged.bam", readgroup=config["samples"][wildcards.sample]['read_groups'].keys(), sample={wildcards.sample})
    output: temp("out/{sample}.aligned_ubam-merged_RG-merged.bam")
    params: n="36", R="'span[hosts=1] rusage[mem=11]'", o="out/logs/merge_RGs.out", eo="out/logs/merge_RGs.err", J="merge_RGs"
    run:
        if len(input) > 1:
            shell("{SAMBAMBA} merge -t {params.n} {output} {input}")
        else:
            shell("mv {input} {output}")

rule SambambaMarkdups:
    input: temp("out/{sample}.aligned_ubam-merged_RG-merged.bam")
    output: "out/{sample}.aligned_ubam-merged_RG-merged_dedup.bam"
    params: n="36", R="'span[hosts=1] rusage[mem=11]'", o="out/logs/markdups.out", eo="out/logs/markdups.err", J="markdups"
    shell: "{SAMBAMBA} markdup -t {params.n} --tmpdir {TMP} --hash-table-size 1200000 \
            --overflow-list-size 1200000 --io-buffer-size 1024 {input} {output}"

#rule MergeRealignedBamRGsAndMarkDuplicates:
#    input: lambda wildcards: expand("out/readgroups/{sample}.realigned_ubam-merged.bam", sample=config["samples"][wildcards.sample])
#    output: "out/{sample}.realigned_ubam-merged_RGmerged_dedup.bam"
#    params: n="36", R="'span[hosts=1] rusage[mem=11]'", o="out/logs/markdups.out", eo="out/logs/markdups.err", J="markdups"
    #shell: "source activate wgs; \
    #        {PICARD} MarkDuplicates -Xmx390g TMP_DIR={TMP} \
    #          $(echo '{input}' | sed -r \'s/[^ ]+/INPUT=&/g') \
    #	      OUTPUT={output} \
    #          METRICS_FILE='{WD}/out/logs/markdups.metrics' \
    #          VALIDATION_STRINGENCY=SILENT \
    #          OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 \
    #          ASSUME_SORT_ORDER='queryname' \
    #      MAX_RECORDS_IN_RAM=95000000 \
    #          CREATE_MD5_FILE=true" 
    #run: 
   # 	if len(list(config['samples'][wildcards.sample])) > 1:
   # 	    shell("{SAMBAMBA} merge -t {params.n} -l 0 /dev/stdout {' '.join({input})} | \
   # 	    {SAMBAMBA} markdup -t {params.n} -l 5 --tmpdir {TMP} --hash-table-size 1200000 \
   # 	    --overflow-list-size 1200000 --io-buffer-size 80000 /dev/stdin {output}")
   # 	else:
   # 	    shell("{SAMBAMBA} markdup -t {params.n} -l 5 --tmpdir={TMP} --hash-table-size=1200000 \
   # 	    --overflow-list-size=1200000 --io-buffer-size=1024 {input} {output}")

rule SortAndFixRealignedBamTags:
    input: bam="out/{sample}.aligned_ubam-merged_RG-merged_dedup.bam", ref_idx=config['ref']+".fai", ref_dict=config['dict']
    output: bam="out/{sample}.aligned_ubam-merged_RG-merged_dedup_fixedtags.bam", \
            idx="out/{sample}.aligned_ubam-merged_RG-merged_dedup_fixedtags.bai"
    benchmark: "out/benchmarks/SortAndFixTags.txt"
    params: n="36", R="'span[hosts=1] rusage[mem=10]'", \
            o="out/logs/sort_fix_tags.out", eo="out/logs/sort_fix_tags.err", \
            J="sort_fix_tags"
    shell: "{SAMBAMBA} sort -t {params.n} -m 350G --tmpdir {TMP} -o /dev/stdout {input.bam} | \
            {PICARD} SetNmMdAndUqTags -Xmx128g TMP_DIR={TMP} \
              INPUT=/dev/stdin \
              OUTPUT={output.bam} \
              CREATE_INDEX=true \
              CREATE_MD5_FILE=true \
	      MAX_RECORDS_IN_RAM=32000000 \
              REFERENCE_SEQUENCE={config[ref]}"

rule CreateSequenceGroupingTSV:
    input: config['dict']
    output: path="out/sequence_groups.tsv"
    params: n="1", R="'span[hosts=1] rusage[mem=10]'", \
            o="out/logs/create_sequence_groups.out", eo="out/logs/create_sequence_groups.err", \
            J="create_seq_groups"
    run:
            with open(config['dict'],'r') as ref_dict_file:
                    sequence_tuple_list = []
                    longest_sequence = 0
                    for line in ref_dict_file:
                            if line.startswith("@SQ"):
                                    line_split = line.split("\t")
                                    # (Sequence_name, Sequence_length)
                                    sequence_tuple_list.append((line_split[1].split("SN:")[1], int(line_split[2].split("LN:")[1])))
                    longest_sequence = sorted(sequence_tuple_list, key=lambda x: x[1], reverse=True)[0][1]
            # initialize the tsv string with the first sequence
            tsv_string = sequence_tuple_list[0][0]
            temp_size = sequence_tuple_list[0][1]
            for sequence_tuple in sequence_tuple_list[1:]:
                    if temp_size + sequence_tuple[1] <= longest_sequence:
                            temp_size += sequence_tuple[1]
                            tsv_string += "\t" + str(sequence_tuple[0])
                    else:
                            tsv_string += "\n" + str(sequence_tuple[0])
                            temp_size = sequence_tuple[1]
            f = open(output.path,'w')
            f.write(tsv_string)
            f.close()

snakemake.utils.makedirs('out/recal')
snakemake.utils.makedirs('out/logs/recal')

rule BaseRecalibrator:
    input: bam=temp("out/{sample}.aligned_ubam-merged_RG-merged_dedup_fixedtags.bam"),groups_file=config['targets']#,bai="out/{sample_alignment}_readyforBQSR.bai"
    output: "out/recal/{sample}.aligned_ubam-merged_RG-merged_dedup_fixedtags_bqsr-{group}.report"
    params: n="8", R="'span[hosts=1] rusage[mem=4]'", \
	    o="out/logs/recal/{group}_recal.out", eo="out/logs/recal/{group}_recal.err", \
	    J="base_recal_{group}"
    wildcard_constraints: group="\d+"
    shell: "source activate wgs; {JAVA} -Dsamjdk.use_async_io=false \
              -Djava.io.tmpdir={TMP} -Xmx30g -jar {GATK37JAR} -T BaseRecalibrator \
              -nct 8 -R {config[ref]} -I {input.bam} \
	      --useOriginalQualities -o {output} \
              -knownSites {config[snp]} -knownSites {config[indel]} \
	      --interval_padding 100 \
	      -L $(sed -n {wildcards.group}p {input.groups_file} | sed -r 's/\t/ -L /g')"

rule GatherBqsrReports:
    input: expand("out/recal/{{sample}}.aligned_ubam-merged_RG-merged_dedup_fixedtags_bqsr-{group}.report", group=NUM_TARGETS_RANGE)
    output: "out/recal/{sample}.aligned_ubam-merged_RG-merged_dedup_fixedtags_bqsr-full.report"
    params: n="8", R="'span[hosts=1] rusage[mem=10]'", \
	    o="out/logs/gather_reports.out", eo="out/logs/gather_reports.err", \
	    J="gather_reports"
    shell: "source activate wgs; {JAVA} -Xmx32g -cp {GATK37JAR} org.broadinstitute.gatk.tools.GatherBqsrReports \
	      $(echo {input} | sed -r 's/[^ ]+/INPUT=&/g') \
	      OUTPUT={output}"

rule ApplyBQSRToMappedReads:
    input: bqsr="out/recal/{sample}.aligned_ubam-merged_RG-merged_dedup_fixedtags_bqsr-full.report",bam=temp("out/{sample}.aligned_ubam-merged_RG-merged_dedup_fixedtags.bam"), \
	   groups_file=config['targets']
## GOT THIS FAR, NEED TO KEEP GOING ##
    output: temp("out/recal/{sample}.aligned_ubam-merged_RG-merged_dedup_fixedtags_{group}_scatteredrecal.bam")
    wildcard_constraints: group="\d+"
    params: n="8", R="'span[hosts=1] rusage[mem=4]'", \
	    o="out/logs/mapped_bqsr_{group}.out", eo="out/logs/mapped_bqsr{group}.err", \
	    J="mapped_bqsr"
    shell: "source activate wgs; {JAVA} -XX:+PrintFlagsFinal -XX:+PrintGCTimeStamps -XX:+PrintGCDateStamps \
	      -XX:+PrintGCDetails -Xloggc:gc_log.log -Dsamjdk.use_async_io=false \
	      -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xmx30g \
	      -jar {GATK37JAR} -T PrintReads -nct 8 --generate_md5 --keep_program_records \
	      -R {config[ref]} -I {input.bam} --useOriginalQualities -o {output} \
	      --BQSR {input.bqsr} -SQQ 10 -SQQ 20 -SQQ 30 -SQQ 40 --emit_original_quals \
	      -L $(sed -n {wildcards.group}p {input.groups_file} | sed -r 's/\t/ -L /g')"

rule ApplyBQSRToUnmappedReads:
    input: bqsr="out/recal/{sample}.aligned_ubam-merged_RG-merged_dedup_fixedtags_bqsr-full.report",bam=temp("out/{sample}.aligned_ubam-merged_RG-merged_dedup_fixedtags.bam"), \
	   groups_file=config['targets']
    output: temp("out/recal/{sample}.aligned_ubam-merged_RG-merged_dedup_fixedtags_unmapped_scatteredrecal.bam")
    params: n="8", R="'span[hosts=1] rusage[mem=4]'", \
	    o="out/logs/unmapped_bqsr.out", eo="out/logs/unmapped_bqsr.err", \
	    J="unmapped_bqsr"
    shell: "source activate wgs; {JAVA} -XX:+PrintFlagsFinal -XX:+PrintGCTimeStamps -XX:+PrintGCDateStamps \
	      -XX:+PrintGCDetails -Xloggc:gc_log.log -Dsamjdk.use_async_io=false \
	      -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xmx32g \
	      -jar {GATK37JAR} -T PrintReads -nct 8 --generate_md5 --keep_program_records \
	      -R {config[ref]} -I {input.bam} --useOriginalQualities -o {output} \
	      --BQSR {input.bqsr} -SQQ 10 -SQQ 20 -SQQ 30 -SQQ 40 --emit_original_quals \
	      -L unmapped"
	     
rule GatherRecalibratedBAMs:
    input: mapped=expand("out/recal/{{sample}}.aligned_ubam-merged_RG-merged_dedup_fixedtags_{group}_scatteredrecal.bam",group=NUM_TARGETS_RANGE), \
	   unmapped=temp("out/recal/{sample}.aligned_ubam-merged_RG-merged_dedup_fixedtags_unmapped_scatteredrecal.bam")
    output: "out/{sample}.aligned_ubam-merged_RG-merged_dedup_fixedtags_recal.bam"
    params: n="4", R="'span[hosts=1] rusage[mem=9]'", \
	    o="out/logs/gather_recal.out", eo="out/logs/gather_recal.err", \
	    J="gather_recal"
    shell: "source activate wgs; {PICARD} GatherBamFiles -Xmx32g \
              $(echo '{input.mapped}' | sed -r \'s/[^ ]+/INPUT=&/g') \
	      INPUT={input.unmapped} OUTPUT={output} \
	      CREATE_INDEX=true CREATE_MD5_FILE=true"


snakemake.utils.makedirs('out/intervals')
snakemake.utils.makedirs('out/logs/intervals')

rule ScatterIntervals:
    input: config['wgs_calling_regions']
    output: "out/intervals/{interval}-scattered.intervals"
    params: n="2", R="'span[hosts=1] rusage[mem=9]'", \
	    o="out/logs/intervals/{interval}.out", eo="out/logs/intervals/{interval}.err", \
	    J="generate_intervals_{interval}"
    shell: "source activate wgs; {JAVA} -Xmx16g -jar {GATK4JAR} SplitIntervals \
	      -R {config[ref]} -L {input} -scatter {NUM_VARIANT_INTERVALS} \
	      -O {WD}/out/intervals"

rule CallGermlineVariantsWithHTCaller:
    input: bam="out/{sample}.aligned_ubam-merged_RG-merged_dedup_fixedtags_recal.bam", interval_list="out/intervals/{interval}-scattered.intervals"
    output: temp("out/intervals/{sample}.aligned_ubam-merged_RG-merged_dedup_fixedtags_recal_HTC_{interval}.vcf.gz")
    params: n="1", R="'span[hosts=1] rusage[mem=18]'", \
	    o="out/logs/intervals/vcf_{interval}.out", eo="out/logs/intervals/vcf_{interval}.err", \
	    J="generate_vcf_{interval}"
    shell: "source activate wgs; {JAVA} -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xmx16g \
	      -jar {GATK37JAR} -T HaplotypeCaller -R {config[ref]} \
	      -o {output} -I {input.bam} -L {input.interval_list} \
	      -ERC GVCF --max_alternate_alleles 3 \
	      -variant_index_parameter 128000 -variant_index_type LINEAR \
	      -contamination 0 --read_filter OverclippedRead"

if "paired" in config:
    rule MuTect2:
        input: tumor="out/{sample}.{alignment}_recal.bam",normal={PAIRED_SAMPLE},interval_list="out/intervals/{interval}-scattered.intervals" # FLIP 'tumor' and 'normal' IF SAMPLE IS CONTROL 
        output: vcf=temp("out/intervals/{sample}.{alignment}_recal_m2_{interval}.vcf"),
                idx=temp("out/intervals/{sample}.{alignment}_recal_m2_{interval}.vcf.idx")
        benchmark: "out/benchmarks/{interval}.MuTect2.txt"
        params: pe="1", mem="h_vmem=14G,virtual_free=12G", err="~/error/error.MuTect2", out="~/error/output.MuTect2"
        shell: "{JAVA8} -Djava.io.tmpdir={TMP} -Xmx8g -jar {MUTECT2} \
                -T MuTect2 \
                -R {config[ref]} \
                --dbsnp {config[snp]} \
                -L {config[scatter]}/scatter{wildcards.interval}.bed \
                --input_file:tumor {input.tumor} \
                --input_file:normal {input.normal} \
                --interval_padding 100 \
                -o {output.vcf}"
    rule GatherMutectVCFs:
        input: expand(temp("out/intervals/{{sample}}.{{alignment}}_recal_m2_{interval}.vcf"),interval=[str(x).zfill(4) for x in range(NUM_VARIANT_INTERVALS)])
        output: "out/{sample}.{alignment}_recal_m2.vcf"
        shell: "touch {output}"

rule GatherHTCallerGVCFs:
    input: expand(temp("out/intervals/{{sample}}.aligned_ubam-merged_RG-merged_dedup_fixedtags_recal_HTC_{interval}.vcf.gz"),interval=[str(x).zfill(4) for x in range(NUM_VARIANT_INTERVALS)])
    output: "out/{sample}.aligned_ubam-merged_RG-merged_dedup_fixedtags_recal_HTC.vcf"
    params: n="32", R="'span[hosts=1] rusage[mem=4]'", \
	    o="out/logs/merge_gvcfs.out", eo="out/logs/merge_gvcfs.err", \
	    J="merge_gvcfs"
    shell: "source activate wgs; {JAVA} -Xmx64g -jar {GATK37JAR} -nt 16 -T GenotypeGVCFs -R {config[ref]} \
	      $(echo '{input}' | sed -r 's/[^ ]+/--variant &/g') \
              -o {output}"
    
