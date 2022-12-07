
################################################
## The main entry point of Snakemake workflow ##
# Run FASTQC on ATACseq samples
# Author: Aristides Zonon
# AFFILIATION: Clermont Auvergne University
# CONTACT: aristides.zonon@etu.uca.fr
################################################
configfile: "config/config.yaml",
import glob
##########################################################
#### Target rule with final output of the workflow
#########################################################
# Collecting runs in 'ech' variable


##### load rules #####
rule all:
    input:
        expand("results/fastqc_init/{sample}_fastqc.zip", sample = config["samples"]),
		expand("results/fastqc_init/{sample}_fastqc.html", sample = config["samples"]),
        expand("results/trim/{subsample}_1_trim_paired.fastq.gz", subsample = config["subsamples"]),
        expand("results/trim/{subsample}_1_trim_unpaired.fastq.gz", subsample = config["subsamples"]),
        expand("results/trim/{subsample}_2_trim_paired.fastq.gz", subsample = config["subsamples"]),
        expand("results/trim/{subsample}_2_trim_unpaired.fastq.gz", subsample = config["subsamples"]),
        expand("results/fastqc_post/{sample}_trim_paired_fastqc.zip", sample = config["samples"]),
        expand("results/fastqc_post/{sample}_trim_paired_fastqc.html", sample = config["samples"]),
        expand("results/bowtie2/{subsample}_trim_mapped_sorted_q2.bam", subsample = config["subsamples"]),
        expand("results/bowtie2/{subsample}_trim_mapped_sorted_q2.bam.bai", subsample = config["subsamples"]),
        expand("results/bowtie2/{subsample}.log", subsample = config["subsamples"]),
		expand("results/picard/{subsample}_trim_mapped_sorted_nodup.bam", subsample = config["subsamples"]),
        expand("results/picard/{subsample}_trim_mapped_sorted_nodup.bam.bai", subsample = config["subsamples"]),
        expand("results/picard/{subsample}_trim_mapped_sorted_dups.txt", subsample = config["subsamples"])
        expand("results/picard/{subsample}.log", subsample = config["subsamples"]),
		expand("results/deeptools/{subsample}_trim_mapped_sorted_q2_nodup.bedgraph", subsample = config["subsamples"])
		expand("{subsample}_trim_mapped_sorted_q2_nodup", subsample = config["subsamples"])

rule fastqc_init:
    input:
        "data/data/subset/{sample}.fastq.gz"
    output:
        html = "results/fastqc_init/{sample}_fastqc.html",
        zip = "results/fastqc_init/{sample}_fastqc.zip"
    conda:
        "envs/qc.yaml"
    threads: 2
    log:
        "log/fastqc_init_{sample}.log"
    message:
        "Executing Initial Quality Control on : {input}"
    shell:
        """
        mkdir -p results/fastqc_init
        fastqc {input} -o "results/fastqc_init" -t {threads} 2>{log}
        """

rule trim:
    input:
        fwd = "data/data/subset/{subsample}_1.fastq.gz",
        rev = "data/data/subset/{subsample}_2.fastq.gz"
    output:
        fwd_paired = "results/trim/{subsample}_1_trim_paired.fastq.gz",
        fwd_unpaired = "results/trim/{subsample}_1_trim_unpaired.fastq.gz",
        rev_paired = "results/trim/{subsample}_2_trim_paired.fastq.gz",
        rev_unpaired = "results/trim/{subsample}_2_trim_unpaired.fastq.gz"
    conda:
        "envs/trim.yaml"
    threads: 6
    log: 
        "log/trim_{subsample}.log"
    message: 
        "Executing Trimming on : {input}"
    shell:
        """
        mkdir -p results/trim
        java -jar /opt/apps/trimmomatic-0.38/trimmomatic-0.38.jar PE -threads {threads} \
            -trimlog results/trim/trim.log -summary results/trim/stats \
            {input.fwd} {input.rev} \
            {output.fwd_paired} {output.fwd_unpaired} {output.rev_paired} {output.rev_unpaired} \
                ILLUMINACLIP:data/data/NexteraPE-PE.fa:2:30:10:2:keepBothReads \
               LEADING:3 \
                TRAILING:3 \
                SLIDINGWINDOW:4:15 \
                MINLEN:33
        """

rule fastqc_post:
    input:
        "results/trim/{sample}_trim_paired.fastq.gz"
    output:
        html = "results/fastqc_post/{sample}_trim_paired_fastqc.html",
        zip = "results/fastqc_post/{sample}_trim_paired_fastqc.zip"
    conda:
        "envs/qc.yaml"
    threads: 2
    log:
        "log/fastqc_post_{sample}.log"
    message:
        "Executing Post Trimming Quality Control on : {input}"
    shell:
        """
        mkdir -p results/fastqc_post
        fastqc {input} -o "results/fastqc_post" -t {threads} 2>{log}
        """

rule bowtie2 :
    input:
        fwd = "results/trim/{subsample}_1_trim_paired.fastq.gz",
        rev = "results/trim/{subsample}_2_trim_paired.fastq.gz",
        idx_bt2 = "data/data/reference/Mus_musculus_GRCm39/fasta/all"
    output:
        bam = "results/bowtie2/{subsample}_trim_mapped_sorted_q2.bam",
        bai = "results/bowtie2/{subsample}_trim_mapped_sorted_q2.bam.bai",
        log = “results/bowtie2/{subsample}_mapped_sorted.log”
    conda:
        "envs/bowtie2.yaml"
    threads: 8
    log:
        "log/bowtie2_{subsample}.log"
    message:
        "Executing Mapping on : {input}"
    shell:
        """
        mkdir -p results/bowtie2
        bowtie2  --very-sensitive -p 8 -k 10  -x {databank}  \
            -1 {input.fwd}  -2 {input.rev}  \
              |  samtools view -q 2 -bS  -  |  samtools sort - -o {output.bam}
        samtools index -b {output.bam}
        samtools idxstats {output.bam} > {log}
        """

rule picard : 
    input:
        "results/bowtie2/{subsample}_trim_mapped_sorted_q2.bam"
    output:
        bam = "results/picard/{subsample}_trim_mapped_sorted_q2_nodup.bam",
        bai = "results/picard/{subsample}_trim_mapped_sorted_q2_nodup.bam.bai",
        txt = "results/picard/{subsample}_trim_mapped_sorted_q2_dups.txt",
		log = "results/picard/{subsample}.log"
    conda:
        "envs/picard.yaml"
    threads: 6
    log:
        "log/picard_{subsample}.log"
    message:
        "Executing Mapping cleaning on : {input}"
    shell:
        """
        mkdir -p results/picard
        java -jar /opt/apps/picard-2.18.25/picard.jar MarkDuplicates I="{input}" O={output.bam} M={output.txt} REMOVE_DUPLICATES=true
        samtools index -b {output.bam}
        samtools idxstats {output.bam} > {log}
        """

rule array : 
    input:
        "results/picard/{subsample}_trim_mapped_sorted_q2_nodup.bam"
    output:
        "results/deeptools/array.npz"
    conda:
        "envs/deeptools.yaml"
    threads: 1
    log:
        "log/array_{subsample}.log"
    message:
        "Executing Data mining on : {input}"
    shell:
        """
        mkdir -p results/deeptools
        multiBamSummary bins --bamfiles {input} -o {output}
        """

rule deeptools : 
    input:
        npz = "results/deeptools/array.npz",
        bam = "results/picard/{subsample}_trim_mapped_sorted_q2_nodup.bam",
        all = "results/picard/*.bam"
    output:
        heatmap = "results/deeptools/heatmap.png",
        bedgraph = "results/deeptools/{subsample}_trim_mapped_sorted_q2_nodup.bedgraph",
        tsv = "results/deeptools/table.tsv",
        hist = "results/deeptools/hist.png"
    conda:
        "envs/deeptools.yaml"
    threads: 6
    log:
        "log/deeptools_{subsample}.log"
    message:
        "Executing Data mining on : {input}"
    shell:
        """
        plotCorrelation -in {input.npz} -p heatmap -c pearson -o {output.heatmap} --removeOutliers
        bamCoverage -b {input.bam} -o {output.bedgraph} -of bedgraph
        bamPEFragmentSize -b {input.all} --samplesLabel 0h_R1 0h_R2 0h_R3 24h_R1 24h_R2 24h_R3 --table {output.tsv} -o {output.hist}
        """

rule macs2 : 
    input:
        "results/picard/{subsample}_trim_mapped_sorted_q2_nodup.bam"
    output:
        shortname = "{subsample}_trim_mapped_sorted_q2_nodup",
        dir = "results/macs2"
    conda:
        "envs/macs2.yaml"
    threads: 6
    log:
        "log/macs2_{subsample}.log"
    message:
        "Executing Identification of DNA access sites on : {input}"
    shell:
        """
        mkdir -p results/macs2
        macs2 callpeak -t {input} -f BAM  -n {output.shortmane} --outdir {output.dir}
        """

rule bedtools : 
    input:
        zero = "results/macs2/ss_50k_0h_*.narrowPeak",
        vingtquatre = "results/macs2/ss_50k_24h_*.narrowPeak"
    output:
        Peakzero = "results/bedptools/50k_0h_trim_mapped_sorted_q2_nodup_combined_peaks.narrowPeak",
        Peakvingtquatre = "results/bedtools/50k_24h_trim_mapped_sorted_q2_nodup_combined_peaks.narrowPeak",
        bedzero = "results/bedtools/50k_0h_trim_mapped_sorted_q2_nodup_combined_merged_peaks.bed",
        bedvingtquatre = "results/bedtools/50k_24h_trim_mapped_sorted_q2_nodup_combined_merged_peaks.bed",
        common = "results/bedtools/50k_trim_mapped_sorted_q2_nodup_combined_merged_common.bed",
        uniquezero = "results/bedtools/50k_0h_trim_mapped_sorted_q2_nodup_combined_merged_unique.bed",
        uniquevingtquatre = "results/bedtools/50k_24h_trim_mapped_sorted_q2_nodup_combined_merged_unique.bed"
    conda:
        "envs/bedtools.yaml"
    threads: 6
    log:
        "log/bedtools_{subsample}.log"
    message:
        "Executing Identification of DNA access sites on : {input}"
    shell:
        """
        mkdir -p results/deeptools
        cat {input.zero} > {output.Peakzero}
        cat {input.vingtquatre} > {output.Peakvingtquatre}
        sort -k1,1 -k2,2n {output.Peakzero} | bedtools merge -i - > {output.bedzero}
        sort -k1,1 -k2,2n {output.Peakvingtquatre} | bedtools merge -i - > {output.bedvingtquatre}
        bedtools intersect -wa -a {input.zero} -b {input.vingtquatre} > {output.common}
        bedtools intersect -wa -v -a {input.zero} -b {input.vingtquatre} > {output.uniquezero}
        bedtools intersect -wa -v -a {input.vingtquatre} -b {input.zero} > {output.uniquevingtquatre}
        """

