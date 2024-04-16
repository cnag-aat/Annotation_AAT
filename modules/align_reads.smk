from datetime import datetime
import re
import os

rule star_index:
  input:
    genome = "genome.fa"
  output: 
    genomedir = directory("star_index"),
    ok = "index.ok"
  params:
    opts = ""
  conda:
    "../envs/star2.7.10a.yaml"
  threads: 8
  shell:
    "mkdir -p {output.genomedir};"
    "STAR --runMode genomeGenerate --genomeDir {output.genomedir} --runThreadN {threads}  --genomeFastaFiles {input.genome} {params.opts};"
    "touch {output.ok};"

rule star:
  input:
    genomedir ="star_index",
    reads = [ "reads.1.fastq.gz", "reads.2.fastq.gz" ]
  output:
    bam = "star/readsAligned.sortedByCoord.out.bam",
  params:
    stardir = "star",
    basename = "reads",
    additional = ""
  conda:
    "../envs/star2.7.10a.yaml"
  threads: 4
  shell:
    "mkdir -p {params.stardir};"
    "cd {params.stardir};"
    "STAR --genomeDir {input.genomedir} --readFilesIn {input.reads} --readFilesCommand zcat --runThreadN {threads} --outFileNamePrefix {params.basename} --outSAMstrandField intronMotif --outSAMtype BAM SortedByCoordinate --outSAMattrIHstart 0 --outFilterIntronMotifs RemoveNoncanonical {params.additional};"
    "echo 'STAR run completed.';"

rule minimap2:
  input:
    genome = "genome.fa",
    reads = "cdna_reads.fastq.gz"
  output:
    bam =  "cDNA/cdna_reads.sorted.bam",
  params:
    basename = "cdna_reads",
    minimap_opts  = ""
  conda:
    "../envs/minimap2.24.yaml"
  threads: 8
  shell:
    "minimap2 -x splice{params.minimap_opts} -t {threads} -a {input.genome} {input.reads} |  samtools sort -@ {threads}  -O BAM -o {output.bam};"
    "samtools index -c {output.bam};"
    "echo 'MINIMAP2 run completed.';"

rule bam2sam:
  input:
    bam = "cDNA/cdna_reads.sorted.bam"
  output:
    sam = "cDNA/cdna_reads.sorted.sam"
  conda:
    "../envs/ESPRESSO1.3.0.yaml"
  threads: 12
  shell:
    "samtools view -@{threads} {input.bam} > {output.sam}; "