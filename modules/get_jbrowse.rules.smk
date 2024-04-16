from datetime import datetime
import re
import os

date = datetime.now().strftime('%Y%m%d.%H%M%S')

rule browse_seq:
  input:
    genome= "genome.fa"
  output:
    checkpoint = "seq.ok"
  params:
    jbrowse_dir = "jbrowse",
  conda:
    "../envs/JBrowse1.16.yaml"
  threads: 8
  shell:
    "cd {params.jbrowse_dir};"
    "$CONDA_PREFIX/bin/prepare-refseqs.pl --fasta {input.genome}  --out $PWD;"
    "touch {output.checkpoint};"
    "sleep 500;"

rule browse_tracks:
  input:
    gff = "geneid_predictions.gff3"
  output:
    checkpoint = "browse.Geneid.done"
  params:
    jbrowse_dir = "jbrowse",
    processed_gff = "geneid_predicions",
    type = "mRNA",
    label = "Geneid",
    className = "pred", 
    setup = ""
  conda:
    "../envs/JBrowse1.16.yaml"
  threads: 8
  shell:
    "cd {params.jbrowse_dir};"
    "{params.setup}"
    "$CONDA_PREFIX/bin/flatfile-to-json.pl --gff {params.processed_gff} --out $PWD " +\
    " --getSubfeatures --getPhase --getLabel --type {params.type} " +\
    " --trackLabel {params.label} --className {params.className};"
    "sleep 500;"
    "touch {output.checkpoint};"

rule get_tar:
  input:
    checkpoints = "browse.Geneid.done"
  output:
    tar = "tracks.tar.gz",
  params:
    jbdir = "jbrowse",
    dir = "tracks"
  shell:
    "tar -C {params.jbdir} -zcvf {output.tar} {params.dir};"
    "rm -r {params.jbdir}{params.dir};"

rule get_GCcontent:
  input:
    genome = "genome.fa",
    lengths = "assembly.genome"
  output:
    wig = "genome.GC.wig",
    bw = "genome.GC.bw"
  params:
    jbrowse_dir = "jbrowse/",
    scripts_dir = "../scripts/",
    window = "50"
  shell:
    "cd {params.jbrowse_dir};"
    "{params.scripts_dir}gcpct.pl {input.genome} > mean_gc.txt;"
    "{params.scripts_dir}gcCounter -w {params.window} {input.genome} > {output.wig};"
    "{params.scripts_dir}wigToBigWig -clip {output.wig} {input.lengths} {output.bw};"

rule get_bw:
  input:
    bam = "RNAseq.bam",
    glen = "assembly.genome"
  output:
    bw = "RNAseq.bw"
  params:
    opt = " -split ",
    dir = "RNA/",
    bg = "RNAseq.bg",
  conda:
    "../envs/bedtools2.30.0.yaml"
  shell:
    "mkdir -p {params.dir};"
    "genomeCoverageBed {params.opt} -ibam {input.bam} -bg | sortBed -i -  > {params.bg};"
    "bedGraphToBigWig {params.bg} {input.glen} {output.bw};"
    "rm {params.bg};"

