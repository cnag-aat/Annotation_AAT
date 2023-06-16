from datetime import datetime
import re
import os

date = datetime.now().strftime('%Y%m%d.%H%M%S')

rule browse_seq:
    input:
      genome= "genome.fa"
    output:
      outdir = directory("seq")
    params:
      jbrowse_dir = "jbrowse",
    conda:
        "../envs/JBrowse1.16.yaml"
    threads: 8
    shell:
        "cd {params.jbrowse_dir};"
        "$CONDA_PREFIX/bin/prepare-refseqs.pl --fasta {input.genome}  --out $PWD;"

rule browse_tracks:
    input:
        gff = "geneid_predictions.gff3"
    output:
        outdir = directory("tracks/Geneid")
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

rule get_tar:
  input:
    dir = ["tracks"]
  output:
    tar = "tracks.tar.gz"
  params:
  shell:
    "tar -zcvf {output.tar} {input.dir};"
    "rm -r {input.dir};"

rule get_GCcontent:
  input:
    genome = "genome.fa",
    lengths = "assembly.genome"
  output:
    wig = "genome.GC.wig",
    bw = "genome.GC.bw"
  params:
    jbrowse_dir = "jbrowse/",
    scripts_dir = "/software/assembly/pipelines/Annotation_AAT/scripts/",
    window = "50"
  shell:
    "cd {params.jbrowse_dir};"
    "{params.scripts_dir}gcpct.pl {input.genome} > mean_gc.txt;"
    "{params.scripts_dir}gcCounter -w {params.window} {input.genome} > {output.wig};"
    "{params.scripts_dir}wigToBigWig -clip {output.wig} {input.lengths} {output.bw};"



