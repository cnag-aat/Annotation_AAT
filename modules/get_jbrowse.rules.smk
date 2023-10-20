from datetime import datetime
import re
import os

date = datetime.now().strftime('%Y%m%d.%H%M%S')

rule browse_seq:
    input:
      genome= "genome.fa"
    output:
     # outdir = "seq",
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
        "touch {output.checkpoint};"

rule get_tar:
  input:
    checkpoints = "browse.Geneid.done"
  output:
    tar = "tracks.tar.gz",
  params:
    jbdir = "jbrowse",
    dir = "tracks"
  run:
    if isinstance(params.dir, list):
      nfiles = []
      for file in params.dir:
        nfile = file.replace(params.jbdir, "")
        nfiles.append(nfile)
        #print (nfile)
    else:
      nfiles = params.dir.replace(params.jbdir, "")

    shell("tar -C {params.jbdir} -zcvf {output.tar} {nfiles};"
          "rm -r {params.dir};"
    )

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



