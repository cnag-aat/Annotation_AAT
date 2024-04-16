from datetime import datetime
import re
import os

rule get_chunks_fasta:
  input:
    fasta = "assembly.fa"
  output:
    expand("assembly.{i}.fa", i=range(1, 31)),
  params:
    numberchunks = 30,
    scripts_dir = "scripts/",
    dirChunks = "chunks/"
  conda:
    '../envs/ann_base.yaml'
  threads: 2
  shell:
      "mkdir -p {params.dirChunks};"
      "cd {params.dirChunks};" 
      "{params.scripts_dir}/fasta2chunks.pl -f {input.fasta} -n {params.numberchunks};"
      "sleep 5m;"

rule merge_gffs:
  input:
    touch = "gene_predictions.ok"
  output:
    out = "gene_predictions.gff3",
  params:
    array_inputs = expand("gene_predictions.gff3.{i}", i=range(1, 31))
  threads: 1
  shell: 
    "cat {params.array_inputs} > {output.out};"
    "sleep 5m;"

rule predictions4EVM:
  input:
    gff = "gene_predictions.gff3"
  output:
    EVM_out = "gene_predictions_preEVM.gff3",
  params:
    reformat_cmd = "$CONDA_PREFIX/opt/evidencemodeler-1.1.1/EvmUtils/misc/augustus_GFF3_to_EVM_GFF3.pl {input.gff} | sed 's/Augustus/augustus_hints/g'",
    create_weights = "echo \"ABINITIO_PREDICTION\tAugustus\t2\">> weights_1.txt;",
    EVM_dir = "step04_EVM.V01/",
    link_out = "step04_EVM.V01/gene_predictions.gff3"
  conda:
    "../envs/evm.yaml"
  threads: 1
  shell:
    "{params.reformat_cmd} > {output.EVM_out};"
    "ln -s {output.EVM_out} {params.link_out};"
    "cd {params.EVM_dir};"
    "{params.create_weights}"
    "sleep 4m;"
