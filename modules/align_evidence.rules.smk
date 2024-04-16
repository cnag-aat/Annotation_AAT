from datetime import datetime
import re
import os

rule pasa:
  input:
    genome = 'assembly.fa',
    config_file =  'alignAssembly.v2.5.2.pasa2.config',
    transcripts = 'transcripts_fake.fasta',
    trans_gtf =  'TACO_assembled.gtf',
  output:
    pasa_out = "pasadb.pasa_assemblies.gff3",
    EVM_out = "step04_EVM.V01/transcripts.gff3"
  params:
    aligners = 'gmap', 
    dirPasa = 'step03_annotation_pipeline.V01/protein_and_transcript_mappings/pasa/',
    additional_pasa_opts = " -C ",
    EVM_dir = "step04_EVM.V01/",
    create_weights = "echo \"TRANSCRIPT\tassembler-pasadb\t10\">> weights_1.txt;",
  conda:
    "../envs/pasa2.5.2.yaml"
  threads: 8
  shell:
    "cd {params.dirPasa};"
    "PATH=$PASAHOME/bin:$PATH;"
    "ln -s {input.transcripts} transcripts.fa;"
    "$PASAHOME/Launch_PASA_pipeline.pl -c {input.config_file} -R -g {input.genome} " +\
    " --ALIGNERS {params.aligners} -t transcripts.fa --CPU {threads} {params.additional_pasa_opts};" 
    "ln -s {output.pasa_out} {output.EVM_out};"
    "cd {params.EVM_dir};"
    "{params.create_weights}"

rule Transdecoder:
  input:
    genome = 'assembly.fa',
    transcripts_gff3 = "pasadb.pasa_assemblies.gff3",
    transcripts_fasta = "pasadb.assemblies.fasta"
  output:
    out = "pasadb.assemblies.fasta.transdecoder.genome.gff3",
    EVM_out = "step04_EVM.V01/transdecoder_predictions.gff3"
  params:
    dirPasa = 'step03_annotation_pipeline.V01/protein_and_transcript_mappings/pasa/',
    EVM_dir = "step04_EVM.V01/",
    create_weights = "echo \"OTHER_PREDICTION\ttransdecoder\t4\">> weights_1.txt;",
  conda:
    "../envs/pasa2.5.2.yaml"
  threads: 1
  shell:
    "cd {params.dirPasa};"
    "PATH=$PASAHOME/bin:$PATH;"
    "$PASAHOME/scripts/pasa_asmbls_to_training_set.dbi --pasa_transcripts_fasta {input.transcripts_fasta}" +\
    " --pasa_transcripts_gff3 {input.transcripts_gff3};"
    "cd {params.EVM_dir};"
    "ln -s {output.out} {output.EVM_out};"
    "{params.create_weights}"

rule miniprot:
  input:
    genome = 'assembly.fa',
    proteins = 'proteins.fa'
  output:
    out_cds = "step03_annotation_pipeline.V01/protein_and_transcript_mappings/miniprot/proteins_miniprot_cds.gff3",
  #  out_gene = "step03_annotation_pipeline.V01/protein_and_transcript_mappings/miniprot/proteins_miniprot_gene.gff3",
   # EVM_out = "step04_EVM.V01/proteins.gff3"
  params: 
    miniprot_path = "/software/assembly/src/miniprot/miniprot/",
    miniprot_dir = "step03_annotation_pipeline.V01/protein_and_transcript_mappings/miniprot/",
    additional_opts = " ",  
    create_weights = "echo \"PROTEIN\tminiprot\t8\">> weights_1.txt;",
    EVM_dir = "step04_EVM.V01/",
  threads: 20
  shell: 
    "cd {params.miniprot_dir};"
    "export PATH={params.miniprot_path}:$PATH;"    
    "miniprot -t {threads} -d genome.mpi {input.genome};"
    "miniprot -t {threads} {params.additional_opts} --gff genome.mpi {input.proteins}  > {output.out_gene};"
    "grep -v '#' {output.out_gene} | gawk '$3==\"CDS\"' |" +\
    " perl -ane 'BEGIN{{$c=0;}}chomp; @F =split /\t/, $_; printf \"$F[0]\t$F[1]\t$F[2]\t$F[3]\t$F[4]\t$F[5]\t$F[6]\t$F[7]\tID=%07d;$F[8]\n\", $c; $c++;' > {output.out_cds};"
    "rm genome.mpi;"
    "cd {params.EVM_dir};"
    "ln -s {output.out_cds} {output.EVM_out};"
    "{params.create_weights}"