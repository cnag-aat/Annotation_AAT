from datetime import datetime
import re
import os

rule PASA_update:
  input:
    genome = "assembly.fa",
    update_config = "annotCompare_v2.5.2.pasa2.config",
    current_annot = "step04_EVM_V01/evm.best.gff3",
    transcripts = 'transcripts_fake.fasta',
  output:
    annotation_updt = "step05_Annotation_Update.V01/pasa_test.sqlite.PASA_update1.gff3"
  params:
    i = str(1), 
    update_dir = "step05_Annotation_Update.V01/",
    pasadb = "pasa_test.sqlite"
  conda:
    "../envs/pasa2.5.2.yaml"
  threads: 1
  shell:
    "cd {params.update_dir};"
    "PATH=$PASAHOME/bin:$PATH;"
    "$PASAHOME/Launch_PASA_pipeline.pl -c {input.update_config} -A -g {input.genome} " +\
    " -t {input.transcripts} -L --annots {input.current_annot} --CPU {threads};"
    "mv {params.pasadb}.gene_structures_post_PASA_updates.*.gff3 {output.annotation_updt};"

rule process_update:
  input:
    genome = "assembly.fa",
    pasa_updates= "step05_Annotation_Update.V01/pasa_test.sqlite.PASA_update2.gff3",
    EVM_out = "step04_EVM_V01/evm.best.gff3",
  output:
    annotation = "step05_Annotation_Update.V01/TEST1A.gff3",
    stats = "step05_Annotation_Update.V01/TEST1A.stats.txt",
    pep_annot = "step05_Annotation_Update.V01/TEST1A.pep.fa"
  params:
    project = "TEST1",
    version = "A",
    update_dir = "step05_Annotation_Update.V01/",
    scripts_dir = "../scripts/",
    pasadb = "pasa_test.sqlite"
  conda:
    "../envs/perl5.32.1.yaml"
  threads: 1
  shell:
    "cd {params.update_dir};"
    "{params.scripts_dir}process_pasa_update.V3.pl {input.pasa_updates} {input.EVM_out} " +\
    " > {params.pasadb}.processed_updates.gff3;"
    "{params.scripts_dir}assignIDs_generic_fast.pl -a {params.pasadb}.processed_updates.gff3 " +\
    " -project {params.project} -V {params.version};"
    "cat *.{params.project}.{params.version}.gff3 | {params.scripts_dir}/nmd_filter.v2.pl > {output.annotation};"
    "ln -s {input.genome} genome.fa;"
    "export PATH={params.scripts_dir}:$PATH;"
    "CDS2seq.v2.pl -gff {output.annotation} -seq genome.fa;"
    "transcript2seq.pl -gff {output.annotation} -seq genome.fa;"
    "get_longest_peptide.pl {params.project}{params.version}.pep.fa " +\
    " > {params.project}{params.version}.longestpeptide.fa ;"
    "gff_stats.JGG.pl -a {output.annotation} -f genome.fa > {output.stats};"