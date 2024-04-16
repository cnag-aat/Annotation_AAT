from datetime import datetime
import re
import os

rule tRNAscan:
  input:
    genome = "genome.masked.fa"
  output:
    out = "tRNAscan.out", 
    stats = "tRNAscan.stats"
  params:
  threads: 4
  conda: 
    "../envs/tRNAscan-SEv2.0.11.yaml"
  shell:
      "tRNAscan-SE --thread {threads} --brief -o {output.out} -m {output.stats} {input.genome};"

rule cmsearch:
  input:
    genome = "genome.masked.fa",
    rfam = "Rfam.cm"
  output:
    out = "cmsearch.tbl"
  params:
  conda: 
    "../envs/tRNAscan-SEv2.0.11.yaml"
  threads: 16
  shell:
    "cmsearch --cpu {threads} --tblout {output.out} {input.rfam} {input.genome};"

rule lncRNAannotation:
  input:
    coding = "step05_Annotation_Update.V01/TEST1A.gff3",
    PASA = "step04_EVM.V01/transcripts.gff3",
    genome = "genome.fa",
   # RM_gff = "Repeats.gff3",
    proteins = "proteins.fasta",
    pep_annot = "step05_Annotation_Update.V01/TEST1A.pep.fa"
  output:
    out = "step06_ncRNA_annotation.V01/lncRNA_annotation.nonclassified.gff3",
    prot_chunks = expand("step06_ncRNA_annotation.V01/{dirs}/proteins.{i}.fa", dirs = ["prot_chunks", "annot_chunks"], i=range(1,21))
  params:
    ncRNA_DIR = "step06_ncRNA_annotation.V01",
    scripts_dir = "../scripts/",
    project = "TEST",
    version = "A",
    chunks = 50,
    blast_proteins = "step06_ncRNA_annotation.V01/prot_chunks/",
    blast_annotated = "step06_ncRNA_annotation.V01/annot_chunks/"
  conda:
    "../envs/ann_base.yaml"
  threads: 1
  shell:
    "cd {params.ncRNA_DIR};"
    "{params.scripts_dir}/annotate_lncRNAs.V02.pl {input.coding} {input.PASA} {params.project} {params.version} > lncRNA_annotation.prev2clust.gff3;"
    "{params.scripts_dir}/assignIDs_generic_fast.pl -a lncRNA_annotation.prev2clust.gff3 -project {params.project}.lnc -v {params.version};"
    "gawk \'$3!=\"CDS\"\' *.{params.project}.lnc.{params.version}.gff3 | perl -ane \'chomp;(s/;product=[^\\n]+//); print \"$_\\n\";\' > {output.out};"
    "mv *.{params.project}.lnc.{params.version}.gff3 {params.project}.lnc.{params.version}.gff3;"
    "gawk \'$3==\"transcript\"\' {output.out} > lncRNA_transcripts.nonclassified.gff3;"
    "{params.scripts_dir}/transcript2seq.pl -gff {output.out} -seq {input.genome};"
    "makeblastdb -in lncRNA_annotation.nonclassified.transcripts.fa -dbtype nucl -out lncRNA_trans_db;"
    "mkdir -p {params.blast_proteins};"
    "cd {params.blast_proteins};"
    "ln -s {input.proteins} proteins.fa;"
    "{params.scripts_dir}/fasta2chunks.pl -f proteins.fa -n {params.chunks};"
    "cd ..;"
    "mkdir -p {params.blast_annotated};"
    "cd {params.blast_annotated};"
    "ln -s {input.pep_annot} proteins.fa;"
    "{params.scripts_dir}/fasta2chunks.pl -f proteins.fa -n {params.chunks};"   

rule Blast_prot:
  input:
    "proteins_blast/proteins.1.fa"
  output:    
    "proteins_blast/blast.1.out"
  params:
    db = "lncRNA_trans_db",
    evalue = "1e-06",
  conda:
    "../envs/blast_bedtools.yaml"
  threads: 4
  shell:
    "tblastn -outfmt \"6 qseqid sseqid pident length mismatch gapopen qlen qstart qend slen sstart send evalue bitscore\" -num_threads {threads} " 
    " -query {input} -db {params.db} -evalue {params.evalue} > {output};"

rule ncAnnotation:
  input:
    genome = "genome.fa",
    blast = "proteins_blast/blast.1.out",
    tRNAscan = "tRNAscan.out",
    cmsearch_out = "cmsearch.tbl",
    lncRNA = "step06_ncRNA_annotation.V01/lncRNA_annotation.nonclassified.gff3"
  output:
    out_blast = "proteins.blastp.out",
    snc_out =  "TEST.snc.A.gff3",
    out_gff = "TEST1ncA.gff3"
  params:
    ncRNA_DIR = "step06_ncRNA_annotation.V01",
    scripts_dir = "../scripts/",
    project = "TEST",
    version = "A"
  conda:
    "../envs/ann_base.yaml"  
  shell:
    "cd {params.ncRNA_DIR};"
    "cat {input.blast} > {output.out_blast} ;"
    "cat {output.out_blast} | cut -f 2 | perl -ane \'if ($_ =~ s/T([0-9]+)(\\n)//) {{print \"$_\\n\";}}\' | sort |uniq > {params.ncRNA_DIR}/pseudogenes.ids;"
    "{params.scripts_dir}/annotate_ncRNAs.V2.pl {input.cmsearch_out} {input.tRNAscan} {params.project} {params.version} >  {output.snc_out};"
    "gawk \'$3==\"ncRNA\"\' {output.snc_out} > {params.ncRNA_DIR}smallncRNA_annotation.transcripts.gff3;"
    "intersectBed -s -f 0.80 -wo -a {params.ncRNA_DIR}lncRNA_transcripts.nonclassified.gff3 -b {params.ncRNA_DIR}smallncRNA_annotation.transcripts.gff3 >  {params.ncRNA_DIR}lncRNA_BT_small_0.80.out;"   
    "cat {params.ncRNA_DIR}lncRNA_BT_small_0.80.out | perl -ane \'if ($F[8] =~ /ID=([^;]+)/) {{$tr=\"$1\";}} $tr =~ s/T([0-9]+)$//; if ($F[17] =~ /ID=([^;]+)/){{$id=$1;}} print \"$tr\\t$id\\n\";\' > {params.ncRNA_DIR}small.ids;"
    "cat {input.lncRNA} |   perl -ne \'BEGIN{{open(ID_PS,\"<pseudogenes.ids\");while(<ID_PS>){{chomp; $ids_pseudo{{$_}}++;}}close ID_PS;open(ID_SMALL, \"<small.ids\"); while (<ID_SMALL>){{chomp; @line_small=split /\t/, $_; $ids_small{{$line_small[0]}}=$line_small[1];}} close ID_SMALL;}}chomp; if(m/ID=([^;]+)/){{$tr=$1; $tr =~ s/T([0-9]+)$//;if (exists $ids_pseudo{{$tr}}){{$_=~s/Description=lncRNA/Description=lncRNA;Class=pseudogene/;}}elsif (exists $ids_small{{$tr}}) {{$_=~s/Description=lncRNA/Description=lncRNA;Class=small_ncRNA/;}}}}print \"$_\\n\";\' > lncRNA_annotation.classified.gff3;"
    "{params.scripts_dir}merge_ncRNA_lncRNA.pl {output.snc_out} lncRNA_annotation.classified.gff3 {params.project} {params.version} > {output.out_gff};"
    "{params.scripts_dir}transcript2seq.pl -gff {output.out_gff} -seq {input.genome};"