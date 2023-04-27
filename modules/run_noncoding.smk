from datetime import datetime
import os
import re
import subprocess

date = datetime.now().strftime('%Y%m%d.%H%M%S')

##Include modules:

module noncoding_workflow:
  snakefile: "../modules/non_coding.rules.smk"

module preprocess_workflow:
  snakefile: "../modules/preprocessing.rules.smk"

##Prepare run

if config["Parameters"]["run_non_coding"]:
  ncRNA_DIR = config["Outputs"]["ncRNA_annotation_dir"]
  genome_chunks = int(config["Chunks"]["genome_chunks"])
  protein_chunks = int(config["Chunks"]["protein_chunks"])
  dirGenomeChunks = config["Outputs"]["dir_genome_chunks"] 
  if not os.path.exists(dirGenomeChunks):
    os.makedirs(dirGenomeChunks)
  blast_dirs = []
  blast_dirs.append("evidence_blast")
  blast_dirs.append("annotated_blast")

  use rule cmsearch from noncoding_workflow with:
    input:
      genome = masked_reference,
      rfam = config["ncRNA_annotation"]["Rfam"]
    output:
      out = config["Outputs"]["out_cmsearch"]
    conda:
      "../envs/tRNAscan-SEv2.0.11.yaml"      
    log: 
      logs_dir + str(date) + ".j%j.cmsearch.out",
      logs_dir + str(date) + ".j%j.cmsearch.err"
    benchmark:
      logs_dir + str(date) + ".cmsearch.benchmark.txt"
    threads: config["ncRNA_annotation"]["cmsearch_CPUs"]

  # use rule get_chunks_fasta from preprocess_workflow with:
  #   input:
  #     fasta = genome 
  #   output:
  #     expand(dirGenomeChunks + base + ".{i}.fa", i=range(1, genome_chunks+1))
  #   params:
  #     numberchunks = masked_chunks,
  #     scripts_dir = scripts_dir,
  #     dirChunks = dirGenomeChunks
  #   conda:
  #     '../envs/ann_base.yaml'
  #   log:
  #     logs_dir + str(date) + ".j%j.getChunks." + base + ".out",
  #     logs_dir + str(date) + ".j%j.getChunks." + base + ".err",
  #   benchmark:
  #     logs_dir + str(date) + ".getChunks." + base + ".benchmark.txt"
  #   threads: 2  

  use rule tRNAscan from noncoding_workflow with:
    input:
      genome = masked_reference
    output:
      out = config["Outputs"]["out_tRNAscan"] , 
      stats = config["Outputs"]["out_tRNAscan"] + ".stats"
    conda:
      "../envs/tRNAscan-SEv2.0.11.yaml"       
    log:
      logs_dir + str(date) + ".j%j.tRNAscan.out",
      logs_dir + str(date) + ".j%j.tRNAscan.err",
    benchmark:     
      logs_dir + str(date) + ".tRNAscan.benchmark.txt",
    threads: 1  

  use rule lncRNAannotation from noncoding_workflow with:
    input:
      coding = rules.process_update.output.annotation,
      PASA = rules.pasa.output.EVM_out,
      genome = genome,
     # RM_gff = config["Inputs"]["RM_gff"],
      proteins = config["Inputs"]["proteins"],
      pep_annot = rules.process_update.output.pep_annot
    output:
      out = ncRNA_DIR + "lncRNA_annotation.nonclassified.gff3",
      prot_chunks = expand(ncRNA_DIR + "{dirs}/proteins.{i}.fa", dirs = blast_dirs, i=range(1, protein_chunks+1))
    params:
      ncRNA_DIR = ncRNA_DIR,
      scripts_dir = scripts_dir,
      project = config["update"]["project_name"][0],
      version = config["ncRNA_annotation"]["ncRNA_version"],
      chunks = config["Chunks"]["protein_chunks"],
      blast_proteins = ncRNA_DIR + blast_dirs[0],
      blast_annotated = ncRNA_DIR + blast_dirs[1]
    log:
      logs_dir + str(date) + ".j%j.lncRNA_annotation.out",
      logs_dir + str(date) + ".j%j.lncRNA_annotation.err",
    benchmark:     
      logs_dir + str(date) + ".lncRNA_annotation.benchmark.txt",
    conda:
      "../envs/ann_base.yaml"
    threads: 1  

  use rule Blast_prot from noncoding_workflow with:
    input:
      ncRNA_DIR + "{dirs}/proteins.{i}.fa"
    output:    
      ncRNA_DIR + "{dirs}/blast.{i}.out"
    params:
      db = ncRNA_DIR + "lncRNA_trans_db"
    log:
      logs_dir + str(date) + ".j%j.{dirs}.{i}.out",
      logs_dir + str(date) + ".j%j.{dirs}.{i}.err"
    benchmark:
      logs_dir + str(date) + ".{dirs}.{i}.benchmark.txt"
    conda:
      "../envs/ann_base.yaml"
    threads: config["ncRNA_annotation"]["blast_threads"]

  use rule ncAnnotation from noncoding_workflow with:
    input:
      genome = genome,
      blast = expand(ncRNA_DIR + "{dirs}/blast.{r}.out", dirs = blast_dirs, r=range(1, protein_chunks+1)),
#      blast_next = expand(ncRNA_DIR + blast_dirs[1] + "blast.{r}.out", r=range(1, protein_chunks+1)),
      tRNAscan = config["Outputs"]["out_tRNAscan"],
      cmsearch_out = rules.cmsearch.output.out,
      lncRNA = rules.lncRNAannotation.output.out 
    output:
      out_blast = ncRNA_DIR + "proteins.blastp.out",
 #     out_next = ncRNA_DIR + "annotated.blastp.out",
      snc_out =  ncRNA_DIR + project_name[0] + ".snc." + config["ncRNA_annotation"]["ncRNA_version"] + ".gff3",
      out_gff = ncRNA_DIR + project_name[0] + "nc" + config["ncRNA_annotation"]["ncRNA_version"] + ".gff3"
    params:
      ncRNA_DIR = ncRNA_DIR,
      scripts_dir = scripts_dir,
      project = project_name[0],
      version = config["ncRNA_annotation"]["ncRNA_version"]
    log:
      logs_dir + str(date) + ".j%j.obtain_ncRNA_annotation.out",
      logs_dir + str(date) + ".j%j.obtain_ncRNA_annotation.err"
    benchmark:
      logs_dir + str(date) + ".obtain_ncRNA_annotation.benchmark.txt"
    conda:
      "../envs/ann_base.yaml"