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
protein_chunks = int(config["Chunks"]["protein_chunks"])
blast_dirs = []
blast_dirs.append("evidence_blast")
blast_dirs.append("annotated_blast")

##Execute rules
use rule cmsearch from noncoding_workflow with:
  input:
    genome = config["Outputs"]["genome_masked"],
    rfam = config["ncRNA_annotation"]["Rfam"]
  output:
    out = config["Outputs"]["ncRNA_annotation_dir"] + "cmsearch.tbl"
  conda:
    "../envs/tRNAscan-SEv2.0.11.yaml"      
  log: 
    logs_dir + str(date) + ".j%j.cmsearch.out",
    logs_dir + str(date) + ".j%j.cmsearch.err"
  benchmark:
    logs_dir + str(date) + ".cmsearch.benchmark.txt"
  threads: config["ncRNA_annotation"]["cmsearch_CPUs"]

use rule tRNAscan from noncoding_workflow with:
  input:
    genome = config["Outputs"]["genome_masked"]
  output:
    out = config["Outputs"]["ncRNA_annotation_dir"] + "tRNAscan-SE/tRNAscan.out", 
    stats = config["Outputs"]["ncRNA_annotation_dir"] + "tRNAscan-SE/tRNAscan.out.stats"
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
    proteins = config["Inputs"]["proteins"],
    pep_annot = rules.process_update.output.pep_annot
  output:
    out = config["Outputs"]["ncRNA_annotation_dir"] + "lncRNA_annotation.nonclassified.gff3",
    prot_chunks = expand(config["Outputs"]["ncRNA_annotation_dir"] + "{dirs}/proteins.{i}.fa", dirs = blast_dirs, i=range(1, protein_chunks+1))
  params:
    ncRNA_DIR = config["Outputs"]["ncRNA_annotation_dir"],
    scripts_dir = scripts_dir,
    project = config["Parameters"]["project_name"][0],
    version = config["ncRNA_annotation"]["ncRNA_version"],
    chunks = config["Chunks"]["protein_chunks"],
    blast_proteins = config["Outputs"]["ncRNA_annotation_dir"] + blast_dirs[0],
    blast_annotated = config["Outputs"]["ncRNA_annotation_dir"] + blast_dirs[1]
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
    config["Outputs"]["ncRNA_annotation_dir"] + "{dirs}/proteins.{i}.fa"
  output:    
    config["Outputs"]["ncRNA_annotation_dir"] + "{dirs}/blast.{i}.blastout.txt"
  params:
    db = config["Outputs"]["ncRNA_annotation_dir"] + "lncRNA_trans_db",
    evalue = config["BLAST"]["evalue"],
  log:
    logs_dir + str(date) + ".j%j.{dirs}.{i}.out",
    logs_dir + str(date) + ".j%j.{dirs}.{i}.err"
  benchmark:
    logs_dir + str(date) + ".{dirs}.{i}.benchmark.txt"
  conda:
    "../envs/blast_bedtools.yaml"
  threads: config["BLAST"]["blastCores"]

use rule ncAnnotation from noncoding_workflow with:
  input:
    genome = genome,
    blast = expand(config["Outputs"]["ncRNA_annotation_dir"] + "{dirs}/blast.{r}.blastout.txt", dirs = blast_dirs, r=range(1, protein_chunks+1)),
    tRNAscan = config["Outputs"]["ncRNA_annotation_dir"] + "tRNAscan-SE/tRNAscan.out",
    cmsearch_out = rules.cmsearch.output.out,
    lncRNA = rules.lncRNAannotation.output.out 
  output:
    out_blast = config["Outputs"]["ncRNA_annotation_dir"] + "proteins.blastp.out",
    snc_out = config["Outputs"]["ncRNA_annotation_dir"] + config["Parameters"]["project_name"][0] + ".snc." + config["ncRNA_annotation"]["ncRNA_version"] + ".gff3",
    out_gff = config["Outputs"]["ncRNA_annotation_dir"] + config["Parameters"]["project_name"][0] + "nc" + config["ncRNA_annotation"]["ncRNA_version"] + ".gff3"
  params:
    ncRNA_DIR = config["Outputs"]["ncRNA_annotation_dir"],
    scripts_dir = scripts_dir,
    project = config["Parameters"]["project_name"][0],
    version = config["ncRNA_annotation"]["ncRNA_version"]
  log:
    logs_dir + str(date) + ".j%j.obtain_ncRNA_annotation.out",
    logs_dir + str(date) + ".j%j.obtain_ncRNA_annotation.err"
  benchmark:
    logs_dir + str(date) + ".obtain_ncRNA_annotation.benchmark.txt"
  conda:
    "../envs/ann_base.yaml"