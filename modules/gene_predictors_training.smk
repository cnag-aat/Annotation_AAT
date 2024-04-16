from datetime import datetime
import os
import re
import subprocess

date = datetime.now().strftime('%Y%m%d.%H%M%S')

##Include modules:
module training_workflow:
  snakefile: "../modules/train_gene_predictors.rules.smk"

##Get variables and files:
training_dir = config["Outputs"]["annotation_basedir"] + "train_gene_predictors/"
tlogs_dir = training_dir + "logs/"
if not os.path.exists(tlogs_dir):
  os.makedirs(tlogs_dir)

##RUN TRAINING RULES
use rule get_training_candidates from training_workflow with:
  input:
    transdecoder = EVM_dir + "transdecoder_predictions.gff3",
    Repeats = config["Outputs"]["Repeat_gff"],
    genome = genome
  output:
    out =  training_dir + "get_candidates/good_candidates.1trans.gff3"
  params:
    outdir = training_dir + "get_candidates/",
    scripts_dir = scripts_dir,
    project = config["Parameters"]["project_name"][0],
    blastdb = config["BLAST"]["blastdb"],
    eval = config["BLAST"]["evalue"]
  log:
    tlogs_dir + str(date) + ".j%j.get_candidates.out",
    tlogs_dir + str(date) + ".j%j.get_candidates.err",
  benchmark:
    tlogs_dir + str(date) + ".get_candidates.benchmark.txt"
  conda:
    "../envs/blast_bedtools.yaml"
  threads: config["BLAST"]["blastCores"]

if config["augustus"]["species"] and not os.path.exists(os.path.dirname(config["Outputs"]["augustus_prediction"]) + "/" + config["augustus"]["species"] + ".trained"):
  use rule train_augustus from training_workflow with:
    input:
      candidates = training_dir + "get_candidates/good_candidates.1trans.gff3",
      genome = masked_reference
    output:
      trained = os.path.dirname(config["Outputs"]["augustus_prediction"]) + "/" + config["augustus"]["species"] + ".trained"
    params:
      outdir = training_dir + 'train_augustus/',
      species = config["augustus"]["species"],
      config_path = config["augustus"]["config_path"]
    log:
      tlogs_dir + str(date) + ".j%j.train_augustus.out",
      tlogs_dir + str(date) + ".j%j.train_augustus.err",
    benchmark:
      tlogs_dir + str(date) + ".train_augustus.benchmark.txt"
    conda:
      "../envs/augustus3.5.0.yaml"
    threads: config["augustus"]["optimize_threads"]