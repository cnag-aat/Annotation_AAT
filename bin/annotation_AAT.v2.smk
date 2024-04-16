from datetime import datetime
import os
import re
import subprocess

##Author: Jessica Gomez-Garrido
##CNAG
##email: jessica.gomez@cnag.eu

report: "../report/workflow.rst"

date = datetime.now().strftime('%Y%m%d.%H%M%S')

working_dir = config["Parameters"]["base_dir"]

shell.prefix("TMPDIR=" + working_dir + "tmp/;")
scripts_dir = config["Inputs"]["scripts_dir"]
logs_dir = working_dir + "logs/"
if not os.path.exists(logs_dir):
  os.makedirs(logs_dir)

if not os.path.exists(working_dir + "tmp/"):
  os.makedirs(working_dir + "tmp/")

keepfiles = config["Parameters"]["keep_intermediate"]
base = config["Parameters"]["base_name"]

##0. Define path for files and variables

genome = config["Inputs"]["genome"]
masked_reference = config["Outputs"]["genome_masked"]

EVM_dir = config["Outputs"]["EVM_dir"]
if not os.path.exists(EVM_dir):
  os.makedirs(EVM_dir)

jbrowse_dir = working_dir + "jbrowse/"
jlogs_dir = jbrowse_dir + "logs/"

targets = []
rnabrowser_inputs = []
if config["Repeat Annotation"]["Repeat_library"] or config["Repeat Annotation"]["species_repeat_database"] or config["Parameters"]["run_redmask"]:
  targets.append(masked_reference)

if config["Inputs"]["illumina_dir"] or config["Inputs"]["cDNA_dir"] or config["Inputs"]["dRNA_dir"] or config["Inputs"]["PB_dir"]:
  targets.append(config["Outputs"]["GTF models"])
  targets.append(config["Outputs"]["junctions"])

if config["Parameters"]["run_pasa"]:
  targets.append(EVM_dir + "transcripts.gff3")

if config["Parameters"]["run_transdecoder"]:
  targets.append(EVM_dir + "transdecoder_predictions.gff3")

if config["Parameters"]["run_miniprot"]:
  targets.append(config["Outputs"]["miniprot_cds"])

if config["Parameters"]["run_augustus"]:
  targets.append(config["Outputs"]["augustus_prediction"])

if config["Parameters"]["run_augustus_hints"]:
  targets.append(config["Outputs"]["augustus_hints_prediction"])

if config["Parameters"]["run_geneid"]:
  targets.append(config["Outputs"]["geneid_prediction"])

if config["Parameters"]["run_geneid_introns"]:
  targets.append(config["Outputs"]["geneid_introns_prediction"])

if config["Parameters"]["run_genemark"]:
  targets.append(config["Outputs"]["genemark_prediction"])

if config["Parameters"]["run_genemark-ET"]:
  targets.append(config["Outputs"]["genemark_ET_prediction"])

if config["Parameters"]["run_EVM"]:
  targets.append(config["Outputs"]["evm_out"])

if config["Parameters"]["run_update"]:
  targets.append(config["Outputs"]["update_dir"] + config["Parameters"]["project_name"][0] + config["Parameters"]["project_name"][1] + ".gff3")

if config["Parameters"]["run_non_coding"]:
  targets.append (config["Outputs"]["ncRNA_annotation_dir"] + config["Parameters"]["project_name"][0] + "nc" + config["ncRNA_annotation"]["ncRNA_version"] + ".gff3")

if config["Parameters"]["get_JBrowse"]:
  if not os.path.exists(jlogs_dir):
    os.makedirs(jlogs_dir)
  targets.append(jbrowse_dir + "tracks.tar.gz")
  targets.append(jbrowse_dir + "seq.tar.gz")
  if config["Inputs"]["genome_lengths"]:
    targets.append(jbrowse_dir + base + ".GC.bw")
    if config["Inputs"]["illumina_dir"] or config["Inputs"]["cDNA_dir"] or config["Inputs"]["dRNA_dir"] or config["Inputs"]["PB_dir"]:
      targets.append(jbrowse_dir + "RNA.tar.gz")

#1- Define rule all
rule all:
  input:
    targets
  log:
    logs_dir + str(date) + ".j%j.rule_all.out",
    logs_dir + str(date) + ".j%j.rule_all.err"

#2- Include modules
include: "../modules/run_repeats.smk"
include: "../modules/run_RNA_4annot.smk"
include: "../modules/run_evidence_based.smk"
include: "../modules/gene_predictors_training.smk"
include: "../modules/run_gene_predictors.v2.smk"
include: "../modules/run_combiners.v2.smk"
include: "../modules/run_noncoding.smk"
include: "../modules/get_jbrowse.smk"