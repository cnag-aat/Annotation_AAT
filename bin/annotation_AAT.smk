from datetime import datetime
import os
import re
import subprocess

report: "../report/workflow.rst"

date = datetime.now().strftime('%Y%m%d.%H%M%S')

working_dir = config["Outputs"]["base_dir"]
shell.prefix("TMPDIR=" + working_dir + "tmp/; echo tmpdir is set to $TMPDIR;")
scripts_dir = config["Inputs"]["scripts_dir"]
logs_dir = working_dir + "logs/"
if not os.path.exists(logs_dir):
  os.makedirs(logs_dir)

if not os.path.exists(working_dir + "tmp/"):
  os.makedirs(working_dir + "tmp/")

#keepfiles = config["Parameters"]["keep_intermediate"]
base = config["Parameters"]["base_name"]

##0. Define path for files and variables

genome = config["Inputs"]["genome"]
masked_reference = config["Inputs"]["genome_masked"]

EVM_dir = config["Outputs"]["EVM_dir"]
if not os.path.exists(EVM_dir):
  os.makedirs(EVM_dir)

targets = []

if config["Parameters"]["run_augustus"]:
  targets.append(config["Outputs"]["augustus_prediction"])
  targets.append(config["Outputs"]["augustus_preEVM"])

if config["Parameters"]["run_augustus_hints"]:
  targets.append(config["Outputs"]["augustus_hints_prediction"])
  targets.append(config["Outputs"]["augustus_hints_preEVM"])

if config["Parameters"]["run_geneid"]:
  targets.append(config["Outputs"]["geneid_prediction"])
  targets.append(config["Outputs"]["geneid_preEVM"])

if config["Parameters"]["run_geneid_introns"]:
  targets.append(config["Outputs"]["geneid_introns_prediction"])
  targets.append(config["Outputs"]["geneid_introns_preEVM"])

if config["Parameters"]["run_genemark"]:
  targets.append(config["Outputs"]["genemark_prediction"])
  targets.append(config["Outputs"]["genemark_preEVM"])

if config["Parameters"]["run_genemark-ET"]:
  targets.append(config["Outputs"]["genemark_ET_prediction"])
  targets.append(config["Outputs"]["genemark_ET_preEVM"])

if config["Parameters"]["run_pasa"]:
  targets.append(EVM_dir + "transcripts.gff3")

if config["Parameters"]["run_transdecoder"]:
  targets.append(EVM_dir + "transdecoder_predictions.gff3")

if config["Parameters"]["run_miniprot"]:
  targets.append(config["Outputs"]["miniprot_cds"])

if config["Parameters"]["run_EVM"]:
  targets.append(config["Outputs"]["evm_out"])

if config["Parameters"]["run_update"]:
  targets.append(config["Outputs"]["update_dir"] + config["update"]["project_name"][0] + config["update"]["project_name"][1] + ".gff3")

if config["Parameters"]["run_non_coding"]:
  targets.append(config["Outputs"]["out_cmsearch"])
  targets.append( config["Outputs"]["out_tRNAscan"])
  targets.append (config["Outputs"]["ncRNA_annotation_dir"] + config["update"]["project_name"][0] + "nc" + config["ncRNA_annotation"]["ncRNA_version"] + ".gff3")

#1- Define rule all
rule all:
  input:
    targets
  log:
    logs_dir + str(date) + ".j%j.rule_all.out",
    logs_dir + str(date) + ".j%j.rule_all.err"

#2- Include modules
include: "../modules/run_gene_predictors.smk"
include: "../modules/run_evidence_based.smk"
include: "../modules/run_combiners.smk"
include: "../modules/run_noncoding.smk"