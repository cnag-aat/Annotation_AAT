from datetime import datetime
import os
import re
import subprocess

date = datetime.now().strftime('%Y%m%d.%H%M%S')

##Include modules:

module combiners_workflow:
  snakefile: "../modules/combine_data.rules.smk"

module update_workflow:
  snakefile: "../modules/update_annot.rules.smk"

##PREPARE EVM RUN

if config["Parameters"]["run_EVM"]:
  additional_evm_opts =  config["EVM"]["additional_evm_options"]
  if additional_evm_opts == None:
    additional_evm_opts = ""

  total_weight_files = len(miniprot_weights)
  sortcommand = ""
  for i in range(1,total_weight_files+1):
    sortcommand += "sort weights_" + str(i) + ".txt | uniq > weights_" + str(i) + ".clean.txt;"

  evm_predictions = []
  if (config["Parameters"]["run_augustus"]):
    evm_predictions.append(config["Outputs"]["augustus_preEVM"])
  if (config["Parameters"]["run_augustus_hints"]):
    evm_predictions.append(config["Outputs"]["augustus_hints_preEVM"])
  if (config["Parameters"]["run_geneid"]):
    evm_predictions.append(config["Outputs"]["geneid_preEVM"])
  if (config["Parameters"]["run_geneid_introns"]):
    evm_predictions.append(config["Outputs"]["geneid_introns_preEVM"])
  if (config["Parameters"]["run_genemark"]):
    evm_predictions.append(config["Outputs"]["genemark_preEVM"])
  if (config["Parameters"]["run_genemark-ET"]):
    evm_predictions.append(config["Outputs"]["genemark_ET_preEVM"])
  if (config["Parameters"]["run_transdecoder"]):
    evm_predictions.append(EVM_dir + "transdecoder_predictions.gff3")

  other_evm_inputs = []
  evidence_options = ""
  if (config["Parameters"]["run_miniprot"]):
    other_evm_inputs.append(EVM_dir + "proteins.gff3")
    evidence_options += " --protein_alignments " + EVM_dir + "proteins.gff3"
  if (config["Parameters"]["run_pasa"]):
    other_evm_inputs.append(EVM_dir + "transcripts.gff3")
    evidence_options += " --transcript_alignments " + EVM_dir + "transcripts.gff3"

  use rule prepare_evm from combiners_workflow with:
    input:
      genome = genome,
      predictions = evm_predictions,
      other_evm_inputs = other_evm_inputs
    output:
      checkpoint = EVM_dir + "evm_ready",
      predictions_out = EVM_dir + "predictions.gff3",
      genome_out = EVM_dir + "genome.fa",
      weights = expand(EVM_dir + "weights_{w}.clean.txt", w=range(1, total_weight_files+1)),
    params:
      EVM_dir = EVM_dir,
      sortcommand = sortcommand
    log:
      logs_dir +  str(date) + ".prepare_evm.j%j.out",
      logs_dir +  str(date) + ".prepare_evm.j%j.err"
    benchmark:
      logs_dir +  str(date) + ".prepare_evm.benchmark.txt"
    threads: 1

  use rule EVM from combiners_workflow with:
    input:
      weight_file = EVM_dir + "weights_{w}.clean.txt",
      lgenome = EVM_dir + "genome.fa",
      predictions = EVM_dir + "predictions.gff3",
      lproteins = rules.miniprot.output.EVM_out,
      ltranscripts = expand(EVM_dir + "transcripts.gff3", empty = [] if config["Parameters"]["run_pasa"] else [None])
    output:
      evm_models = EVM_dir + "evm_weights_{w}.gff3",
      out = EVM_dir + "evm_weights_{w}.out",
    params:
      set = "{w}",
      evidence_opts = evidence_options,
      additional_evm_opts = additional_evm_opts,
      scripts_dir = scripts_dir
    conda:
      "../envs/evm.yaml"
    log:
      logs_dir +  str(date) + ".evm.{w}.j%j.out",
      logs_dir +  str(date) + ".evm.{w}.j%j.err"
    benchmark:
      logs_dir +  str(date) + ".evm.{w}.benchmark.txt"
    threads: config["EVM"]["evm_cores"]

  use rule select_EVM from combiners_workflow with:
    input:
      models = lambda wildcards: expand(rules.EVM.output.evm_models, w=range(1, total_weight_files+1)),
      base_trans = expand(config["Inputs"]["trans_gtf"] if config["Parameters"]["run_pasa"] else config["Outputs"]["spaln_gene"])
    output:
      EVM_out = config["Outputs"]["evm_out"]
    params:
      EVM_dir = EVM_dir,
      scripts_dir = scripts_dir,
      total_weight_files = total_weight_files
    log:
      logs_dir + str(date) + ".select_evm.j%j.out",
      logs_dir + str(date) + ".select_evm.j%j.err",
    benchmark:
      logs_dir + str(date) + ".select_evm.benchmark.txt"
    conda:
      "../envs/bedtools2.30.0.yaml"
    threads: 1

if config["Parameters"]["run_update"]:
  dirUpdate = config["Outputs"]["update_dir"]
  if not os.path.exists(dirUpdate):
    os.makedirs(dirUpdate)
  project_name = config["update"]["project_name"]
  updt_rounds = config["update"]["update_rounds"]
  round = {}
  round["1"] = config["Outputs"]["evm_out"]
  i = 2
  while i <= updt_rounds:
    n = i-1
    round[str(i)] = dirUpdate + pasadb + ".PASA_update" + str(n) + ".gff3"
    i += 1
  
  use rule PASA_update from update_workflow with:
    input:
      genome = genome,
      update_config = config["Inputs"]["update_config"],
      current_annot = lambda wildcards: round[wildcards.round],
      transcripts = config["Inputs"]["transcripts"]
    output:
      annotation_updt = dirUpdate + pasadb + ".PASA_update{round}.gff3"
    params:
      i = "{round}", 
      update_dir = dirUpdate,
      pasadb = config["pasa"]["pasadb"]
    conda:
      "../envs/pasa2.5.2.yaml"
    log:
      logs_dir + str(date) + ".pasa_updt.{round}.j%j.out",
      logs_dir + str(date) + ".pasa_updt.{round}.j%j.err",
    benchmark:
      logs_dir + str(date) + ".pasa_updt.{round}.benchmark.txt"
    threads: 1

  use rule process_update from update_workflow with:
    input:
      genome = genome,
      pasa_updates = dirUpdate + pasadb + ".PASA_update" + str(updt_rounds) + ".gff3",
      EVM_out = config["Outputs"]["evm_out"],
    output:
      annotation = dirUpdate + config["update"]["project_name"][0] + config["update"]["project_name"][1] + ".gff3",  
      stats = dirUpdate + config["update"]["project_name"][0] + config["update"]["project_name"][1] + ".stats.txt",
      pep_annot = dirUpdate + config["update"]["project_name"][0] + config["update"]["project_name"][1] + ".pep.fa"
    params:
      project = config["update"]["project_name"][0],
      version = config["update"]["project_name"][1],
      update_dir = dirUpdate,
      scripts_dir = scripts_dir,
      pasadb = pasadb
    conda:
      "../envs/perl5.32.1.yaml"
    log:
      logs_dir + str(date) + ".process_updt.j%j.out",
      logs_dir + str(date) + ".process_updt.j%j.err",
    benchmark:
      logs_dir + str(date) + ".process_updt.benchmark.txt"
    threads: 1