from datetime import datetime
import os
import re
import subprocess

date = datetime.now().strftime('%Y%m%d.%H%M%S')

##Include modeules:

module preprocess_workflow:
  snakefile: "../modules/preprocessing.rules.smk"

module align_evidence_workflow:
  snakefile: "../modules/align_evidence.rules.smk"

##Get variables and files:

if genome != None:
  name = os.path.splitext(os.path.basename(genome))[0]

if config["Parameters"]["run_pasa"] or config["Parameters"]["run_transdecoder"] :
  dirPasa = config["Outputs"]["pasa_dir"]
  if not os.path.exists(dirPasa + "/logs"):
    os.makedirs(dirPasa + "/logs")
  pasadb = config["pasa"]["pasadb"]

pasa_weights = config["pasa"]["pasa_weights"]
create_weights_pasa = ""
f = 1
for i in range(len(pasa_weights)):
  create_weights_pasa = create_weights_pasa + "echo \"TRANSCRIPT\tassembler-" + pasadb + "\t" + str(pasa_weights[i]) + "\">> weights_" + str(f) + ".txt;";
  f = f + 1

additional_pasa_opts =  config["pasa"]["add_option"]
if additional_pasa_opts == None:
  additional_pasa_opts = ""
if config["pasa"]["create_database"] == True:
  additional_pasa_opts += " -C"

use rule pasa from align_evidence_workflow with:
  input:
    genome = genome,
    config_file =  config["Inputs"]["pasa_config"],
    transcripts = config["Inputs"]["transcripts"],
    trans_gtf =  config["Inputs"]["trans_gtf"],
  output:
    pasa_out = dirPasa + pasadb + ".pasa_assemblies.gff3",
    EVM_out = EVM_dir + "transcripts.gff3",
    transcripts_fasta =  dirPasa + pasadb + ".assemblies.fasta"
  params:
    aligners = config["pasa"]["aligners"],
    dirPasa = dirPasa,
    EVM_dir = EVM_dir,
    additional_pasa_opts = additional_pasa_opts,
    create_weights = create_weights_pasa, 
  conda:
    '../envs/pasa2.5.2.yaml'
  log:
    dirPasa + "logs/" + str(date) + ".j%j.PASA." + name + ".out",
    dirPasa + "logs/" + str(date) + ".j%j.PASA." + name + ".err",
  benchmark:
    dirPasa + "logs/" + str(date) + ".PASA." + name + ".benchmark.txt"     
  threads: 1

transdecoder_weights = config["transdecoder"]["transdecoder_weights"]
create_weights_transdecoder = ""
f = 1
for i in range(len(transdecoder_weights)):
  create_weights_transdecoder = create_weights_transdecoder + "echo \"OTHER_PREDICTION\ttransdecoder\t" + str(transdecoder_weights[i]) + "\">> weights_" + str(f) + ".txt;";
  f = f + 1
    
use rule Transdecoder from align_evidence_workflow with:
  input:
    genome = genome,
    transcripts_gff3 =  dirPasa + pasadb + ".pasa_assemblies.gff3", 
    transcripts_fasta =  dirPasa + pasadb + ".assemblies.fasta"
  output:
    out = dirPasa + pasadb +  ".assemblies.fasta.transdecoder.genome.gff3",
    EVM_out = EVM_dir + "transdecoder_predictions.gff3"
  params:
    dirPasa = dirPasa,
    EVM_dir = EVM_dir,
    create_weights = create_weights_transdecoder
  conda:
    "../envs/pasa2.5.2.yaml"
  log:
    dirPasa + "logs/" + str(date) + ".j%j.transdecoder." + name + ".out",
    dirPasa + "logs/" + str(date) + ".j%j.transdecoder." + name + ".err",
  benchmark:
    dirPasa + "logs/" + str(date) + ".transdecoder." + name + ".benchmark.txt" 
  threads: 1

if config["Parameters"]["run_miniprot"]:
  dirMiniprotOutput = os.path.dirname(config["Outputs"]["miniprot_gene"])
  if not os.path.exists(dirMiniprotOutput + "/logs"):
    os.makedirs(dirMiniprotOutput + "/logs")

miniprot_weights = config["miniprot"]["miniprot_weights"]
create_weights_miniprot = ""
f = 1
for i in range(len(miniprot_weights)):
  create_weights_miniprot = create_weights_miniprot + "echo \"PROTEIN\tminiprot\t" + str(miniprot_weights[i]) + "\">> weights_" + str(f) + ".txt;";      
  f = f + 1

additional_miniprot_opts =  config["miniprot"]["additional_miniprot_options"]
if additional_miniprot_opts == None:
  additional_miniprot_opts = ""

use rule miniprot from align_evidence_workflow with:
  input:
    genome = genome,
    proteins = config["Inputs"]["proteins"]
  output:
    out_cds = config["Outputs"]["miniprot_cds"],
    out_gene = config["Outputs"]["miniprot_gene"],
    EVM_out = EVM_dir + "proteins.gff3"
  params: 
    miniprot_path = config["miniprot"]["miniprot_path"],
    miniprot_dir = dirMiniprotOutput,
    additional_opts = additional_miniprot_opts,  
    create_weights = create_weights_miniprot,
    EVM_dir = EVM_dir,
  log:
    dirMiniprotOutput + "/logs/" + str(date) + ".j%j.miniprot." + name + ".out",
    dirMiniprotOutput + "/logs/" + str(date) + ".j%j.miniprot." + name + ".err",
  benchmark:
    dirMiniprotOutput + "/logs/" + str(date) + ".miniprot." + name + ".benchmark.txt" 
  threads: config["miniprot"]["miniprot_cores"]
