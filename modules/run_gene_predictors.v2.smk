from datetime import datetime
import os
import re
import subprocess

date = datetime.now().strftime('%Y%m%d.%H%M%S')

##Include modules:
module preprocess_workflow:
  snakefile: "../modules/preprocessing.rules.smk"

module gene_predictors_workflow:
  snakefile: "../modules/gene_predictors.rules.smk"

module combiners_workflow:
  snakefile: "../modules/combine_data.rules.smk"

##Get variables and files:
if masked_reference != None:
  name_masked = os.path.splitext(os.path.basename(masked_reference))[0]

masked_chunks = int(config["Chunks"]["masked_chunks"])
dirMaskedChunks = config["Outputs"]["dir_masked_chunks"] 

# predictions_list = {}
# predictions_links = {}
weights = {}
reformat = {}

geneid_weights = config["geneid"]["weights"]
create_weights_geneid = ""
if config["Parameters"]["run_geneid"]:
  f = 1
  for i in range(len(geneid_weights)):
    create_weights_geneid = create_weights_geneid + "echo \"ABINITIO_PREDICTION\tgeneid_v1.4\t" + str(geneid_weights[i]) + "\">> weights_" + str(f) + ".txt;";      
    f = f + 1
weights[config["Outputs"]["geneid_preEVM"]] = create_weights_geneid


dirAugustusOutput = os.path.dirname(config["Outputs"]["augustus_prediction"]) + "/"
augustus_weights = config["augustus"]["weights"]
create_weights_aug = ""
if config["Parameters"]["run_augustus"]:
  if not os.path.exists(dirAugustusOutput + "/logs/"):
    os.makedirs(dirAugustusOutput + "/logs/")
  f = 1
  for i in range(len(augustus_weights)):
    create_weights_aug = create_weights_aug + "echo \"ABINITIO_PREDICTION\tAugustus\t" + str(augustus_weights[i]) + "\">> weights_" + str(f) + ".txt;";      
    f = f + 1
weights[config["Outputs"]["augustus_preEVM"]] = create_weights_aug
reformat[config["Outputs"]["augustus_preEVM"]] = config["EVM"]["evm_path"] + "EvmUtils/misc/augustus_GFF3_to_EVM_GFF3.pl " + config["Outputs"]["augustus_prediction"]

genemark_weights = config["genemark"]["weights"]
create_weights_gmk = ""
if config["Parameters"]["run_genemark"]:
  f = 1
  for i in range(len(genemark_weights)):
    create_weights_gmk = create_weights_gmk + "echo \"ABINITIO_PREDICTION\tGeneMark.hmm3\t" + str(genemark_weights[i]) + "\">> weights_" + str(f) + ".txt;";      
    f = f + 1
weights[config["Outputs"]["genemark_preEVM"]] = create_weights_gmk

geneid_introns_weights = config["geneid_introns"]["weights"]
create_weights_geneid_introns = ""
if config["Parameters"]["run_geneid_introns"]:
  f = 1
  for i in range(len(geneid_introns_weights)):
    create_weights_geneid_introns = create_weights_geneid_introns + "echo \"ABINITIO_PREDICTION\tgeneid_introns\t" + str(geneid_introns_weights[i]) + "\">> weights_" + str(f) + ".txt;";      
    f = f + 1
weights[config["Outputs"]["geneid_introns_preEVM"]] = create_weights_geneid_introns

genemark_et_weights = config["genemark_ET"]["genemark_ET_weights"]
create_weights_gmk_et = ""
if config["Parameters"]["run_genemark-ET"]:
  f = 1
  for i in range(len(genemark_et_weights)):
    create_weights_gmk_et = create_weights_gmk_et + "echo \"ABINITIO_PREDICTION\tGeneMark-ET\t" + str(genemark_et_weights[i]) + "\">> weights_" + str(f) + ".txt;";      
    f = f + 1
weights[config["Outputs"]["genemark_ET_preEVM"]] = create_weights_gmk_et

dirAugustusHintsOutput = os.path.dirname(config["Outputs"]["augustus_hints_prediction"]) + "/"
augustus_hints_weights = config["augustus_hints"]["augustus_hints_weights"]
create_weights_aug_hints = ""
if config["Parameters"]["run_augustus_hints"]:
  if not os.path.exists(dirAugustusHintsOutput + "/logs/"):
    os.makedirs(dirAugustusHintsOutput + "/logs")
  f = 1
  for i in range(len(augustus_hints_weights)):
    create_weights_aug_hints = create_weights_aug_hints + "echo \"ABINITIO_PREDICTION\taugustus_hints\t" + str(augustus_hints_weights[i]) + "\">> weights_" + str(f) + ".txt;";      
    f = f + 1
weights[config["Outputs"]["augustus_hints_preEVM"]] = create_weights_aug_hints
reformat[config["Outputs"]["augustus_hints_preEVM"]] = config["EVM"]["evm_path"]  + "EvmUtils/misc/augustus_GFF3_to_EVM_GFF3.pl " + config["Outputs"]["augustus_hints_prediction"] + " | sed 's/Augustus/augustus_hints/g'"

#Execute rules
use rule geneid from gene_predictors_workflow with:
  input:
    fasta = masked_reference,
    geneid_parameters = config["geneid"]["parameters"] if config["geneid"]["parameters"] else []
  output:
    out = config["Outputs"]["geneid_prediction"],
    EVM_out = config["Outputs"]["geneid_preEVM"]
  params:
    geneid_options = config["geneid"]["options"],
    path =  config["geneid"]["path"],
    EVM_dir = EVM_dir,
    link_out = EVM_dir + "geneid_predictions.gff3",
    create_weights = weights[config["Outputs"]["geneid_preEVM"]],
  log:
    logs_dir + str(date) + ".j%j.Geneid.out",
    logs_dir + str(date) + ".j%j.Geneid.err",
  benchmark:
    logs_dir + str(date) + ".Geneid.benchmark.txt" 
  threads: 2  

use rule genemark from gene_predictors_workflow with:
  input:
    fasta = masked_reference,
  output:
    out = config["Outputs"]["genemark_prediction"],
    EVM_out = config["Outputs"]["genemark_preEVM"]
  params:
    scripts_dir = scripts_dir, 
    maxgap = config["genemark"]["max_gap"], 
    mincontig = config["genemark"]["min_contig"], 
    maxcontig = config["genemark"]["max_contig"], 
    add_opts =  config["genemark"]["additional_genemark_options"],
    EVM_dir = EVM_dir,
    link_out = EVM_dir + "genemark_predictions.gff3",
    create_weights_gmk = weights[config["Outputs"]["genemark_preEVM"]],
    rmcmd = "cd ../..; rm -r $TMPDIR;"
  envmodules:
    "GeneMark-ET"
  log:
    logs_dir + str(date) + ".j%j.genemark.out",
    logs_dir + str(date) + ".j%j.genemark.err"
  benchmark:
    logs_dir + str(date) + ".genemark.benchmark.txt"
  threads:  config["genemark"]["cores"],

use rule get_chunks_fasta from preprocess_workflow with:
  input:
    fasta = masked_reference
  output:
    expand(dirMaskedChunks + name_masked + ".{i}.fa", i=range(1, masked_chunks+1)),
  params:
    numberchunks = masked_chunks,
    scripts_dir = scripts_dir,
    dirChunks = dirMaskedChunks
  conda:
    '../envs/ann_base.yaml'
  log:
    logs_dir + str(date) + ".j%j.getChunks." + name_masked + ".out",
    logs_dir + str(date) + ".j%j.getChunks." + name_masked + ".err",
  benchmark:
    logs_dir + str(date) + ".getChunks." + name_masked + ".benchmark.txt"
  threads: 2

use rule augustus from gene_predictors_workflow with:
  input:
    fasta = masked_reference,
    trained = dirAugustusOutput + config["augustus"]["species"] + ".trained" if  config["augustus"]["species"] else []
  output:
    predictions = config["Outputs"]["augustus_prediction"]
  params:
    config_path = config["augustus"]["config_path"],
    species = config["augustus"]["species"],
    alternatives_from_sampling =  config["augustus"]["alternatives_from_sampling"],
    alternatives_from_evidence =  config["augustus"]["alternatives_from_evidence"],
    uniqueGeneId =  config["augustus"]["uniqueGeneId"],
    gff3 =  config["augustus"]["gff3"],
    sample =  config["augustus"]["sample"],
    noInFrameStop =  config["augustus"]["noInFrameStop"],
    maxtracks =  config["augustus"]["maxtracks"],
    singlestrand =  config["augustus"]["singlestrand"],
    strand =  config["augustus"]["strand"],
    min_intron_len =  config["augustus"]["min_intron_len"],
    additional_aug_opts = config["augustus"]["additional_augustus_options"],
  conda: 
    "../envs/augustus3.5.0.yaml"
  log:
    dirAugustusOutput + "logs/" + str(date) + ".j%j.Augustus.out",
    dirAugustusOutput + "logs/" + str(date) + ".j%j.Augustus.err",
  benchmark:
    dirAugustusOutput + "logs/" + str(date) + ".Augustus.benchmark.txt"      
  threads: 1

use rule augustus_jobarray from gene_predictors_workflow with:
  input:
    fasta = dirMaskedChunks + name_masked + "." + str(masked_chunks) + ".fa",
    trained = dirAugustusOutput + config["augustus"]["species"] + ".trained" if  config["augustus"]["species"] else []
  output:
    touch = dirAugustusOutput + "augustus_run.ok"
  params:
    config_path = config["augustus"]["config_path"],
    prefix = dirMaskedChunks + name_masked,
    out_prediction_augustus = config["Outputs"]["augustus_prediction"],
    species = config["augustus"]["species"],
    alternatives_from_sampling =  config["augustus"]["alternatives_from_sampling"],
    alternatives_from_evidence =  config["augustus"]["alternatives_from_evidence"],
    uniqueGeneId =  config["augustus"]["uniqueGeneId"],
    gff3 =  config["augustus"]["gff3"],
    sample =  config["augustus"]["sample"],
    noInFrameStop =  config["augustus"]["noInFrameStop"],
    maxtracks =  config["augustus"]["maxtracks"],
    singlestrand =  config["augustus"]["singlestrand"],
    strand =  config["augustus"]["strand"],
    min_intron_len =  config["augustus"]["min_intron_len"],
    additional_aug_opts = config["augustus"]["additional_augustus_options"],
  conda: 
    "../envs/augustus3.5.0.yaml"
  threads: 1

use rule merge_gffs from preprocess_workflow with:
  input:
    touch = "{path}/{name}_run.ok"
  output:
    out = "{path}/{name}_gene_prediction.gff3"
  params:
    array_inputs = lambda wildcards: expand(wildcards.path + "/" + wildcards.name + "_gene_prediction.gff3.{i}", i=range(1, masked_chunks+1)),
  log:
    "{path}/logs/" + str(date) + ".j%j.merge.{name}.out",
    "{path}/logs/" + str(date) + ".j%j.merge.{name}.err",
  benchmark:
    "{path}/logs/" + str(date) + ".merge.{name}.benchmark.txt"
  threads: 1  

use rule predictions4EVM from preprocess_workflow with:
  input:
    gff = "{path}/{name}_gene_prediction.gff3"
  output:
    EVM_out = "{path}/{name}_preEVM.gff3"
  params:
    reformat_cmd = lambda wildcards: reformat[wildcards.path+ "/" + wildcards.name + "_preEVM.gff3"],
    create_weights = lambda wildcards: weights[wildcards.path + "/" + wildcards.name + "_preEVM.gff3"],
    EVM_dir = EVM_dir,
    link_out = EVM_dir + "{name}_predictions.gff3"
  conda:
    "../envs/evm2.1.yaml"
  log:
    "{path}/logs/" + str(date) + ".j%j.preEVM.{name}.out",
    "{path}/logs/" + str(date) + ".j%j.preEVM.{name}.err",
  benchmark:
    "{path}/logs/" + str(date) + ".preEVM.{name}.benchmark.txt"   
  threads: 1

use rule get_coding_junctions from combiners_workflow with:
  input:
    junctions = config["Outputs"]["junctions"], 
    coding = config["Outputs"]["geneid_prediction"] if config["geneid"]["parameters"] else config["Outputs"]["miniprot_gene"],
  output:
    coding_junctions = config["Outputs"]["annotation_basedir"] + "coding_junctions.sorted.gff3"
  conda: 
    "../envs/bedtools2.30.0.yaml"
  log:
    logs_dir + str(date) + ".j%j.getCodingHints.out",
    logs_dir + str(date) + ".j%j.getCodingHints.err",
  benchmark:
    logs_dir + str(date) + ".getCodingHints.benchmark.txt"   
  threads: 1 

use rule geneid_introns from gene_predictors_workflow with:
  input:
    fasta = masked_reference,
    geneid_parameters = config["geneid"]["parameters"] if config["geneid"]["parameters"] else [],
    junctions = config["Outputs"]["annotation_basedir"] + "coding_junctions.sorted.gff3"
  output:
    out = config["Outputs"]["geneid_introns_prediction"],
    EVM_out = config["Outputs"]["geneid_introns_preEVM"]
  params:
    scripts_dir = scripts_dir, 
    geneid_options = config["geneid_introns"]["options"],
    path =  config["geneid"]["path"],
    EVM_dir = EVM_dir,
    link_out = EVM_dir + "geneid_introns_predictions.gff3",
    create_weights = weights[config["Outputs"]["geneid_introns_preEVM"]],
    rmcmd = "cd ../..; rm -r $TMPDIR/geneid_introns"
  log:
    logs_dir + str(date) + ".j%j.Geneid_introns.out",
    logs_dir + str(date) + ".j%j.Geneid_introns.err",
  benchmark:
    logs_dir + str(date) + ".Geneid_introns.benchmark.txt" 
  threads: 2

use rule genemark_ET from gene_predictors_workflow with:
  input:
    fasta = masked_reference,
    hints = config["Outputs"]["annotation_basedir"] + "coding_junctions.sorted.gff3",
  output:
    out = config["Outputs"]["genemark_ET_prediction"],
    EVM_out = config["Outputs"]["genemark_ET_preEVM"]
  params:
    scripts_dir = scripts_dir, 
    maxgap = config["genemark"]["max_gap"], 
    mincontig = config["genemark"]["min_contig"], 
    maxcontig = config["genemark"]["max_contig"], 
    add_opts =  config["genemark_ET"]["additional_genemark_ET_options"],
    EVM_dir = EVM_dir,
    link_out = EVM_dir + "genemark-ET_predictions.gff3",
    create_weights_gmk = weights[config["Outputs"]["genemark_ET_preEVM"]],
    rmcmd = "cd ../..; rm -r $TMPDIR;"
  envmodules:
    "GeneMark-ET"
  log:
    logs_dir + str(date) + ".j%j.genemark-ET.out",
    logs_dir + str(date) + ".j%j.genemark-ET.err"
  benchmark:
    logs_dir + str(date) + ".genemark-ET.benchmark.txt"
  threads:  config["genemark"]["cores"]

use rule augustus_hints from gene_predictors_workflow with:
  input:
    fasta = masked_reference,
    hints = config["Outputs"]["annotation_basedir"] + "coding_junctions.sorted.gff3",
    extrinsic_file = config["augustus_hints"]["extrinsic_file_augustus_hints"],
    trained = dirAugustusOutput + config["augustus"]["species"] + ".trained" if  config["augustus"]["species"]  else []
  output:
    predictions =config["Outputs"]["augustus_hints_prediction"]
  params:
    scripts_dir = scripts_dir,
    config_path = config["augustus"]["config_path"],
    species = config["augustus"]["species"],
    alternatives_from_sampling =  config["augustus"]["alternatives_from_sampling"],
    alternatives_from_evidence =  config["augustus"]["alternatives_from_evidence"],
    uniqueGeneId =  config["augustus"]["uniqueGeneId"],
    gff3 =  config["augustus"]["gff3"],
    sample =  config["augustus"]["sample"],
    noInFrameStop =  config["augustus"]["noInFrameStop"],
    maxtracks =  config["augustus"]["maxtracks"],
    singlestrand =  config["augustus"]["singlestrand"],
    strand =  config["augustus"]["strand"],
    min_intron_len =  config["augustus"]["min_intron_len"],
    additional_aug_opts = config["augustus_hints"]["additional_augustus_hints_options"],
    rmcmd = "cd ../..; rm -r $TMPDIR/augustus_hints"
  conda: 
    "../envs/augustus3.5.0.yaml"
  log:
    dirAugustusHintsOutput + "logs/" + str(date) + ".j%j.Augustus.out",
    dirAugustusHintsOutput + "logs/" + str(date) + ".j%j.Augustus.err",
  benchmark:
    dirAugustusHintsOutput + "logs/" + str(date) + ".Augustus.benchmark.txt"      
  threads: 1

use rule augustus_hints_jobarray from gene_predictors_workflow with:
  input:
    fasta = dirMaskedChunks + name_masked + "." + str(masked_chunks) + ".fa",
    hints = config["Outputs"]["annotation_basedir"] + "coding_junctions.sorted.gff3",
    extrinsic_file = config["augustus_hints"]["extrinsic_file_augustus_hints"],
    trained = dirAugustusOutput + config["augustus"]["species"] + ".trained" if config["augustus"]["species"] else []
  output:
    touch = dirAugustusHintsOutput + "augustus_hints_run.ok",
  params:
    scripts_dir = scripts_dir,
    config_path = config["augustus"]["config_path"],
    prefix = dirMaskedChunks + name_masked,
    out_prediction_augustus = config["Outputs"]["augustus_hints_prediction"],
    species = config["augustus"]["species"],
    alternatives_from_sampling =  config["augustus"]["alternatives_from_sampling"],
    alternatives_from_evidence =  config["augustus"]["alternatives_from_evidence"],
    uniqueGeneId =  config["augustus"]["uniqueGeneId"],
    gff3 =  config["augustus"]["gff3"],
    sample =  config["augustus"]["sample"],
    noInFrameStop =  config["augustus"]["noInFrameStop"],
    maxtracks =  config["augustus"]["maxtracks"],
    singlestrand =  config["augustus"]["singlestrand"],
    strand =  config["augustus"]["strand"],
    min_intron_len =  config["augustus"]["min_intron_len"],
    additional_aug_opts = config["augustus_hints"]["additional_augustus_hints_options"],
    rmcmd = "cd ../..; rm -r $dir"
  conda: 
    "../envs/augustus3.5.0.yaml"
  threads: 1   

##Define rules priority
if masked_chunks > 1:
  ruleorder: merge_gffs > augustus
  ruleorder: merge_gffs > augustus_hints

