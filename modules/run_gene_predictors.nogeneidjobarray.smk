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
predictions_list = {}
predictions_links = {}
weights = {}
reformat = {}

if config["Parameters"]["run_augustus"]:
  out_prediction_augustus = config["Outputs"]["augustus_prediction"]
  dirAugustusOutput = os.path.dirname(out_prediction_augustus)
  if not os.path.exists(dirAugustusOutput + "/logs/"):
    os.makedirs(dirAugustusOutput + "/logs/")

  augustus_weights = config["augustus"]["augustus_weights"]
  create_weights_aug = ""
  f = 1
  for i in range(len(augustus_weights)):
    create_weights_aug = create_weights_aug + "echo \"ABINITIO_PREDICTION\tAugustus\t" + str(augustus_weights[i]) + "\">> weights_" + str(f) + ".txt;";      
    f = f + 1
  weights[config["Outputs"]["augustus_preEVM"]] = create_weights_aug
  reformat[config["Outputs"]["augustus_preEVM"]] = "$CONDA_PREFIX/opt/evidencemodeler-1.1.1/EvmUtils/misc/augustus_GFF3_to_EVM_GFF3.pl " + out_prediction_augustus

  additional_aug_opts =  config["augustus"]["additional_augustus_options"]
  if additional_aug_opts == None:
    additional_aug_opts = ""

if config["Parameters"]["run_augustus_hints"]:
  out_prediction_augustus_hints = config["Outputs"]["augustus_hints_prediction"]
  dirAugustusHintsOutput = os.path.dirname(out_prediction_augustus_hints )
  if not os.path.exists(dirAugustusHintsOutput + "/logs/"):
    os.makedirs(dirAugustusHintsOutput + "/logs")
  
  augustus_hints_weights = config["augustus_hints"]["augustus_hints_weights"]
  create_weights_aug_hints = ""
  f = 1
  for i in range(len(augustus_hints_weights)):
    create_weights_aug_hints = create_weights_aug_hints + "echo \"ABINITIO_PREDICTION\taugustus_hints\t" + str(augustus_hints_weights[i]) + "\">> weights_" + str(f) + ".txt;";      
    f = f + 1
  weights[config["Outputs"]["augustus_hints_preEVM"]] = create_weights_aug_hints
  reformat[config["Outputs"]["augustus_hints_preEVM"]] = "$CONDA_PREFIX/opt/evidencemodeler-1.1.1/EvmUtils/misc/augustus_GFF3_to_EVM_GFF3.pl " + out_prediction_augustus_hints+ " | sed 's/Augustus/augustus_hints/g'"

  additional_aug_hints_opts =  config["augustus_hints"]["additional_augustus_hints_options"]
  if additional_aug_hints_opts == None:
    additional_aug_hints_opts = ""
  additional_aug_hints_opts += " --softmasking=0"

if config["Parameters"]["run_geneid"]:
  out_prediction_geneid = config["Outputs"]["geneid_prediction"]
  dirGeneidOutput = os.path.dirname(out_prediction_geneid)
  if not os.path.exists(dirGeneidOutput + "/logs/"):
    os.makedirs(dirGeneidOutput + "/logs/")

  geneid_weights = config["geneid"]["geneid_weights"]
  create_weights_geneid = ""
  f = 1
  for i in range(len(geneid_weights)):
    create_weights_geneid = create_weights_geneid + "echo \"ABINITIO_PREDICTION\tgeneid_v1.4\t" + str(geneid_weights[i]) + "\">> weights_" + str(f) + ".txt;";      
    f = f + 1
  weights[config["Outputs"]["geneid_preEVM"]] = create_weights_geneid
  reformat[config["Outputs"]["geneid_preEVM"]] = "grep -v '#' " + out_prediction_geneid

  geneid_opts =  config["geneid"]["geneid_options"]
  if geneid_opts == None:
    geneid_opts = ""

if config["Parameters"]["run_geneid_introns"]:
  out_prediction_geneid_introns = config["Outputs"]["geneid_introns_prediction"]
  dirGeneidOutputIntrons = os.path.dirname(out_prediction_geneid_introns)
  if not os.path.exists(dirGeneidOutputIntrons + "/logs/"):
    os.makedirs(dirGeneidOutputIntrons + "/logs/")

  geneid_introns_weights = config["geneid_introns"]["geneid_introns_weights"]
  create_weights_geneid_introns = ""
  f = 1
  for i in range(len(geneid_introns_weights)):
    create_weights_geneid_introns = create_weights_geneid_introns + "echo \"ABINITIO_PREDICTION\tgeneid_introns\t" + str(geneid_introns_weights[i]) + "\">> weights_" + str(f) + ".txt;";      
    f = f + 1
  weights[config["Outputs"]["geneid_introns_preEVM"]] = create_weights_geneid_introns
  reformat[config["Outputs"]["geneid_introns_preEVM"]] = "grep -v '#' " + out_prediction_geneid_introns + " | sed 's/geneid_v1.4/geneid_introns/g' "

  geneid_introns_opts =  config["geneid_introns"]["geneid_introns_options"]
  if geneid_introns_opts == None:
    geneid_introns_opts = ""

if config["Parameters"]["run_genemark"]:
  out_prediction_genemark = config["Outputs"]["genemark_prediction"]

  genemark_weights = config["genemark"]["genemark_weights"]
  create_weights_gmk = ""
  f = 1
  for i in range(len(genemark_weights)):
    create_weights_gmk = create_weights_gmk + "echo \"ABINITIO_PREDICTION\tGeneMark.hmm3\t" + str(genemark_weights[i]) + "\">> weights_" + str(f) + ".txt;";      
    f = f + 1
  weights[config["Outputs"]["genemark_preEVM"]] = create_weights_gmk

  additional_gmk_opts =  config["genemark"]["additional_genemark_options"]
  if additional_gmk_opts == None:
    additional_gmk_opts = ""

if config["Parameters"]["run_genemark-ET"]:
  genemark_et_weights = config["genemark-ET"]["genemark_ET_weights"]
  create_weights_gmk_et = ""
  f = 1
  for i in range(len(genemark_et_weights)):
    create_weights_gmk_et = create_weights_gmk_et + "echo \"ABINITIO_PREDICTION\tGeneMark-ET\t" + str(genemark_et_weights[i]) + "\">> weights_" + str(f) + ".txt;";      
    f = f + 1
  weights[config["Outputs"]["genemark_ET_preEVM"]] = create_weights_gmk_et
  
  additional_gmk_et_opts =  config["genemark-ET"]["additional_genemark_ET_options"]
  if additional_gmk_et_opts == None:
    additional_gmk_et_opts = ""

#Execute rules:
if config["Parameters"]["run_geneid_introns"] or config["Parameters"]["run_augustus_hints"] or config["Parameters"]["run_genemark-ET"]:
  use rule get_coding_junctions from combiners_workflow with:
    input:
      junctions = config["Inputs"]["junctions"], 
      coding = out_prediction_geneid
    output:
      coding_junctions = config["Outputs"]["output_dir"] + "coding_junctions.sorted.gff3"
    conda: 
      "../envs/bedtools2.30.0.yaml"
    log:
      logs_dir + str(date) + ".j%j.getCodingHints.out",
      logs_dir + str(date) + ".j%j.getCodingHints.err",
    benchmark:
      logs_dir + str(date) + ".getCodingHints.benchmark.txt"   
    threads: 1 

if config["Parameters"]["run_augustus"] or config["Parameters"]["run_augustus_hints"]:
  if masked_chunks > 1:
    dirMaskedChunks = config["Outputs"]["dir_masked_chunks"] 
    if not os.path.exists(dirMaskedChunks):
      os.makedirs(dirMaskedChunks)
   # print (dirMaskedChunks + name_masked)
    use rule get_chunks_fasta from preprocess_workflow with:
      input:
        fasta = masked_reference
      output:
        expand(dirMaskedChunks + name_masked + ".{i}.fa", i=range(1, masked_chunks+1)),
       # "hello.txt"
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

    if config["Parameters"]["run_augustus"]:
      use rule augustus_jobarray from gene_predictors_workflow with:
        input:
          fasta = dirMaskedChunks + name_masked + "." + str(masked_chunks) + ".fa"
        output:
          touch = dirAugustusOutput + "/augustus_run.ok"
        params:
          prefix = dirMaskedChunks + name_masked,
          out_prediction_augustus = out_prediction_augustus,
          species = config["augustus"]["aug_species"],
          alternatives_from_sampling =  config["augustus"]["aug_alternatives_from_sampling"],
          alternatives_from_evidence =  config["augustus"]["aug_alternatives_from_evidence"],
          uniqueGeneId =  config["augustus"]["aug_uniqueGeneId"],
          gff3 =  config["augustus"]["aug_gff3"],
          sample =  config["augustus"]["aug_sample"],
          noInFrameStop =  config["augustus"]["aug_noInFrameStop"],
          maxtracks =  config["augustus"]["aug_maxtracks"],
          singlestrand =  config["augustus"]["aug_singlestrand"],
          strand =  config["augustus"]["aug_strand"],
          min_intron_len =  config["augustus"]["aug_min_intron_len"],
          additional_aug_opts = additional_aug_opts
        conda: 
          "../envs/augustus3.5.0.yaml"
        threads: 1

    if config["Parameters"]["run_augustus_hints"]:
      use rule augustus_hints_jobarray from gene_predictors_workflow with:
        input:
          fasta = dirMaskedChunks + name_masked + "." + str(masked_chunks) + ".fa",
          hints = config["Outputs"]["output_dir"] + "coding_junctions.sorted.gff3",
          extrinsic_file = config["augustus_hints"]["extrinsic_file_augustus_hints"]
        output:
          touch = dirAugustusHintsOutput + "/augustus_hints_run.ok",
        params:
          scripts_dir = scripts_dir,
          prefix = dirMaskedChunks + name_masked,
          out_prediction_augustus = out_prediction_augustus_hints,
          species = config["augustus"]["aug_species"],
          alternatives_from_sampling =  config["augustus"]["aug_alternatives_from_sampling"],
          alternatives_from_evidence =  config["augustus"]["aug_alternatives_from_evidence"],
          uniqueGeneId =  config["augustus"]["aug_uniqueGeneId"],
          gff3 =  config["augustus"]["aug_gff3"],
          sample =  config["augustus"]["aug_sample"],
          noInFrameStop =  config["augustus"]["aug_noInFrameStop"],
          maxtracks =  config["augustus"]["aug_maxtracks"],
          singlestrand =  config["augustus"]["aug_singlestrand"],
          strand =  config["augustus"]["aug_strand"],
          min_intron_len =  config["augustus"]["aug_min_intron_len"],
          additional_aug_opts = additional_aug_hints_opts,
          rmcmd = "cd ../..; rm -r $dir"
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

  else:
    if config["Parameters"]["run_augustus"]:
      use rule augustus from gene_predictors_workflow with:
        input:
          fasta = masked_reference
        output:
          predictions = out_prediction_augustus
        params:
          species = config["augustus"]["aug_species"],
          alternatives_from_sampling =  config["augustus"]["aug_alternatives_from_sampling"],
          alternatives_from_evidence =  config["augustus"]["aug_alternatives_from_evidence"],
          uniqueGeneId =  config["augustus"]["aug_uniqueGeneId"],
          gff3 =  config["augustus"]["aug_gff3"],
          sample =  config["augustus"]["aug_sample"],
          noInFrameStop =  config["augustus"]["aug_noInFrameStop"],
          maxtracks =  config["augustus"]["aug_maxtracks"],
          singlestrand =  config["augustus"]["aug_singlestrand"],
          strand =  config["augustus"]["aug_strand"],
          min_intron_len =  config["augustus"]["aug_min_intron_len"],
          additional_aug_opts = additional_aug_opts,
          rmcmd = "cd ../..; rm -r $TMPDIR/augustus_hints"
        conda: 
          "../envs/augustus3.5.0.yaml"
        log:
          logs_dir + str(date) + ".j%j.Augustus.out",
          logs_dir + str(date) + ".j%j.Augustus.err",
        benchmark:
          logs_dir + str(date) + ".Augustus.benchmark.txt"      
        threads: 1

    if config["Parameters"]["run_augustus_hints"]:
      use rule augustus_hints from gene_predictors_workflow with:
        input:
          fasta = masked_reference,
          hints= config["Outputs"]["output_dir"] + "coding_junctions.sorted.gff3",
          extrinsic_file = config["augustus_hints"]["extrinsic_file_augustus_hints"]
        output:
          predictions = out_prediction_augustus_hints
        params:
          scripts_dir = scripts_dir,
          species = config["augustus"]["aug_species"],
          alternatives_from_sampling =  config["augustus"]["aug_alternatives_from_sampling"],
          alternatives_from_evidence =  config["augustus"]["aug_alternatives_from_evidence"],
          uniqueGeneId =  config["augustus"]["aug_uniqueGeneId"],
          gff3 =  config["augustus"]["aug_gff3"],
          sample =  config["augustus"]["aug_sample"],
          noInFrameStop =  config["augustus"]["aug_noInFrameStop"],
          maxtracks =  config["augustus"]["aug_maxtracks"],
          singlestrand =  config["augustus"]["aug_singlestrand"],
          strand =  config["augustus"]["aug_strand"],
          min_intron_len =  config["augustus"]["aug_min_intron_len"],
          additional_aug_opts = additional_aug_opts,
        conda: 
          "../envs/augustus3.5.0.yaml"
        log:
          logs_dir + str(date) + ".j%j.AugustusHints.out",
          logs_dir + str(date) + ".j%j.AugustusHints.err",
        benchmark:
          logs_dir + str(date) + ".AugustusHints.benchmark.txt"      
        threads: 1
    
  if config["Parameters"]["run_geneid"]:
    use rule geneid from gene_predictors_workflow with:
      input:
        fasta = masked_reference,
        geneid_parameters = config["Inputs"]["geneid_parameters"]
      output:
        out_prediction_geneid = out_prediction_geneid
      params:
        geneid_options = geneid_opts,
        path =  config["geneid"]["geneid_path"]
      log:
        logs_dir + str(date) + ".j%j.Geneid.out",
        logs_dir + str(date) + ".j%j.Geneid.err",
      benchmark:
        logs_dir + str(date) + ".Geneid.benchmark.txt" 
      threads: 2   

  if config["Parameters"]["run_geneid_introns"]:
    use rule geneid_introns from gene_predictors_workflow with:
      input:
        fasta = masked_reference,
        geneid_parameters = config["Inputs"]["geneid_parameters"],
        junctions = config["Outputs"]["output_dir"] + "coding_junctions.sorted.gff3",
      output:
        predictions = out_prediction_geneid_introns
      params:
        scripts_dir = scripts_dir, 
        geneid_options = geneid_introns_opts,
        path =  config["geneid"]["geneid_path"],
        rmcmd = "cd ../..; rm -r $TMPDIR/geneid_introns"
      log:
        logs_dir + str(date) + ".j%j.Geneid_introns.out",
        logs_dir + str(date) + ".j%j.Geneid_introns.err",
      benchmark:
        logs_dir + str(date) + ".Geneid_introns.benchmark.txt" 
      threads: 2

  if config["Parameters"]["run_genemark"]:
    use rule genemark from gene_predictors_workflow with:
      input:
        fasta = masked_reference,
      output:
        out = config["Outputs"]["genemark_prediction"],
        EVM_out = config["Outputs"]["genemark_preEVM"]
      params:
        scripts_dir = scripts_dir, 
        maxgap = config["genemark"]["gmk_max_gap"], 
        mincontig = config["genemark"]["gmk_min_contig"], 
        maxcontig = config["genemark"]["gmk_max_contig"], 
        add_opts = additional_gmk_opts,
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
      threads:  config["genemark"]["gmk_cores"],
      
  if config["Parameters"]["run_genemark-ET"]:
    use rule genemark_ET from gene_predictors_workflow with:
      input:
        fasta = masked_reference,
        hints = config["Outputs"]["output_dir"] + "coding_junctions.sorted.gff3",
      output:
        out = config["Outputs"]["genemark_ET_prediction"],
        EVM_out = config["Outputs"]["genemark_ET_preEVM"]
      params:
        scripts_dir = scripts_dir, 
        maxgap = config["genemark"]["gmk_max_gap"], 
        mincontig = config["genemark"]["gmk_min_contig"], 
        maxcontig = config["genemark"]["gmk_max_contig"], 
        add_opts = additional_gmk_et_opts,
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
      threads:  config["genemark"]["gmk_cores"],
  
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
      "../envs/evm.yaml"
    log:
      "{path}/logs/" + str(date) + ".j%j.preEVM.{name}.out",
      "{path}/logs/" + str(date) + ".j%j.preEVM.{name}.err",
    benchmark:
      "{path}/logs/" + str(date) + ".preEVM.{name}.benchmark.txt"   
    threads: 1

