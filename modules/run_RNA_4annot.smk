from datetime import datetime
import os
import re
import subprocess

date = datetime.now().strftime('%Y%m%d.%H%M%S')

##Include modules:
module align_reads_workflow:
  snakefile: "../modules/align_reads.smk"

module preprocess_workflow:
  snakefile: "../modules/process_reads.rules.smk"

module get_models_workflow:
  snakefile: "../modules/get_models.smk"

# file wildcards 
illumina_files = config["Wildcards"]["illumina_fastqs"]
cDNA_files = config["Wildcards"]["cDNA_fastqs"]
dRNA_files = config["Wildcards"]["dRNA_fastqs"]
pb_files = config["Wildcards"]["isoseq_fastqs"]
pb_fasta_files = config["Wildcards"]["isoseq_fastas"]

##PREPARE RNA RUN
rna_dir = config["Outputs"]["RNA_outdir"]
rnalogs_dir = rna_dir + "logs/"
if not os.path.exists(rnalogs_dir):
  os.makedirs(rnalogs_dir)

illumina_processed = ""
starlogs_dir = config["Illumina RNA"]["star_dir"] + "logs/"
stringtie_in = {}
extensions = {}
taco_models_in = {}
taco_opts = {}
sams_list = []
jbrowse_targets={}
jbrowse_targets["RNA"] = []

taco_models_in[config["Outputs"]["TACO_dir"]] = []
taco_opts[config["Outputs"]["TACO_dir"]] = config["Model RNA"]["TACO global options"]
if illumina_files != None:
  illumina_processed= config["Illumina RNA"]["star_dir"] + "trimmed/"
  if not os.path.exists(illumina_processed):
    os.makedirs(illumina_processed)
  taco_models_in[config["Illumina RNA"]["star_dir"]] = []
  taco_models_in[config["Outputs"]["TACO_dir"]].append(config["Illumina RNA"]["star_dir"] + "TACO_assembled.gtf")
  for name in illumina_files.split(','):
    stringtie_in[config["Illumina RNA"]["star_dir"] + name + ".stringtie.gtf"] = config["Illumina RNA"]["stringtie_opts"]
    extensions[config["Illumina RNA"]["star_dir"] + name] = "Aligned.sortedByCoord.out.bam"
    taco_models_in[config["Illumina RNA"]["star_dir"] ].append(config["Illumina RNA"]["star_dir"] + name + ".stringtie.gtf")
    sams_list.append(config["Illumina RNA"]["star_dir"] + name + "Aligned.sortedByCoord.out.sam")
    jbrowse_targets["RNA"].append(jbrowse_dir + "RNA/" + os.path.basename(os.path.dirname(config["Illumina RNA"]["star_dir"])) + "/" + name + "Aligned.sortedByCoord.out.bw")
  taco_opts[config["Illumina RNA"]["star_dir"]] = config["Illumina RNA"]["TACO_opts"]
  if not os.path.exists(starlogs_dir):
    os.makedirs(starlogs_dir)
  
minimap_in = {}
minimap_opts = {}
if cDNA_files != None:
  taco_models_in[config["cDNA RNA"]["minimap_dir"] ] = []
  taco_models_in[config["Outputs"]["TACO_dir"]].append(config["cDNA RNA"]["minimap_dir"] + "TACO_assembled.gtf")
  for name in cDNA_files.split(','):
    stringtie_in[config["cDNA RNA"]["minimap_dir"] + name + ".stringtie.gtf"] = config["cDNA RNA"]["stringtie_opts"]
    extensions[config["cDNA RNA"]["minimap_dir"] + name] = ".sorted.bam"
    minimap_in[config["cDNA RNA"]["minimap_dir"] + name + ".sorted.bam"] = config["Inputs"]["cDNA_dir"] + name + ".fastq.gz"
    minimap_opts[config["cDNA RNA"]["minimap_dir"] + name] =  config["cDNA RNA"]["minimap2_opts"]
    taco_models_in[config["cDNA RNA"]["minimap_dir"]].append(config["cDNA RNA"]["minimap_dir"] + name + ".stringtie.gtf")
    sams_list.append(config["cDNA RNA"]["minimap_dir"] + name + ".sorted.sam")
    jbrowse_targets["RNA"].append(jbrowse_dir + "RNA/" + os.path.basename(os.path.dirname(config["cDNA RNA"]["minimap_dir"])) + "/" + name + ".sorted.bw")
  taco_opts[config["cDNA RNA"]["minimap_dir"] ] = config["cDNA RNA"]["TACO_opts"],
  if not os.path.exists(config["cDNA RNA"]["minimap_dir"] + "/logs"):
    os.makedirs(config["cDNA RNA"]["minimap_dir"] + "/logs")
  
if dRNA_files != None:
  taco_models_in[config["dRNA RNA"]["minimap_dir"] ] = []
  taco_models_in[config["Outputs"]["TACO_dir"]].append(config["dRNA RNA"]["minimap_dir"] + "TACO_assembled.gtf")
  for name in dRNA_files.split(','):
    stringtie_in[config["dRNA RNA"]["minimap_dir"] + name + ".stringtie.gtf"] = config["dRNA RNA"]["stringtie_opts"]
    extensions[config["dRNA RNA"]["minimap_dir"] + name] = ".sorted.bam"
    minimap_in[config["dRNA RNA"]["minimap_dir"] + name + ".sorted.bam"] = config["Inputs"]["dRNA_dir"] + name + ".fastq.gz"
    minimap_opts[config["dRNA RNA"]["minimap_dir"] + name] =  config["dRNA RNA"]["minimap2_opts"]
    taco_models_in[config["dRNA RNA"]["minimap_dir"]].append(config["dRNA RNA"]["minimap_dir"] + name + ".stringtie.gtf")
    sams_list.append(config["dRNA RNA"]["minimap_dir"] + name + ".sorted.sam")
    jbrowse_targets["RNA"].append(jbrowse_dir + "RNA/" + os.path.basename(os.path.dirname(config["dRNA RNA"]["minimap_dir"])) + "/" + name + ".sorted.bw")
  taco_opts[config["dRNA RNA"]["minimap_dir"] ] = config["dRNA RNA"]["TACO_opts"],
  if not os.path.exists(config["dRNA RNA"]["minimap_dir"] + "/logs"):
    os.makedirs(config["dRNA RNA"]["minimap_dir"] + "/logs")
  
if pb_files != None or pb_fasta_files != None:
  taco_models_in[config["Isoseq"]["minimap_dir"] ] = []
  taco_models_in[config["Outputs"]["TACO_dir"]].append(config["Isoseq"]["minimap_dir"] + "TACO_assembled.gtf")
  taco_opts[config["Isoseq"]["minimap_dir"] ] = config["Isoseq"]["TACO_opts"],
  if not os.path.exists(config["Isoseq"]["minimap_dir"] + "/logs"):
    os.makedirs(config["Isoseq"]["minimap_dir"] + "/logs")
  name = []
  if pb_files != None:
     for i in pb_files.split(','):
      name.append(i)
      minimap_in[config["Isoseq"]["minimap_dir"] + i + ".sorted.bam"] = config["Inputs"]["PB_dir"] + i + ".fastq.gz"
  if pb_fasta_files != None:
    for i in pb_fasta_files.split(','):
      name.append(i)
      minimap_in[config["Isoseq"]["minimap_dir"] + i + ".sorted.bam"] = config["Inputs"]["PB_dir"] + i + ".fasta"
  for name in name:
    stringtie_in[config["Isoseq"]["minimap_dir"] + name + ".stringtie.gtf"] = config["Isoseq"]["stringtie_opts"]
    extensions[config["Isoseq"]["minimap_dir"] + name] = ".sorted.bam"
    minimap_opts[config["Isoseq"]["minimap_dir"] + name] =  config["Isoseq"]["minimap2_opts"]
    taco_models_in[config["Isoseq"]["minimap_dir"]].append(config["Isoseq"]["minimap_dir"] + name + ".stringtie.gtf")
    sams_list.append(config["Isoseq"]["minimap_dir"] + name + ".sorted.sam")
    jbrowse_targets["RNA"].append(jbrowse_dir + "RNA/" + os.path.basename(os.path.dirname(config["Isoseq"]["minimap_dir"])) + "/" + name + ".sorted.bw")

##RUN RULES
use rule star_index from align_reads_workflow with: 
  input:
    genome = config["Inputs"]["genome"]
  output: 
    genomedir = directory(config["Illumina RNA"]["star_genome_dir"]),
    ok = rna_dir + "index.ok"
  params:
    opts = config["Illumina RNA"]["star_index_additional_opts"]
  log:
    rnalogs_dir + str(date) + ".star_index.j%j.out",
    rnalogs_dir + str(date) + ".star_index.j%j.err",
  benchmark:
    rnalogs_dir + str(date) + ".star_index.benchmark.txt"
  conda:
    "../envs/star2.7.10a.yaml"
  threads:  config["Illumina RNA"]["starCores"]

use rule trim_galore from preprocess_workflow with:
  input:
    read1 = lambda wildcards: config["Inputs"]["illumina_dir"]  + "{file}.1.fastq.gz" \
            if config["Illumina RNA"]["PE"] == True else \
            lambda wildcards: config["Inputs"]["illumina_dir"]  + "{file}.fastq.gz",
    read2 = lambda wildcards: config["Inputs"]["illumina_dir"]  + "{file}.2.fastq.gz" \
            if config["Illumina RNA"]["PE"] == True else ""
  output:
    trim1 = illumina_processed + "{file}.trimmed.1.fastq.gz",
    trim2 = illumina_processed + "{file}.trimmed.2.fastq.gz" if config["Illumina RNA"]["PE"] == True else None,
  params:
    outdir = illumina_processed,
    opts = config["Trim_Galore"]["options"] 
  log:
    starlogs_dir + str(date) + ".{file}.trim_galore.j%j.out",
    starlogs_dir + str(date) + ".{file}.trim_galore.j%j.err",
  benchmark:
    starlogs_dir + str(date) + ".{file}.trim_galore.benchmark.txt",
  conda:
    "../envs/trim_galore0.6.7.yaml"
  threads: config["Trim_Galore"]["Trim_Illumina_cores"]

use rule star from align_reads_workflow with:
  input:
    genomedir = config["Illumina RNA"]["star_genome_dir"],
    reads = [illumina_processed + "{file}.trimmed.1.fastq.gz",  illumina_processed + "{file}.trimmed.2.fastq.gz"]
            if config["Illumina RNA"]["PE"] == True else lambda wildcards:[ illumina_processed  + "{file}.trimmed.fastq.gz" ]
  output:
    bam = config["Illumina RNA"]["star_dir"] + "{file}Aligned.sortedByCoord.out.bam",
  params:
    stardir = config["Illumina RNA"]["star_dir"],
    basename = "{file}",
    additional = config["Illumina RNA"]["star_additional_opts"]
  log:
    starlogs_dir + str(date) + ".{file}.star.j%j.out",
    starlogs_dir + str(date) + ".{file}.star.j%j.err",
  benchmark:
    starlogs_dir + str(date) + ".{file}.star.benchmark.txt",
  conda:
    "../envs/star2.7.10a.yaml"
  threads:  config["Illumina RNA"]["starCores"]

use rule minimap2 from align_reads_workflow with:
  input:
    genome = config["Inputs"]["genome"],
    reads = lambda wildcards: minimap_in[wildcards.dir + "/" + wildcards.file + ".sorted.bam"],
  output:
    bam = "{dir}/{file}.sorted.bam",
  params:
    basename = "{file}",
    minimap_opts = lambda wildcards:  minimap_opts[wildcards.dir + "/" + wildcards.file]
  log:
    "{dir}/logs/" + str(date) + ".{file}.minimap2.j%j.out",
    "{dir}/logs/" + str(date) + ".{file}.minimap2.j%j.err"
  benchmark:
    "{dir}/logs/" + str(date) + ".{file}.minimap2.benchmark.txt",
  conda:
    "../envs/minimap2.24.yaml"    
  threads: config["Model RNA"]["minimapCores"]

use rule stringtie from get_models_workflow with:
  input:
    bam = lambda wildcards: wildcards.dir + "/" + wildcards.file + extensions[wildcards.dir + "/" + wildcards.file],
  output:
    models = "{dir}/{file}" + ".stringtie.gtf"
  params:
    basename = "{file}",
    stringtie_opts = lambda wildcards: stringtie_in[wildcards.dir + "/" +  wildcards.file + ".stringtie.gtf"],
  threads:
    config["Illumina RNA"]["starCores"]
  log:
    "{dir}/logs/" + str(date) + ".{file}.stringtie.j%j.out",
    "{dir}/logs/" + str(date) + ".{file}.stringtie.j%j.err",
  benchmark:
    "{dir}/logs/" + str(date) + ".{file}.stringtie.benchmark.txt",
  conda:
    "../envs/stringtie2.2.1.yaml"

use rule join_models from get_models_workflow with:
  input:
    genome = config["Inputs"]["genome"],
    models = lambda wildcards: taco_models_in[wildcards.dir]
  output:
    joined_models = "{dir}TACO_assembled.gtf",
  params:
    TACO_opts = lambda wildcards: taco_opts[wildcards.dir],
    outdir = "{dir}"
  threads: 
    config["Model RNA"]["TACOCores"]
  log:
    "{dir}logs/" + str(date) + ".TACO.j%j.out",
    "{dir}logs/" + str(date) + ".TACO.j%j.err",
  benchmark:
    "{dir}logs/" + str(date) + ".TACO.benchmark.txt",
  conda:
    "../envs/taco0.7.3.yaml"

use rule bam2sam from align_reads_workflow with:
  input:
    bam = "{dir}/{file}.bam",
  output:
    sam = "{dir}/{file}.sam"
  log:
    "{dir}/logs/" + str(date) + ".{file}.bam2sam.j%j.out",
    "{dir}/logs/" + str(date) + ".{file}.bam2sam.j%j.err",
  benchmark:
    "{dir}/logs/" + str(date) + ".{file}.bam2sam.benchmark.txt",
  conda:
    "../envs/ESPRESSO1.3.0.yaml"
  threads:
    config["Model RNA"]["EspressoCores"]

use rule ESPRESSO from get_models_workflow with:
  input:
    genome = config["Inputs"]["genome"],
    sams =  sams_list,
    tsv = config["Inputs"]["RNA_Samples_tsv"],
  output:
    junctions = config["Outputs"]["junctions"]
  params:
    ESPRESSO_path = config["Model RNA"]["ESPRESSO path"],
    out_dir = config["Model RNA"]["ESPRESSO outdir"]
  log:
    rnalogs_dir + str(date) + ".espresso.j%j.out",
    rnalogs_dir + str(date) + ".espresso.j%j.err",
  benchmark:
    rnalogs_dir + str(date) + ".espresso.benchmark.txt",
  conda:
    "../envs/ESPRESSO1.3.0.yaml"
  threads: 
    config["Model RNA"]["EspressoCores"]
