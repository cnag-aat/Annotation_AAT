from datetime import datetime
import os
import re
import subprocess

date = datetime.now().strftime('%Y%m%d.%H%M%S')

##Include modules:
module repeats_workflow:
  snakefile: "../modules/repeats.rules.smk"

##PREPARE REPEAT ANNOTATION RUN
repeats_dir = config["Outputs"]["Repeat_Annotation_dir"]
rlogs_dir = repeats_dir + "logs/"
if not os.path.exists(rlogs_dir):
  os.makedirs(rlogs_dir)

rruns = 0
if config["Repeat Annotation"]["species_repeat_database"]:
  rruns += 1
if config["Repeat Annotation"]["Repeat_library"]:
  rruns +=1
if config["Parameters"]["run_redmask"]:
  rruns +=1

rmasked = {}
rdb = {}
rm_input_genome = config["Inputs"]["genome"]
repeat_masker_out_list = []

if config["Repeat Annotation"]["species_repeat_database"]:
  rspecies = config["Repeat Annotation"]["species_repeat_database"]
  if ' '  in rspecies:
    rspecies = re.sub('\s+', '_', rspecies) 
  
  if rruns == 1:
    repeat_base = base
    rmasked[repeat_base] = rm_input_genome
    rdb[repeat_base]= " -species " + config["Repeat Annotation"]["species_repeat_database"]
  else:
    repeat_base = base + ".rmask_" + rspecies
    rmasked[repeat_base] = rm_input_genome
    rdb[repeat_base]= " -species " + config["Repeat Annotation"]["species_repeat_database"]
    rm_input_genome = repeats_dir + repeat_base + ".masked.fa"
    rruns -=1
  repeat_masker_out_list.append(repeats_dir + repeat_base + ".out")
  
if config["Repeat Annotation"]["Repeat_library"]:
  repeat_masker_out_list.append(repeats_dir + repeat_base + ".out")
  if rruns == 1:
    repeat_base = os.path.basename(masked_reference).rsplit('.', 2)[0]
    rmasked[repeat_base] = rm_input_genome
    rdb[repeat_base]= " -lib " + config["Repeat Annotation"]["Repeat_library"]
  else:
    lib_name = os.path.splitext(os.path.basename(config["Repeat Annotation"]["Repeat_library"]))[0]
    repeat_base = repeat_base + ".rmask_" + lib_name
    rmasked[repeat_base] = rm_input_genome
    rdb[repeat_base]= " -lib " + config["Repeat Annotation"]["Repeat_library"]
    rm_input_genome = repeats_dir + repeat_base + ".masked.fa"
    rruns -= 1

redmask_input = ""
redmask_out_pref = ""
redmask_in_pref = ""
redmask_dir = repeats_dir + "redmask_run/"
if config["Parameters"]["run_redmask"]:
  repeat_base = os.path.basename(masked_reference).rsplit('.', 2)[0]
  redmask_input = rm_input_genome
  redmask_in_pref = os.path.basename(rm_input_genome).rsplit('.', 1)[0]
  redmask_out_pref = repeat_base
  
##RUN REPEAT ANNOTATION RULES
use rule repeat_masker from repeats_workflow with:
  input: 
   genome = lambda wildcards: rmasked[wildcards.masked]
  output:
    masked = "{dir}/{masked}.masked.fa",
    out = "{dir}/{masked}.out"
  params:
    db = lambda wildcards: rdb[wildcards.masked],
    repdir = repeats_dir,
    inbasename = lambda wildcards: os.path.basename(rmasked[wildcards.masked])
  log:
    "{dir}/logs/" + str(date) + ".j%j.repeat_masker.{masked}.out",
    "{dir}/logs/" + str(date) + ".j%j.repeat_masker.{masked}.err",
  benchmark:
    "{dir}/logs/" + str(date) + ".repeat_masker.{masked}.benchmark.txt"
  conda:
    "../envs/repeatmasker-4.1.5-0.yaml"
  threads: config["Repeat Annotation"]["rmaskCores"]

use rule redmask from repeats_workflow with:
  input:
    genome = redmask_input
  output:
    fasta = repeats_dir + redmask_in_pref + ".redmask.msk.fa",
    #fasta = repeats_dir + redmask_out_pref + ".masked.fa",
    bed = repeats_dir + redmask_in_pref + ".redmask.rep.bed"
  params:
    word_length = config["RedMask"]["wordlen"],
    min_kmer = config["RedMask"]["minkmer"],
    additional_redmask_opts = config["RedMask"]["add_option"],
    run_dir = redmask_dir
  log:
    rlogs_dir + str(date) + ".j%j.redmask.out",
    rlogs_dir + str(date) + ".j%j.redmask.err",
  benchmark:
    rlogs_dir + str(date) + ".redmask.benchmark.txt"
  conda:
    "../envs/redmask.yaml"
  threads: 1

use rule filter_prot from repeats_workflow with:
  input:
    genome = redmask_input,
    repeat_bed = repeats_dir + redmask_in_pref + ".redmask.rep.bed",
  output:
    repeat_fasta = repeats_dir + redmask_in_pref + ".redmask.rep.fasta",
    blast = repeats_dir + redmask_in_pref + ".redmask.rep.blast.out",
    masked_genome = repeats_dir + redmask_out_pref + ".masked.fa",
    redmask_bed = repeats_dir + "repeats.noblast.bed"
  params:
    blastdb = config["BLAST"]["blastdb"],
    evalue = config["BLAST"]["evalue"],
    repdir = repeats_dir
  log:
    rlogs_dir + str(date) + ".j%j.filter_prot.out",
    rlogs_dir + str(date) + ".j%j.filter_prot.err",
  benchmark:
    rlogs_dir + str(date) + ".filter_prot.benchmark.txt"
  conda:
    "../envs/blast_bedtools.yaml"
  threads: config["BLAST"]["blastCores"]

use rule get_repeats_gff from repeats_workflow with:
  input:
    redmask_bed = repeats_dir + "repeats.noblast.bed" if config["Parameters"]["run_redmask"] else [],
    repeatmasker_out = repeat_masker_out_list
  output:
    intermediate_gffs = [repeats_dir + "repeats.noblast.gff3", repeats_dir + "repeatmasker.gff3"],
    repeats_gff = config["Outputs"]["Repeat_gff"]
  params:
    repdir = "step01_Repeat_Annotation.V01/"
  log:
    rlogs_dir + str(date) + ".j%j.get_repeats_gff.out",
    rlogs_dir + str(date) + ".j%j.get_repeats_gff.err",
  benchmark:
    rlogs_dir + str(date) + ".get_repeats_gff.benchmark.txt"
  threads: 1  
